!
!****************************************************************
!
   SUBROUTINE WINDP1 (WIND10     ,THETAW     ,IDWMIN     ,      &
                      IDWMAX     ,FPM        ,UFRIC      ,      &
		      WX2        ,WY2        ,SPCDIR     ,      &
		      UX2        ,UY2        ,SPCSIG     ,      &
		      IG         )  
!
!****************************************************************
!
!        Computation of parameters derived from the wind for several
!        subroutines such as :
!        SWIND1, SWIND2, SWIND3, CUTOFF
!
!        Output of this subroutine :
!
!        WIND10 , THETAW, IDWMIN, IDWMAX , UFRIC, FPM
!
!     METHOD
!
!     a. For SWIND1 and SWIND2 :
!
!        SIGMA_FPM = 0.13 * GRAV * 2 * PI / WIND10
!
!     b. For SWIND3 (wind input according to Snyder (1981) ***
!
!       - wind friction velocity according to Wu (1982):
!           *                                                -3
!         U  =  UFRIC = wind10 sqrt( (0.8 + 0.065 wind10 ) 10  )
!
!     c. For SWIND4 (wind input according to Janssen 1991)
!
!       - wind friction velocity:
!
!          UFRIC = sqrt ( CDRAG) * U10
!
!          for U10 < 7.5 m/s  ->  CDRAG = 1.2873.e-3
!
!          else wind friction velocity according to Wu (1982):
!
!         *                                                -3
!         U  =  UFRIC = wind10 sqrt( (0.8 + 0.065 wind10 ) 10  )
!
!     d.
!        The Pierson Moskowitz radian frequency for a fully developed
!        sea state spectrum for all third generation wind input
!        models is equal to:
!
!                        grav
!        SIGMA_FPM  =  ---------
!                      28 UFRIC
!        -----------------------------------------------------------
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
   USE ALL_VARS, ONLY : MT                                                       
   
   IMPLICIT NONE
!
   REAL    :: SPCDIR(MDC,6),SPCSIG(MSC)        
   INTEGER :: IDWMIN ,IDWMAX
   INTEGER :: ID,IDDUM,IG                                         
!
   REAL    :: WIND10 ,THETAW ,UFRIC ,FPM ,CDRAG ,SDMEAN       
   REAL    :: WX2(MT), WY2(MT), UX2(MT), UY2(MT)         

   REAL         AWX, AWY, RWX, RWY                     
!
!  compute absolute wind velocity                                    
   IF(VARWI)THEN
     AWX = WX2(IG)
     AWY = WY2(IG)
   ELSE
     AWX = U10 * COS(WDIC)
     AWY = U10 * SIN(WDIC)
   END IF
!  compute relative wind velocity                                      
   IF(ICUR == 0)THEN
     RWX = AWX
     RWY = AWY
   ELSE
     RWX = AWX - UX2(IG)
     RWY = AWY - UY2(IG)
   END IF
!  compute absolute value of relative wind velocity                   
   WIND10 = SQRT(RWX**2+RWY**2)
   IF(WIND10 > 0.)THEN
     THETAW = ATAN2 (RWY,RWX)
     THETAW = MOD ( (THETAW + PI2_W) , PI2_W )
   ELSE
     THETAW = 0.
   END IF
!
!  *** compute the minimum and maximum counter for the active  ***
!  *** wind field :                                            ***
!  ***                                   .                     ***
!  *** IDWMAX =135 degrees             .<--mean wind direction ***
!  ***              \ o o o | o o o  .     (THETAW)            ***
!  ***                \ o o | o o  .  o                        ***
!  ***                  \ o | o  .  o o                        ***
!  ***                    \ |  .  o o o                        ***
!  ***      ----------------\--------------------              ***
!  ***                      | \ o o o o                        ***
!  ***                      |   \ o o o                        ***
!  ***                      |     \ o o                        ***
!  ***                             IDWMIN = 325 degrees        ***
!  ***                                                         ***
!
!  move ThetaW to the right interval, shifting + or - 2*PI             
   SDMEAN = 0.5 * (SPCDIR(1,1) + SPCDIR(MDC,1))                   
   IF(THETAW < SDMEAN - PI_W) THETAW = THETAW + 2.*PI_W
   IF(THETAW > SDMEAN + PI_W) THETAW = THETAW - 2.*PI_W
!
   IF((THETAW-0.5*PI_W) <= SPCDIR(1,1))THEN                   
     IF((THETAW+1.5*PI_W) >= SPCDIR(MDC,1))THEN
       IDWMIN = 1
     ELSE
       IDWMIN = NINT ( (THETAW + 1.5*PI_W - SPCDIR(1,1)) / DDIR ) + 1    
     END IF
   ELSE
     IDWMIN = NINT ( (THETAW - 0.5*PI_W - SPCDIR(1,1)) / DDIR ) + 1      
   END IF
!
   IF((THETAW + 0.5 * PI_W) >= SPCDIR(MDC,1))THEN                  
     IF((THETAW - 1.5 * PI_W) <= SPCDIR(1,1))THEN                  
       IDWMAX = MDC
     ELSE
       IDWMAX = NINT ( (THETAW - 1.5 * PI_W - SPCDIR(1,1)) / DDIR ) + 1  
     END IF
   ELSE
     IDWMAX = NINT ( (THETAW + 0.5 * PI_W - SPCDIR(1,1)) / DDIR ) + 1    
   END IF
!
   IF(IDWMIN > IDWMAX) IDWMAX = MDC + IDWMAX
!
!  *** compute the Pierson Moskowitz frequency ***
!
   IF(IWIND == 1 .OR. IWIND == 2)THEN
!
!    *** first and second generation wind wave model ***
!
     IF(WIND10 < PWIND(12) ) WIND10 = PWIND(12)
     FPM = 2. * PI_W * PWIND(13) * GRAV_W / WIND10
!
!    *** determine U friction in case predictor is obtained ***
!    *** with second genaration wave model                  ***
!
     IF(WIND10 > 7.5)THEN
       UFRIC = WIND10 * SQRT (( 0.8 + 0.065 * WIND10 ) * 0.001 )
     ELSE
       CDRAG = 0.0012873
       UFRIC = SQRT ( CDRAG ) * WIND10
     END IF
!
!    Reformulation of the wind speed in terms of friction velocity.      
!    This formulation is based on Bouws (1986) and described in Delft    
!    Hydraulics report H3515 (1999)                                      
!
     WIND10 = WIND10 * SQRT(((0.8 + 0.065 * WIND10) * 0.001) /    &
                            ((0.8 + 0.065 * 15.   ) * 0.001))          
!
   ELSE IF(IWIND >= 3)THEN
!
!    *** Calculate the wind friction velocity  ***
!    *** according to Wu (1982)                ***
!
     IF(WIND10 > 7.5)THEN
       UFRIC = WIND10 * SQRT (( 0.8 + 0.065 * WIND10 ) * 0.001 )
     ELSE
       CDRAG = 0.0012873
       UFRIC = SQRT ( CDRAG ) * WIND10
     END IF
!
!    *** Wind friction velocity and PM-frequency ***
!
     UFRIC = MAX ( 1.E-15 , UFRIC)
     FPM =  GRAV_W / ( 28.0 * UFRIC )

   END IF

   RETURN
   END SUBROUTINE WINDP1
 
!
!****************************************************************
!
   SUBROUTINE WINDP2 (IDWMIN  ,IDWMAX  ,SIGPKD  ,FPM     ,      &
                      ETOTW   ,                                 &
		      AC2     ,SPCSIG  ,                        &
		      WIND10                                      )    
!  (This subroutine has not been tested yet)
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
   USE ALL_VARS, ONLY : MT                                                        
!
!
!  2. Purpose
!
!     Computation of wind sea energy spectrum for the second
!     generation wind growth model. Output of the subroutine
!     is ETOTW (total wind sea energy spectrum)
!
!  3. Method
!
!     Compute the wind sea energy spectrum : ETOTW
!
!             +90  inf
!             |   |
!     ETOTW = |   |        E(s,d) ds dd
!            -90  0.7*FPM
!
!     Compute the total energy density ETOTW for F > 0.7 FPM
!
!          ^  |
!       E()|  |          *
!             |        *   *
!             |              *
!             |       *      | *      / ETOTW(e)
!             |              | o  * /
!             |      *       | o o/o *
!             |              | o o o o o *
!             |     *        | o o o o o o o o*
!            0---------------|---------------------------
!             0          0.7*FPM              SIGMA --> s
!
!                   SIGMA MAX
!                  |
!      ETOTW(d) =  |  E(s,d)ds      ISFPM = FPM
!                 0.7 FPM
!
!      and over the interval +/- 90 degrees (according to
!      the mean wind direction as computed in WINDP1 )
!
!            IDWMAX = "135"
!               o
!      ETOTW =  o  ETOTW(d)dd
!
!         IDWMIN ="325"
!
!
   REAL    SPCSIG(MSC)                                                 
!
!     9. STRUCTURE
!
!     ------------------------------------------------------------
!     Determine the counter for the FPM frequency
!     integrate the wind sea energy spectrum
!     add the energy of the high frequency tail to the spectrum
!     ------------------------------------------------------------
!
!***********************************************************************
!
!
   INTEGER  IDWMIN  ,IDWMAX  ,IDDUM   ,ID      ,IS      ,ISFPM
!
   REAL     ETOTW   ,FPM     ,SIG     ,ATOTD

   REAL     AC2(MDC,MSC,0:MT)
!
!  *** compute wind sea energy spectrum for IS > 0.7 FPM       ***
!  *** minimum FPM is equal : 2 * pi * 0.13 * grav / pwind(12) ***
!  *** is equal 8 rad/s = 1.27 Hz                              ***
!
   ISFPM = MSC                                                         
   FACINT = 0.                                                         
   DO IS = 1, MSC
     SIG = SPCSIG(IS)                                                  
     IF(FRINTH * SIG > (0.7 * FPM))THEN
       ISFPM =  IS
       FACINT = (FRINTH - 0.7*FPM/SIG) / (FRINTH - 1./FRINTH)
       GOTO 11
     END IF
   ENDDO
11 CONTINUE
!
!  *** calculate the energy in the wind sea part of the spectrum ***
!  *** from ISFPM.                                               ***
!
   ETOTW = 0.
   DO IS = ISFPM, MSC
     SIG = SPCSIG(IS)                                                  
     ATOTD = 0.
     DO IDDUM = IDWMIN, IDWMAX
       ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
       ATOTD = ATOTD + AC2(ID,IS,KCGRD(1))                             
     ENDDO
     IF(IS == ISFPM)THEN
       ETOTW = ETOTW + FACINT * FRINTF * SIG**2 * DDIR * ATOTD
     ELSE
       ETOTW = ETOTW + FRINTF * SIG**2 * DDIR * ATOTD
     ENDIF
   ENDDO
!  add high-frequency tail:
   ETOTW = ETOTW + PWTAIL(6) * SIG**2 * DDIR * ATOTD

   RETURN
   END SUBROUTINE WINDP2
!
!********************************************************************
!
      SUBROUTINE WINDP3 (ISSTOP  ,ALIMW   ,AC2     ,               &
                         GROWW   ,IDCMIN  ,IDCMAX  )                 
!    (This subroutine has not been tested yet)
!
!****************************************************************
!
      USE SWCOMM3                                                        
      USE SWCOMM4                                                        
      USE OCPCOMM4  
      USE ALL_VARS, ONLY : MT                                                     
!
!
!     2. PURPOSE
!
!        Reduce the energy density in spectral direction direct after
!        solving the tri-diagonal matrix if the energy density level is
!        larger then the upper bound limit given by a Pierson Moskowitz
!        spectrum. This is only carried out if a particular wave
!        component is 'growing'.
!
!        If the energy density in a bin is larger than the upper bound
!        limit (for instance when crossing wind seas are present) then
!        the energy density level is a lower limit.
!
!     3. METHOD
!
!        The upper bound limit is given by:
!
!                       2                  -4
!               ALPHA  g              sigma     2     2
!    A(s,t) = ----------- exp ( -5/4 (----)   ) --- cos ( d - dw )
!               sigma^6                FPM       pi
!
!         in which ALPHA is wind sea dependent (see subroutine
!         SWIND2) :
!
!                      /                                 C60   \
!                     |                / E_windsea g^2 \        |
!         ALPHA = MIN | 0.0081 ,  C50 |  ------------   |       |
!                      \               \    U_10^4     /        /
!
!
!     4. PARAMETERLIST
!
!        INTEGERS :
!        ----------
!        IX IY            Grid point in geographical space
!        MDC ,MSC         Counters in spectral space
!        ISSTOP           Maximum frequency that fall within a sweep
!
!        ARRAYS:
!        -------
!        AC2      4D      Action density as funciton of x,y,s,t
!        ALIMW    2D      Contains the action density upper bound
!                         limit regardind spectral action density per
!                         spectral bin (A(s,t))
!        GROWW    2D      Logical array which determines is there is
!                         a) generation  ( E < E_lim -> .TRUE.  ) or
!                         b) dissipation ( E > E_lim -> .FALSE. )
!        IDCMIN   1D      Frequency dependent minimum counter
!        IDCMAX   1D      Frequency dependent minimum counter
!
!     5. SUBROUTINES CALLING
!
!        SWOMPU
!
!     6. SUBROUTINES USED
!
!        NONE
!
!     7. Common blocks used
!
!
!     8. REMARKS
!
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!   For every spectral bin
!     If wave input for a bin is true (GROWW(ID,IS) = .TRUE.) and
!        if wave action is larger then maximum then reduce
!        wave action to limit its value to upper boundary limit
!     else if no growth (see subroutine SWIND1 and SWIND2), check
!       if the energy density level in the wind sea part
!       of the spectrum is lower than upper bound limit. If so
!       set energy level equal to upper bound limit
!     endif
!   ------------------------------------------------------------
!   End of the subroutine WINDP3
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!***********************************************************************
!
      INTEGER     IS      ,ID      ,ISSTOP  ,IDDUM
!
      INTEGER     IDCMIN(MSC)       ,IDCMAX(MSC)
!
      REAL        AC2CEN
!
      REAL        AC2(MDC,MSC,0:MT)    ,ALIMW(MDC,MSC)
!
      LOGICAL     GROWW(MDC,MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'WINDP3')
!
!     *** limit the action density spectrum ***
!
      DO IS = 1, ISSTOP
        DO IDDUM = IDCMIN(IS) , IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
          AC2CEN = AC2(ID,IS,KCGRD(1))
          IF ( GROWW(ID,IS) .AND. AC2CEN .GT. ALIMW(ID,IS) )            &
	       AC2(ID,IS,KCGRD(1)) = ALIMW(ID,IS)
          IF ( .NOT. GROWW(ID,IS) .AND. AC2CEN .LT. ALIMW(ID,IS) )      &
	       AC2(ID,IS,KCGRD(1)) = ALIMW(ID,IS)
!
          IF (TESTFL .AND. ITEST .GE. 50) THEN
             WRITE(PRINTF,300) IS,ID,GROWW(ID,IS),AC2CEN,ALIMW(ID,IS)
 300         FORMAT(' WINDP3 : IS ID GROWW AC2CEN ALIM:',2I4,L4,2E12.4)
          END IF
!
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (TESTFL .AND. ITEST .GE. 50) THEN
        WRITE(PRINTF,4000) KCGRD(1),ISSTOP,MSC,MDC,MCGRD
 4000   FORMAT(' WINDP3 : POINT ISSTOP MSC MDC MCGRD :',5I5)
      END IF
!
      RETURN
      END SUBROUTINE WINDP3
!
!****************************************************************
!
   SUBROUTINE SWIND0 (IDCMIN  ,IDCMAX  ,ISSTOP  ,SPCSIG  ,THETAW  , &
		      UFRIC   ,FPM     ,PLWNDA  ,IMATRA  ,SPCDIR  )
!
!****************************************************************
!
!     Computation of the source term for the wind input for a
!     third generation wind growth model:
!
!     1)  Linear wind input term according to Cavaleri and
!         Malanotte-Rizzoli (1981)
!
!  METHOD
!
!     To ensure wave growth when no wave energy is present in the
!     numerical model a linear growth term is used (see
!     Cavaleri and Malanotte-Rizzoli 1981). Contributions for
!     frequencies lower then FPM have been eliminated buy aa filter :
!
!                  -3
!            1.5*10      *                       4      sigma -4
!  A = Sin = ------- { U  max[0, (cos(d - dw )] } exp{-(-----)  } / Jac
!               2                                        FPM
!              g
!
!        With Jac = Jacobian =  2 pi sigma
!
!        The Pierson Moskowitz radian frequency for a fully developed
!        sea state spectrum is as follows (computed in WINDP2 :
!
!                1       g
!        FPM  = ---- ---------  * 2 pi
!               2 pi  28 UFRIC
!
!****************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   
   IMPLICIT NONE

   REAL    :: SPCDIR(MDC,6) ,SPCSIG(MSC)                                                 
   INTEGER :: IDDUM   ,ID      ,IS      ,ISSTOP
   REAL    :: FPM     ,UFRIC   ,THETA   ,THETAW  ,SWINEA  ,SIGMA   ,  &
              TEMP1   ,TEMP2   ,CTW     ,STW     ,COSDIF  ,TEMP3   ,  &
	      FILTER  ,ARGU
   REAL    :: IMATRA(MDC,MSC)      ,PLWNDA(MDC,MSC,NPTST)        
   INTEGER :: IDCMIN(MSC)          ,IDCMAX(MSC)
!
!     *** calculate linear wind input term ***
!
   CTW = COS(THETAW)                                                   
   STW = SIN(THETAW)                                                   
   FPM =  GRAV_W / ( 28.0 * UFRIC )
   TEMP1 = PWIND(31) / ( GRAV_W**2 * 2. * PI_W ) 
   DO IS = 1, ISSTOP
     SIGMA  = SPCSIG(IS)                                               
!
!    ****            ARGU   =  FPM / SIGMA                     ***
!    **** the value of ARGU was change for MIN () because for  ***
!    **** values of fpm/sigma too small could be some problems ***
!    **** with some computers to handle small numbers          ***
     ARGU   = MIN (2., FPM / SIGMA)                                    
     FILTER = EXP ( - ARGU**4 )                                        
     TEMP2  = TEMP1 / SIGMA
     DO IDDUM = IDCMIN(IS), IDCMAX(IS)
       ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
       IF(SIGMA >= (0.7 * FPM))THEN
         THETA  = SPCDIR(ID,1)                                         
         COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  
         TEMP3  = ( UFRIC *  MAX( 0. , COSDIF))**4                     
         SWINEA = MAX( 0. , TEMP2 * TEMP3 * FILTER )
         IMATRA(ID,IS) = IMATRA(ID,IS) + SWINEA
         IF(TESTFL) PLWNDA(ID,IS,IPTST) = SWINEA                       
       END IF
     ENDDO
   ENDDO

   RETURN
   END SUBROUTINE SWIND0

!
!****************************************************************
!
   SUBROUTINE SWIND3 (SPCSIG  ,THETAW  ,IMATDA  ,KWAVE   ,IMATRA  , &
		      IDCMIN  ,IDCMAX  ,UFRIC   ,FPM     ,PLWNDB  , &
		      ISSTOP  ,SPCDIR  ,IG      )
!
!****************************************************************
!
!     Computation of the source term for the wind input for a
!     third generation wind growth model:
!
!     1)  Exponential input term, (Snyder et al. 1981, which
!         expression has been modified by Komen et al. 1984).
!
!         This input term should be combinated with the dissipation
!         term of Komen et al. (1984)
!
!  METHOD
!
!     The exponential term used is taken from Snyder et al. (1981)
!     and Komen et al. (1984):
!
!     Sin (s,d) =  B*E(s,d)
!        e                             *
!                                 28 U cos( d - dw )
!     B = max(0. , (0.25 rhoaw( -------------------  -1 ) )) sigma
!                                   sigma / kwave
!
!     with :
!
!        *                                                -3
!      U  =  UFRIC = wind10 sqrt( (0.8 + 0.065 wind10 ) 10  )
!
!     UFRIC is computed in WINDP1
!
!
!     The Pierson Moskowitz radian frequency for a fully developed
!     sea state spectrum is as follows (computed in WINDP1):
!
!             1       g
!     FPM  = ---- ---------  * 2 pi
!            2 pi  28 UFRIC
!
!
!****************************************************************
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4   
!   USE ALL_VARS, ONLY : MT,AC2  
   USE VARS_WAVE, ONLY : MT,AC2  
   
   IMPLICIT NONE                                                   

   REAL    :: SPCDIR(MDC,6),SPCSIG(MSC)                                                 
   INTEGER :: IDDUM ,ID    ,IS    ,ISSTOP  ,IG
   REAL    :: FPM   ,UFRIC ,THETA ,THETAW,SIGMA ,SWINEB,TEMP1,          &
              CTW   ,STW   ,COSDIF,TEMP2 ,TEMP3 ,CINV
   REAL    :: IMATDA(MDC,MSC),IMATRA(MDC,MSC),        &
              KWAVE(MSC,ICMAX),PLWNDB(MDC,MSC,NPTST)
   INTEGER :: IDCMIN(MSC),IDCMAX(MSC)

   CTW   = COS(THETAW)                                                 
   STW   = SIN(THETAW)                                                 
   TEMP1 = 0.25 * PWIND(9)
   TEMP2 = 28.0 * UFRIC
   DO IS = 1, ISSTOP
     SIGMA = SPCSIG(IS)                                                
     CINV  = KWAVE(IS,1) / SIGMA
     TEMP3 = TEMP2 * CINV
     DO IDDUM = IDCMIN(IS), IDCMAX(IS)
       ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
       THETA  = SPCDIR(ID,1)                                         
       COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  
       SWINEB = TEMP1 * ( TEMP3 * COSDIF - 1.0 )   
       SWINEB = MAX ( 0. , SWINEB * SIGMA )

       IMATRA(ID,IS) = IMATRA(ID,IS) + SWINEB * AC2(ID,IS,IG)
       IMATDA(ID,IS) = IMATDA(ID,IS) - SWINEB
       IF(TESTFL) PLWNDB(ID,IS,IPTST) = SWINEB                      
     ENDDO
   ENDDO

   RETURN
   END SUBROUTINE SWIND3
 
!
!****************************************************************
!
   SUBROUTINE SWIND4 (IDWMIN  ,IDWMAX  ,SPCSIG  ,WIND10  ,THETAW  ,  &
                      XIS     ,DD      ,KWAVE   ,IMATRA  ,IMATDA  ,  &
		      IDCMIN  ,IDCMAX  ,UFRIC   ,PLWNDB  ,ISSTOP  ,  &
		      ITER    ,USTAR   ,ZELEN   ,SPCDIR  ,IT      ,  &
		      IG      )
!
!******************************************************************
!
!     Computation of the source term for the wind input for a
!     third generation wind growth model:
!
!     1)  Computation of the exponential input term based on a
!         quasi linear theory developped by Janssen (1989, 1991).
!         This formulation should be used in combination with the
!         whitecapping dissipation source term according to
!         Janssen (1991) and Mastenbroek et al. (1993)
!
!  Method
!
!     The exponential term for the wind input used is taken
!     from Janssen (1991):
!
!     Sin (s,d) =  B*E(s,d)
!        e
!                                        *
!                               / kwave U  \    2
!     B = max(0. , (beta rhoaw | --------- | cos( d - dw ) sigma))
!                               \ sigma   /
!
!     The friction velocity is a fucntion of the roughness
!     length Ze. A first gues for U* is given by Wu (1982),
!     which is computed in subroutine WINDP2.
!
!******************************************************************
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4  
!   USE ALL_VARS, ONLY : MT,AC2    
   USE VARS_WAVE, ONLY : MT,AC2    
   
   IMPLICIT NONE                                                  
   REAL    :: SPCDIR(MDC,6),SPCSIG(MSC)                                                 
   INTEGER :: IDWMAX  ,IDWMIN  ,IDDUM   ,ID      ,ISSTOP  ,IS
   REAL    :: THETA  ,THETAW ,DD     ,SWINEB ,WIND10 ,                   &
              ZO     ,ZE     ,BETA1  ,BETA2  ,UFRIC  ,UFRIC2 ,DS     ,   &
	      ZARG   ,ZLOG1  ,ZLOG2  ,ZCN1   ,ZCN2   ,ZCN    ,XIS    ,   &
	      SIGMA  ,SIGMA1 ,SIGMA2 ,WAVEN  ,WAVEN1 ,WAVEN2 ,TAUW   ,   &
	      TAUTOT ,TAUDIR ,COS1   ,COS2   ,CW1    ,RHOA   ,RHOW   ,   &
	      RHOAW  ,ALPHA  ,XKAPPA ,F1     ,TAUWX  ,TAUWY  ,SE1    ,   &
	      CTW    ,STW    ,COSDIF ,                                   &
	      SE2    ,SINWAV ,COSWAV
   REAL    :: IMATDA(MDC,MSC)      ,                                     &
              IMATRA(MDC,MSC)      ,                                      &
	      KWAVE(MSC,ICMAX)     ,                                      &
	      PLWNDB(MDC,MSC,NPTST),                                      &
	      USTAR(MT)         ,                                      &
	      ZELEN(MT)
   INTEGER :: IDCMIN(MSC)          ,IDCMAX(MSC)
   REAL    :: ZTEN ,RATIO ,BETAMX ,TXHFR ,TYHFR ,CW2 ,ZARG1 ,ZARG2
   REAL    :: GAMHF ,SIGMAX ,SIGHF1 ,SIGHF2 ,ZCNHF1 ,ZCNHF2 ,AUX
   REAL    :: ZAHF1 ,ZAHF2 ,FACHFR ,FA ,FB ,FC ,FD ,FE ,FCEN
   REAL    :: FF1 ,FF2 ,FF3 ,DCEN ,TAUNEW ,XFAC2 ,BETA ,ZLOG
   INTEGER :: IT ,IG ,ITER ,J ,II                                         

!
!  *** initialization ***
!
   ALPHA  = PWIND(14)
   XKAPPA = PWIND(15)
   RHOA   = PWIND(16)
   RHOW   = PWIND(17)
   RHOAW  = RHOA / RHOW
   ZTEN   = 10.
   RATIO  = 0.75
   BETAMX = 1.2
   F1     = BETAMX / XKAPPA**2
   CTW    = COS(THETAW)                                                
   STW    = SIN(THETAW)                                                
!
   IF(NSTATC == 1 .AND. IT == 1)THEN                               
!
!     *** nonstationary and first time step (the number of        ***
!     *** iterations however still can increase per time step     ***
!
     ZO     = ALPHA * UFRIC * UFRIC / GRAV_W
     ZE     = ZO / SQRT( 1. - RATIO )
     USTAR(IG) = UFRIC                                           
     ZELEN(IG) = ZE                                              
   ELSE IF(NSTATC == 0 .AND. ICOND == 4 .AND. ITER == 1)THEN     
!
!    *** non-first stationary computations and first iteration   ***
!
     ZO     = ALPHA * UFRIC * UFRIC / GRAV_W
     ZE     = ZO / SQRT( 1. - RATIO )
     USTAR(IG) = UFRIC                                           
     ZELEN(IG) = ZE                                              
   ELSE IF ( NSTATC.EQ.0 .AND. ICOND.NE.4 .AND. ITER .EQ. 2 ) THEN     
!
!    *** first stationary computation (this subroutine is never ***
!    *** excecuted anyway, this subroutine in entered after 1   ***
!    *** iteration) and thus calculate ZO and ZE as a first     ***
!    *** prediction only and only in the second sweep           ***
!
     ZO     = ALPHA * UFRIC * UFRIC / GRAV_W
     ZE     = ZO / SQRT( 1. - RATIO )
     USTAR(IG) = UFRIC
     ZELEN(IG) = ZE
   ELSE
!
!    *** calculate wave stress using the value of the  ***
!    *** velocity U* and roughness length Ze from the  ***
!    *** previous iteration                            ***
!
     UFRIC = USTAR(IG)
     ZE    = ZELEN(IG)
!
     TAUW   = 0.
     TAUWX  = 0.
     TAUWY  = 0.
     TXHFR  = 0.
     TYHFR  = 0.
!
!    *** use old friction velocity to calculate wave stress ***
!
     UFRIC2 = UFRIC * UFRIC
!
     DO IS = 1, MSC-1
       SIGMA1 = SPCSIG(IS)                                             
       SIGMA2 = SPCSIG(IS+1)                                           
       WAVEN1 = KWAVE(IS,1)
       WAVEN2 = KWAVE(IS+1,1)
       DS     = SIGMA2 - SIGMA1
       CW1    = SIGMA1 / WAVEN1
       CW2    = SIGMA2 / WAVEN2
       ZCN1   = ALOG ( GRAV_W * ZE / CW1**2 )
       ZCN2   = ALOG ( GRAV_W * ZE / CW2**2 )
       DO IDDUM = IDWMIN, IDWMAX
         ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
         THETA  = SPCDIR(ID,1)                                         
         COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  
         SINWAV = SPCDIR(ID,3)                                         
         COSWAV = SPCDIR(ID,2)                                         
         COS1   = MAX ( 0. , COSDIF )                                  
         COS2   = COS1 * COS1
         BETA1  = 0.
         BETA2  = 0.
!
!        *** Miles constant Beta ***
!
         IF(COS1 > 0.01)THEN
           ZARG1 = XKAPPA * CW1 / ( UFRIC * COS1 )
           ZARG2 = XKAPPA * CW2 / ( UFRIC * COS1 )
           ZLOG1 = ZCN1 + ZARG1
           ZLOG2 = ZCN2 + ZARG2
           IF(ZLOG1 < 0.) BETA1 = F1 * EXP (ZLOG1) * ZLOG1**4
           IF(ZLOG2 < 0.) BETA2 = F1 * EXP (ZLOG2) * ZLOG2**4
         ENDIF
!
!        *** calculate wave stress by integrating input source ***
!        *** term in x- and y direction respectively           ***
!
         SE1 = WAVEN1**2 * BETA1 * SIGMA1 * AC2(ID,IS  ,IG)
         SE2 = WAVEN2**2 * BETA2 * SIGMA2 * AC2(ID,IS+1,IG)
!
         TAUWX = TAUWX + 0.5 * ( SE1 + SE2 ) * DS * COSWAV * COS2
         TAUWY = TAUWY + 0.5 * ( SE1 + SE2 ) * DS * SINWAV * COS2
       ENDDO
     ENDDO
!
!    *** determine effect of high frequency tail to wave stress ***
!    *** assuming deep water conditions                         ***
!
     GAMHF =  XKAPPA * GRAV_W / UFRIC
     SIGMAX = SPCSIG(MSC)                                              
     SIGHF1 = SIGMAX
     DO J=1, 50
       SIGHF2 = XIS * SIGHF1
       DS     = SIGHF2 - SIGHF1
       ZCNHF1 = ALOG ( ZE * SIGHF1**2 / GRAV_W )
       ZCNHF2 = ALOG ( ZE * SIGHF2**2 / GRAV_W )
       AUX    = 0.0
       DO IDDUM = IDWMIN, IDWMAX
         ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
         THETA  = SPCDIR(ID,1)                                         
         COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  
         SINWAV = SPCDIR(ID,3)                                         
         COSWAV = SPCDIR(ID,2)                                         
         COS1   = MAX ( 0. , COSDIF )                                  
         COS2   = COS1 * COS1
         BETA1  = 0.0
         BETA2  = 0.0
!
         IF(COS1 > 0.01)THEN
!          *** beta is independent of direction ! ***
           ZAHF1 = GAMHF / SIGHF1
           ZAHF2 = GAMHF / SIGHF2
           ZLOG1 = ZCNHF1 + ZAHF1
           ZLOG2 = ZCNHF2 + ZAHF2
           IF(ZLOG1 < 0.) BETA1 = F1 * EXP (ZLOG1) * ZLOG1**4
           IF(ZLOG2 < 0.) BETA2 = F1 * EXP (ZLOG2) * ZLOG2**4
           AUX = AUX + BETA1 + BETA2
         ENDIF
!
!        *** calculate contribution of high frequency tail to ***
!        *** wave stress by integrating input source term in  ***
!        *** x- and y direction respectively                  ***
!
         FACHFR = SIGMAX**6 * AC2(ID,MSC,IG) * COS2 / GRAV_W**2    

         SE1 = FACHFR * BETA1 / SIGHF1
         SE2 = FACHFR * BETA2 / SIGHF2
!
         TXHFR = TXHFR + 0.5 * ( SE1 + SE2 ) * DS * COSWAV
         TYHFR = TYHFR + 0.5 * ( SE1 + SE2 ) * DS * SINWAV
!
!        *** if coeffcient BETA = 0. for a frequency over ***
!        *** all directions is zero skip loop             ***
!
         IF(AUX == 0.) GOTO 5000
!
       ENDDO

       SIGHF1 = SIGHF2
     ENDDO
5000 CONTINUE

     TAUTOT = RHOA * UFRIC2
!    *** wave stress ***
     TAUWX  = TAUWX + TXHFR
     TAUWY  = TAUWY + TYHFR
     IF(ABS(TAUWX) > 0. .OR. ABS(TAUWY) > 0.)THEN
       TAUDIR = ATAN2 ( TAUWX, TAUWY )
     ELSE
       TAUDIR = 0.
     ENDIF
     TAUDIR = MOD ( (TAUDIR + 2. * PI_W) , (2. * PI_W) )
     TAUW   = RHOA * UFRIC2 * DD * SQRT ( TAUWX**2 + TAUWY**2 )
     TAUW   = MIN ( TAUW , 0.999 * TAUTOT )

     DO II = 1, 20
!      *** start iteration process ***
       FA = SQRT ( 1. - TAUW / TAUTOT )
       FB = ZTEN * RHOA * GRAV_W / ALPHA
       FC = FA * ( FB / TAUTOT  - 1. )
       FD = SQRT ( TAUTOT )
       FE = ALOG ( FC + 1. )
!
!      *** calculate function value and derivative in ***
!      *** numerical point considered                 ***
!
       FCEN = FD * FE - SQRT(RHOA) * WIND10 * XKAPPA
       FF1  = 0.5 * FE / FD
       FF2  = 0.5 * TAUW * FC / FA**2 - FA * FB
       FF3  = TAUTOT**1.5 * ( FC + 1. )
       DCEN = FF1 + FF2 / FF3
!
!      *** new total stress ***
!
       TAUNEW = TAUTOT - FCEN / DCEN
!
       IF(TAUNEW <= TAUW) TAUNEW = .5 * (TAUTOT + TAUW)           
       IF(ABS( TAUNEW - TAUTOT ) <= 1.E-5 ) GOTO 3000
!
       TAUTOT = TAUNEW
     ENDDO
3000 CONTINUE
!
     UFRIC  = SQRT ( TAUTOT / RHOA )
!
     ZO     = ALPHA * UFRIC * UFRIC / GRAV_W
     ZE     = ZO / SQRT ( 1. - TAUW / TAUTOT )
!
     USTAR(IG) = UFRIC
     ZELEN(IG) = ZE
!
   ENDIF
!
! ----->
!
!     *** calculate critical height and Miles parameter and  ***
!     *** calculate input source term B for with the updated ***
!     *** values of UFRIC and ZE                             ***
!
   UFRIC2 = UFRIC * UFRIC
!
   DO IS = 1, ISSTOP
     SIGMA  = SPCSIG(IS)                                               
     WAVEN  = KWAVE(IS,1)
     CW1    = SIGMA / WAVEN
     ZCN    = ALOG ( GRAV_W * ZE / CW1**2 )
     DO IDDUM = IDCMIN(IS), IDCMAX(IS)
       ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
         THETA  = SPCDIR(ID,1)                                         
         COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  
         COS1   = MAX ( 0. , COSDIF )                                  
         COS2   = COS1 * COS1
         XFAC2  = ( UFRIC / CW1 )**2
         BETA   = 0.
         IF(COS1 > 0.01)THEN
           ZARG = XKAPPA * CW1 / ( UFRIC * COS1 )
           ZLOG = ZCN + ZARG
           IF(ZLOG < 0.) BETA = F1 * EXP (ZLOG) * ZLOG**4
         ENDIF
!
!        *** compute the factor B and store result in array ***
!
         SWINEB = RHOAW * BETA * XFAC2 * COS2 * SIGMA
         IMATRA(ID,IS) = IMATRA(ID,IS) + SWINEB * AC2(ID,IS,IG)  
         IMATDA(ID,IS) = IMATDA(ID,IS) - SWINEB  
         IF(TESTFL) PLWNDB(ID,IS,IPTST) = SWINEB                      
     ENDDO
   ENDDO

   RETURN
   END SUBROUTINE SWIND4
 
!
!****************************************************************
!
   SUBROUTINE SWIND5 (SPCSIG  ,THETAW  ,ISSTOP  ,UFRIC   ,KWAVE   , &
                      IMATRA  ,IMATDA  ,IDCMIN  ,IDCMAX  ,PLWNDB  , &
		      SPCDIR  ,IG      )
!
!****************************************************************
!
!     Computation of the source term for the wind input for a
!     third generation wind growth model:
!
!     The exponential input term is according to Yan (1987). This
!     input term is valid for the higher frequency part of the
!     spectrum (strongly forced wave components). The expression
!     reduces to the Snyder (1982) expression form for spectral
!     wave components with weak wind forcing and to the Plant (1982)
!     form for more strongly forced wave components:
!
!  Method
!
!     The expression reads -->   with  X = Ustar / C
!
!            / /      2                      \
!     Sin = | | 0.04 X + 0.00544 X + 0.000055 | * cos (theta)
!            \ \                             /
!                       \
!              - 0.00031 | sigma * AC2(d,s,x,y)
!                       /
!****************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4 
!   USE ALL_VARS, ONLY : MT,AC2   
   USE VARS_WAVE, ONLY : MT,AC2   
   
   IMPLICIT NONE                                                    

   REAL    :: SPCDIR(MDC,6),SPCSIG(MSC)                                                 
   INTEGER :: IDDUM  ,ID     ,IS     ,ISSTOP,  IG
   REAL    :: UFRIC  ,THETA  ,THETAW ,SIGMA  ,SWINEB ,                &
              CTW    ,STW    ,COSDIF ,                                &
	      USTAC1 ,USTAC2 ,COF1   ,COF2   ,COF3   ,COF4
   REAL    :: IMATDA(MDC,MSC)   ,IMATRA(MDC,MSC)      , &
              KWAVE(MSC,ICMAX)  ,PLWNDB(MDC,MSC,NPTST)
   INTEGER :: IDCMIN(MSC)       ,IDCMAX(MSC)
   REAL    :: TEMP3

!
!  *** input according to Yan (1987) ***
!
   COF1 = 0.04
   COF2 = 0.00544
   COF3 = 0.000055
   COF4 = 0.00031
!
!  adapted Yan fit for use of Alves and Banner method                  
!
   IF(IWCAP == 7)THEN                                                
     COF1 = 0.04                                                      
     COF2 = 0.00552                                                   
     COF3 = 0.000052                                                  
     COF4 = 0.000302                                                  
   END IF                                                              
!
   CTW  = COS(THETAW)                                                  
   STW  = SIN(THETAW)                                                  
   DO IS = 1, ISSTOP
     SIGMA  = SPCSIG(IS)
     USTAC1 = ( UFRIC * KWAVE(IS,1) ) / SIGMA
     USTAC2 = USTAC1 * USTAC1
     TEMP3  = ( COF1 * USTAC2 + COF2 * USTAC1 + COF3)
     DO IDDUM = IDCMIN(IS), IDCMAX(IS)
       ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
         THETA  = SPCDIR(ID,1)                                         
         COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  
         SWINEB = TEMP3 * COSDIF - COF4                                
         SWINEB = MAX ( 0. , SWINEB * SIGMA )
!         IMATRA(ID,IS) = IMATRA(ID,IS) + SWINEB * AC2(ID,IS,IGC)
         IMATRA(ID,IS) = IMATRA(ID,IS) + SWINEB * AC2(ID,IS,IG)
         IMATDA(ID,IS) = IMATDA(ID,IS) - SWINEB
         IF(TESTFL) PLWNDB(ID,IS,IPTST) = SWINEB                      
     ENDDO
   ENDDO

   RETURN
   END SUBROUTINE SWIND5
!
 
!
!****************************************************************
!
   SUBROUTINE WNDPAR (ISSTOP,IDWMIN,IDWMAX,IDCMIN,IDCMAX,       &
                      DEP2  ,WIND10,THETAW,AC2   ,KWAVE ,       &
		      IMATRA,IMATDA,SPCSIG,CGO   ,ALIMW ,       &
		      GROWW ,ETOTW ,PLWNDA,PLWNDB,SPCDIR,       &
		      ITER  )   
!  (This subroutine has not been tested yet)		      
!
!****************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4 
   USE ALL_VARS, ONLY : MT                                                       
!
   IMPLICIT NONE                                                       
!
!
!  2. Purpose
!
!     Computation of the wind input source term with formulations
!     of a first-generation model (constant porportionality coefficient)
!     and a second-generation model (proportionality coefficient depends
!     on the energy in the wind sea part of the spectrum). The
!     expressions are from Holthuijsen and de Boer (1988) and from
!     the DOLPHIN-B model (Holthuijsen and Booij). During the
!     implementation of the terms modifications to the code have been
!     made after personal communications with Holthuijsen and Booij.
!
!  3. Method
!
!     The source term of the following nature:
!
!     S = A + B E          for E < Elim   | t - tw | < pi/2
!
!         (Elim-E)
!     S = --------         for E > Elim   | t - tw | < pi/2
!            TAU
!
!     S = 0                for E > Elim   | t - tw | > pi/2
!
!     in which the terms A and B are:
!
!         [cf10]          2        2        2              4
!     A = ------ pi (1./g)  [rhoaw]  [cdrag]  (U cos(t-tw))
!          2 pi
!
!     and:
!                              U cos(t-tw)              s
!     B = 5 [cf20] [rhoaw] f { ----------- -  [cf30] } ----
!                                  Cph                 2 pi
!
!     The coefficient TAU in the relaxation model is given by:
!                           2
!                   / 2 pi \      g
!     TAU = [cf40] | -----  | ------------
!                   \  s    /  U cos(d-dw)
!
!     The limiting spectrum is given by:
!
!                    -3
!             ALPHA K                 s   -4    2     2
!     Elim = ----------- exp ( -5/4 (----)   ) --- cos ( t - tw )
!             2  Cg                  spm        pi
!
!     in which:
!
!        ALPHA   wind sea and/or depth dependent proportionality
!                coefficient which controls the energy scale of the
!                limiting spectrum.
!              * In the first-generation model ALPHA is a constant
!                equal to 0.0081 (fully developed)
!              * In the second-generation model ALPHA depends on the
!                energy in the wind sea part of the spectrum. ALPHA
!                is calculated here by:
!                                                [cf60]
!                            ALPHA = [cf50] * Edml
!
!        spm     adapted Pierson-Moskowitz (1964) peak frequency
!
!     The total non-dimensional energy in the wind sea part of the
!     spectrum is calculated by (see subroutine WINDP2):
!
!                   2
!               grav  * ETOTW
!       Edml =  -------------
!                       4
!                 wind10
!
!
!  4. Argument variables
!
! i   SPCDIR: (*,1); spectral directions (radians)                        
!             (*,2); cosine of spectral directions                        
!             (*,3); sine of spectral directions                          
!             (*,4); cosine^2 of spectral directions                      
!             (*,5); cosine*sine of spectral directions                   
!             (*,6); sine^2 of spectral directions                        
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 
!
   REAL    SPCDIR(MDC,6)                                               
   REAL    SPCSIG(MSC)                                                 
!
!  6. Local variables
!
!     IENT  : Number of entries into this subroutine
!
   INTEGER IENT
!
   REAL  :: FPM    ! Pierson-Moskowitz frequency
   REAL  :: SWIND_EXP, SWIND_IMP    ! explicit and implicit part of wind source
!
!        INTEGERS:
!        ---------
!        IDWMIN      Minimum counter for spectral wind direction
!        IDWMAX      Maximum counter for spectral wind direction
!        IX          Counter of gridpoint in x-direction
!        IY          Counter of gridpoint in y-direction
!        IS          Counter of frequency bin
!        ISSTOP      Countrer for the maximum frequency of all directions
!        IDDUM       Dummy counter
!        ID          Counter of directional distribution
!        IDWMIN/IDWMAX  Minimum / maximum counter in wind sector (180 degrees)
!
!        REALS:
!        ---------
!        ALPM        Coefficient for overshoot at deep water
!        ALPMD       Coefficient for overshoot corrected for shallow
!                    water using expression of Bretschnieder (1973)
!        ALIMW       limiting spectrum in terms of action density
!        ARG1, ARG2  Exponent
!        CDRAG       Wind drag coefficient
!        DND         Nondimensional depth
!        DTHETA      Difference in rad between wave and wind direction
!        EDML        Dimensionless energy
!        ETOTW       Total energy of the wind sea part of the spectrum
!        RHOAW       Density of air devided by the density of water
!        SIGPK       Peak frequency in terms of rad /s
!        SIGPKD      Adapted peak frequency for shallow water
!        TAU         Variable for the wind growth equation
!        THETA       Spectral direction
!        THETAW      Mean direction of the relative wind vector
!        TWOPI       Two times pi
!        WIND10      Velocity of the relative wind vector
!
!        one and more dimensional arrays:
!        ---------------------------------
!        AC2       4D    Action density as function of D,S,X,Y and T
!        ALIMW     2D    Limiting action density spectrum
!        DEP2      1D    Depth
!        CGO       2D    Group velocity
!        KWAVE     2D    Wave number
!        LOGSIG    1D    Logaritmic distribution of frequency
!        IMATRA    2D    Coefficients of right hand side of vector
!        IMATDA    2D    Coefficients of the diagonal
!        PLWNDA    3D    Values of source term for test point
!        PLWNDB    3D    Values of source term for test point
!        SPCDIR    1D    Spectral direction of wave component
!        IDCMIN    1D    Minimum counter
!        IDCMAX    1D    Maximum counter in directional space
!        GROWW     2D    Aux. array to determine whether there are
!                        wave generation conditions
!
!        PWIND(1)  = CF10     188.0
!        PWIND(2)  = CF20     0.59
!        PWIND(3)  = CF30     0.12
!        PWIND(4)  = CF40     250.0
!        PWIND(5)  = CF50     0.0023
!        PWIND(6)  = CF60    -0.2233
!        PWIND(7)  = CF70     0.       (not used)
!        PWIND(8)  = CF80    -0.56     (not used)
!        PWIND(9)  = RHOAW    0.00125  (density air / density water)
!        PWIND(10) = EDMLPM   0.0036   (limit energy Pierson Moskowitz)
!        PWIND(11) = CDRAG    0.0012   (drag coefficient)
!        PWIND(12) = UMIN     1.0      (minimum wind velocity)
!        PWIND(13) = PMLM     0.13     (  )
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SOURCE
!
!     6. SUBROUTINES USED
!
!        WINDP2     (compute the total energy in the wind sea part of
!                     the spectrum). Subroutine WINDP2 is called in
!                     SWANCOM1 in subroutine SOURCE
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
!   --------------------------------------------------------------
!   Calculate the adapted peak frequency
!   --------------------------------------------------------------
!   If first-generation model
!     alpha is constant
!   else
!     Calculate energy in wind sea part of spectrum ETOTW
!     Calculate alpha on the basis of ETOTW
!   end
!   --------------------------------------------------------------
!   Take depth effects into account for alpha
!   --------------------------------------------------------------
!   For each frequency and direction
!     compute limiting spectrum and determine whether there is
!     grow or decay
!   end
!   --------------------------------------------------------------
!   Do for each frequency and direction
!     If wind-wave generation conditions are present
!       calculate A + B E
!     else if energy is larger than limiting spectrum
!       calculate dissipation rate with relaxation model
!     endif
!   enddo
!   --------------------------------------------------------------
!   Store results in matrix (IMATRA or IMATDA)
!   --------------------------------------------------------------
!   End of the subroutine WNDPAR
!   --------------------------------------------------------------
!
!     10. SOURCE
!
!***********************************************************************
!
   INTEGER  IS,ID,ITER,IDWMIN,IDWMAX,IDDUM ,ISSTOP
!
   REAL     WIND10,THETA ,THETAW,EDML  ,ARG1  ,ARG2  ,             &
            ALPM  ,ALPMD ,TEMP1 ,TEMP2 ,FACTA ,FACTB ,             &
	    ADUM  ,BDUM  ,CINV  ,SIGTPI,SIGMA ,TWOPI ,TAUINV,      &
	    SIGPK ,SIGPKD,DND   ,ETOTW ,ALIM1D,                    &
	    CTW   ,STW   ,COSDIF,                                  &
	    DIRDIS,AC2CEN,DTHETA
!
   REAL  :: AC2(MDC,MSC,0:MT)
   REAL  :: ALIMW(MDC,MSC)
   REAL  :: IMATDA(MDC,MSC), IMATRA(MDC,MSC)
!  Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   
   REAL  :: KWAVE(MSC,MICMAX)                                          
   REAL  :: PLWNDA(MDC,MSC,NPTST)                                      
   REAL  :: PLWNDB(MDC,MSC,NPTST)
   REAL  :: DEP2(MT)
!  Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   
   REAL  :: CGO(MSC,MICMAX)                                            
!
   INTEGER  IDCMIN(MSC),IDCMAX(MSC)
!
   LOGICAL  GROWW(MDC,MSC)
!
!   SAVE IENT
!   DATA IENT/0/
!   IF (LTRACE) CALL STRACE (IENT,'WNDPAR')
!
!  *** initialization of arrays ***
!
   DO IS = 1, MSC
     DO ID = 1, MDC
       GROWW(ID,IS) = .FALSE.
       ALIMW(ID,IS) = 0.
     ENDDO
   ENDDO
!
!  *** calculate the adapted shallow water peak frequency         ***
!  *** according to Bretschneider (1973) using the nondimensional ***
!  *** depth DND                                                  ***
!
   TWOPI  = 2. * PI_W
!   DND    = MIN( 50. , GRAV_W * DEP2(IGC) / WIND10**2 )
   DND    = MIN( 50. , GRAV_W * DEP2(kcgrd(1)) / WIND10**2 )
   SIGPK  = TWOPI * 0.13 * GRAV_W / WIND10
   SIGPKD = SIGPK / TANH(0.833*DND**0.375)
   FPM    = SIGPKD                                                     
   CTW    = COS(THETAW)                                                
   STW    = SIN(THETAW)                                                
!
   IF(IWIND == 1)THEN
!
!    *** first generation model ***
!
     ALPM = 0.0081
!
   ELSE IF(IWIND == 2)THEN
!
!    *** second generation model ***
!
!    *** Determine the proportionality constant alpha on the basis ***
!    *** of the total energy in the wind sea part of the spectrum  ***
!    *** output of subroutine (WINDP2) is ETOTW                    ***
!
     CALL WINDP2 (IDWMIN  ,IDWMAX  ,SIGPKD  ,FPM     ,             &
                  ETOTW   ,                                        &
		  AC2     ,SPCSIG  ,         WIND10               )    

     EDML = MIN ( PWIND(10) , (GRAV_W**2 * ETOTW) / WIND10**4 )
     EDML = MAX ( 1.E-25 , EDML )
!
     ARG1 = ABS(PWIND(6))
     ALPM = MAX( 0.0081, (PWIND(5) * (1./EDML)**ARG1) )
!
   ENDIF
!
!  *** Take into account depth effects for proportionality ***
!  *** constant alpha through the nondimensional depth DND ***
!
   ALPMD  = 0.0081 + ( 0.013 - 0.0081 ) * EXP ( -1. * DND )
   ALPM   = MIN ( 0.155  ,  MAX ( ALPMD , ALPM ) )
!
!  *** Calculate the limiting spectrum in terms of action density   ***
!  *** for the wind sea part (centered around the local wind        ***
!  *** direction). For conversion of f^-5 --> k^-3 and coefficients ***
!  *** see Kitaigorodskii et al. 1975                               ***
!
   DO IS = 1, ISSTOP
     TEMP1  = ALPM / ( 2. * KWAVE(IS,1)**3 * CGO(IS,1) )
     ARG2   = MIN ( 2. , SIGPKD / SPCSIG(IS) )                         
     TEMP2  = EXP ( (-5./4.) * ARG2**4 )
     ALIM1D = TEMP1 * TEMP2 / SPCSIG(IS)                               
     DO IDDUM = IDWMIN, IDWMAX
       ID     = MOD ( IDDUM - 1 + MDC, MDC ) + 1
       THETA  = SPCDIR(ID,1)                                           
       COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                    
!
!  For better convergence the first guess of the directional spreading 
!  is modified in third generation mode. The new formulation better    
!  fits the directional spreading of the deep water growth curves.     
!
       IF((ITER == 1).AND.(IGEN == 3))THEN                           
         DIRDIS = 0.434917 * (MAX(0., COSDIF))**0.6                    
       ELSE                                                            
         DIRDIS = (2./PI_W) * COSDIF**2                                  
       END IF                                                          
!
       ALIMW(ID,IS) = ALIM1D * DIRDIS
!       AC2CEN       = AC2(ID,IS,IGC)
       AC2CEN       = AC2(ID,IS,kcgrd(1))
       IF(AC2CEN <= ALIMW(ID,IS))THEN
         GROWW(ID,IS) = .TRUE.
       ELSE
         GROWW(ID,IS) = .FALSE.
       ENDIF
     ENDDO
!    *** test output ***
!     IF(TESTFL .AND. ITEST >= 10)THEN
!       WRITE(PRINTF,2002) IS, SPCSIG(IS), KWAVE(IS,1), CGO(IS,1)       
!2002   FORMAT(' WNDPAR: IS SPCSIG KWAVE CGO :',I3,3E12.4)              
!       WRITE(PRINTF,2003) TEMP1, TEMP2, ARG2
!2003   FORMAT(' WNDPAR: TEMP1 TEMP2 ARG2    :',3X,3E12.4)
!     END IF
   ENDDO
!
!  *** Calculate the wind input (linear term A and exponential  ***
!  *** term B) in wave generating conditions or disspation term ***
!  *** if energy in bin is larger than limiting spectrum        ***
!
!
   FACTA = PWIND(1) * PI_W * PWIND(9)**2 * PWIND(11)**2 / GRAV_W**2
!
   DO IS = 1, ISSTOP
     SIGMA   = SPCSIG(IS)                                              
     SIGTPI  = SIGMA * TWOPI
     CINV    = KWAVE(IS,1) / SIGMA
     FACTB = PWIND(2) * PWIND(9) * SIGMA / TWOPI                       
     DO IDDUM = IDCMIN(IS), IDCMAX(IS)
       ID     = MOD ( IDDUM - 1 + MDC, MDC ) + 1
       DTHETA = SPCDIR(ID,1) - THETAW                                  
       COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                    
!       AC2CEN = AC2(ID,IS,IGC)
       AC2CEN = AC2(ID,IS,kcgrd(1))
!
       SWIND_EXP = 0.                                                  
       SWIND_IMP = 0.                                                  

       IF(GROWW(ID,IS))THEN
!        *** term A ***
         IF(SIGMA >= ( 0.7 * SIGPKD ))THEN
           ADUM = FACTA * (WIND10 * COSDIF)**4                         
           ADUM = MAX ( 0. , ADUM / SIGTPI )
         ELSE
           ADUM = 0.
         END IF
!        *** term B; Note that BDUM is multiplied with a factor 5 ***
!        *** as in the DOLPHIN-B model                            ***
!
         BDUM = MAX( 0., ((WIND10 * CINV) * COSDIF-PWIND(3)))          
         BDUM = FACTB * BDUM * 5.
         SWIND_EXP = ADUM + BDUM * AC2CEN                              
!
       ELSE IF(.NOT. GROWW(ID,IS) .AND. AC2CEN > 0.)THEN
!
!        *** for no energy dissipation outside the wind field     ***
!        *** TAUINV is set equal zero (as in the DOLPHIN-B model) ***
!
         IF(COSDIF < 0.) THEN                                    
           TAUINV = 0.
         ELSE
           TAUINV = ( SIGMA**2 * WIND10 * ABS(COSDIF) ) /        &
	            ( PWIND(4) * GRAV_W * TWOPI**2 )
         ENDIF
         SWIND_EXP = TAUINV * ALIMW(ID,IS)
         SWIND_IMP = -TAUINV
         ADUM = ALIMW(ID,IS)
         BDUM = TAUINV
       END IF
!
!      *** store results in IMATDA and IMATRA ***
!
       IMATRA(ID,IS) = IMATRA(ID,IS) + SWIND_EXP
!       IMATDA(ID,IS) = IMATDA(ID,IS) - SWIND_IMP
       IF (TESTFL) PLWNDA(ID,IS,IPTST) = SWIND_EXP                     
       IF (TESTFL) PLWNDB(ID,IS,IPTST) = SWIND_IMP                     
!
!
!      *** test output ***
!
!      Value of ITEST changed from 10 to 110 to reduce test output     
!       IF ( TESTFL .AND. ITEST .GE. 110 ) THEN                        
!         WRITE(PRINTF,2004) IS, ID, GROWW(ID,IS), ADUM, BDUM           
!2004     FORMAT(' WNDPAR: IS ID GROWW ADUM BDUM     :',            &
!                2I3,2X,L1,2X,2E12.4)
!       END IF
     ENDDO
   ENDDO
!
!  *** test output ***
!
!  Value of ITEST changed from 10 to 60 to reduce test output          
!   IF(TESTFL .AND. ITEST >= 60)THEN                              
!     WRITE(PRINTF,*)
!     WRITE(PRINTF,6051) IDWMIN, IDWMAX
!6051 FORMAT(' WNDPAR : IDWMIN IDWMAX     :',2I5)
!     WRITE(PRINTF,6052) THETAW,WIND10,SIGPK,SIGPKD
!6052 FORMAT(' WNDPAR : Tw U10 Spk Spk,d   :',4E12.4)
!     WRITE(PRINTF,7050) ETOTW, EDML, ALPM, ALPMD
!7050 FORMAT(' WNDPAR: ETOW EDML ALPM ALPMD:',4E12.4)
!   ENDIF
!
   RETURN

   END SUBROUTINE WNDPAR

