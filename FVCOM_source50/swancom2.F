
!     SWAN/COMPU      file 2 of 5
!
!     PROGRAM SWANCOM2.FOR
!
!     This file SWANCOM2 of the main program SWAN
!     include the next subroutines (mainly subroutines for
!     the source terms for dissipation and some general stuff):
!
!     DISSIPATION SOURCE TERMS :
!
!     SBOT    (Bottom friction)
!     SVEG    (Dissipation due to vegetation)                             40.55

   SUBROUTINE SBOT (ABRBOT  ,DEP2    ,ECOS    ,ESIN    ,KWAVE   , &
                    SPCSIG  ,UBOT    ,UX2     ,UY2     ,IDCMIN  , &
		    IDCMAX  ,ISSTOP  ,VARFR   ,FRCOEF  ,IG      , &
		    IMATRA                             )
!  (This subroutine has been tested only for case IBOT = 1)
!
!****************************************************************
!
!     Computation of the source terms due to bottom friction
!
!  Method
!
!     In SWAN several bottom dissipation models are computed, i.e.:
!
!     IBOT = 1   Jonswap bottom friction model
!     IBOT = 2   Collins bottom friction model
!     IBOT = 3   Madsen bottom friction model (see Tolman)
!
!     Both methods are implemented in SWAN and the user has to make
!     a choice in the input file.
!
!     1. Jonswap model:
!     -----------------
!
!     The bottom interaction term SEbf(s,d) is supposed to take the
!     Jonswap form (Hasselman et al. 1973):
!                         2
!                    sigma  E(s,d)
!     SEbf = -GAMMA ----------------
!                     2     2
!                    g  sinh  (kD)
!                                                                 2 -3
!     where GAMMA is the decay parameter ,(default GAMMA = 0.067 m s  ).
!     In the Jonswap form the current velocities are not taken into
!     account.
!
!     2. COLLINS model:
!     -----------------
!
!     The energy dissipation due to bottom friction is modelled
!     according the quadratic friction law:
!                    2
!      SE = Tau * |U|
!
!      which for a spectrum can be written as:
!                          2
!                     sigma
!      SE(s,d)= - ---------------- * (Cfw.Ub + Cfc.Uc) * E(s,d)
!                       2
!                 g sinh (K(s) * D)
!
!     Ub is the velocity due to the wave at the bottom
!
!     The current velocity is Uc
!
!     2. MADSEN formulation:
!     ----------------------
!
!     The bottom dissipation applying Madsen formulation is as
!     follows:
!
!                          fw [n - 1/2] UBR E(s,d)
!     [1]    Sdsb(s,d) = -  ------------------------
!                                      D
!
!     in which :
!                            2
!                           s * D
!     [1a]   (n - 1/2) = -------------
!                                2
!                        2 g sinh (kD)
!
!     UBOT(IX,IY) is computed in the subroutine SINTGRL. The friction
!     factor fw is estimated using the formulation of Jonsson (1963,
!     1966a):
!
!                1                1                        Ab,r
!     [2]     -------- + log  { ---------- } = mf + log  { ----- }
!            4 sqrt(fw)     10  4 sqrt(fw)             10   Kn
!
!     with:
!
!               2        //      1
!     [3]   Ab,r  = 2 * // -------------- E(s,d) ds dd
!                      //      2
!                          sinh (kD)
!
!     with: Ab,r is the representative near bottom excursion
!                amplitude
!           Kn   equivalent roughness
!           mf   constant ( mf = -0.08) (determined by Jonsson
!                                        and Carlssen 1976 )
!
!     [2] is only valid for Ab,r/Kn larger than approximately 1.
!     For smaller values a constant value of fw is used (fw = 0.3
!     for Ab,r/Kn < 1.57 )
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4  
!   USE ALL_VARS, ONLY : MT,AC2                                                      
   USE VARS_WAVE, ONLY : MT,AC2                                                      

   IMPLICIT NONE
   
   REAL    SPCSIG(MSC)                                                 
   INTEGER  ID,IS,ISSTOP,J,IDDUM,IG
   REAL :: XDUM,KD,SBOTEO,CFBOT,FACB,CFW,FW,CURR,UC,ABRBOT,ADUM,CDUM,DDUM
   LOGICAL  VARFR
   REAL :: DEP2(MT),ECOS(MDC),ESIN(MDC),IMATDA(MDC,MSC),IMATRA(MDC,MSC), &
           KWAVE(MSC,ICMAX),PLBTFR(MDC,MSC,NPTST),UBOT(MT),  &
	   UX2(MT),UY2(MT),DISSC1(MDC,MSC),FRCOEF(MT)  
   INTEGER :: IDCMIN(MSC),IDCMAX(MSC)
   REAL    :: AKN

   IF(IBOT >= 1 .AND. DEP2(IG) > 0.)THEN
     IF(IBOT == 1)THEN
!
!      *** Jonswap model ***
!
!      PBOT(3) = GAMMA (a) in the Jonswap equation,
!
       CFBOT = PBOT(3) / GRAV_W**2
       
!       CFBOT = CFBOT*10.0               !JQI added
     ELSE IF (IBOT == 2)THEN
!
!      *** Collins model ***
!
!      PBOT(2) = [cfw]
!
       IF(VARFR)THEN                                                
         CFW = FRCOEF(IG)
       ELSE
         CFW = PBOT(2)
       ENDIF
       CFBOT = CFW * UBOT(IG) / GRAV_W
     ELSE IF(IBOT == 3)THEN
!
!      *** Madsen model ***
!
       IF(VARFR)THEN                                                 
         AKN = FRCOEF(IG)
       ELSE
         AKN = PBOT(5)
       ENDIF
!
!        *** PBOT(4) = Mf                      ***
!        *** AKN = PBOT(5) = [kn]  (roughness) ***
!
       IF((ABRBOT/AKN) > 1.57)THEN
         XDUM = PBOT(4) + LOG10 ( ABRBOT / AKN )
!
!        *** solving the implicit equation using a Newton ***
!        *** Rapshon iteration proces : a + log a = b     ***
!        *** the start value for ADUM = 0.3 because 0.3626 ***
!        *** is the minimum value of ADUM with b=-0.08.    ***
!
         ADUM = 0.3
         DO J = 1, 50
           CDUM  = ADUM
           DDUM  = ( ADUM + LOG10(ADUM) - XDUM ) / ( 1.+ ( 1. / ADUM) )
           ADUM  = ADUM - DDUM
           IF(ABS(CDUM - ADUM) < 1.E-4) GOTO 29
         END DO	   
         WRITE(*,*) ' error in iteration fw: Madsen formulation'
29       CONTINUE
!                                                 1               1
!        *** computation of FW -->  A = ----- --> FW = -----
!                                              4 uFW          16 A**2
         FW = 1. / (16. * ADUM**2)
       ELSE
         FW = 0.3
       ENDIF
       CFBOT =  UBOT(IG) * FW / (SQRT(2.) * GRAV_W)
     ENDIF

     DO IS = 1, ISSTOP
       KD      = KWAVE(IS,1) * DEP2(IG)
       IF(KD < 10.)THEN
         FACB = CFBOT * (SPCSIG(IS) / SINH(KD)) **2                   
!
         DO IDDUM = IDCMIN(IS) , IDCMAX(IS)
           ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!
           SBOTEO = FACB                                               
           IF(IBOT == 2 .AND. ICUR == 1 .AND. PBOT(1) > 0.)THEN
!            additional dissipation due to current, seldom used
             CURR = UX2(IG)*ECOS(ID) + UY2(IG)*ESIN(ID)
             UC   = ABS(CURR)
!             PBOT(1) = [cfc]
             SBOTEO = FACB + PBOT(1) * UC * (SPCSIG(IS) / SINH(KD)) **2 
           END IF
!
!          *** store the results in the array IMATDA             ***
!
           IMATRA(ID,IS) = IMATRA(ID,IS) - SBOTEO*AC2(ID,IS,IG)
           IMATDA(ID,IS) = IMATDA(ID,IS) + SBOTEO
!           IF (TESTFL) PLBTFR(ID,IS,IPTST) = -1.* SBOTEO               
!           DISSC1(ID,IS) = DISSC1(ID,IS) + SBOTEO
	 END DO   !IDDUM  
       END IF
     END DO   !IS  
!
   ENDIF
!
   RETURN
   END SUBROUTINE SBOT
 
#  if defined (VEGETATION)
!
!****************************************************************
!
      SUBROUTINE SVEG ( DEP2   ,IMATDA   ,ETOT   ,SMEBRK    ,  KMESPC&
           & ,PLVEGT   , IDCMIN ,IDCMAX   ,ISSTOP ,DISSC1    , NPLA2 &
           & ,      IG , IMATRA )

!
!****************************************************************
!
      USE SWCOMM2
      USE SWCOMM3
      USE SWCOMM4
      USE OCPCOMM4
      USE M_GENARR
!
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
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2017  Delft University of Technology
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
!     40.55: Bastiaan Burger, Martijn Meijer
!     40.61: Marcel Zijlema
!     40.58: Tomo Suzuki, Marcel Zijlema
!
!  1. Updates
!
!     40.55, May. 05: implementation of vegetation dissipation formula
!     40.61, Sep. 06: introduce DISVEG variable for output purposes
!     40.58, Nov. 08: some modifications and corrections
!
!  2. Purpose
!
!     Computation of the source term due to vegetation dissipation
!
!  3. Method
!
!     The energy dissipation due to vegetation is described by a
!     Morrison type equation, modelling the plants as vertical,
!     noncompliant cylinders, neclecting swaying motions induced by
!     waves. Vegetation characteristics that are used as input are
!     drag coefficient, vegetation height, plant density and diameter.
!
!     The formula used in SWAN is due to Dalrymple (1984)
!     (see Mendez and Losada, 2004):
!
!     d(Cg E)
!     ------- = -epsv
!        dx
!
!     with dissipation due to vegetation:
!
!     epsv = 1/2 * 1/sqrt(pi) * rho * Cd * bv * Nv *
!                                                                3
!                              (gk/2sigma)^3 * ((A + B)/C) * Hrms
!
!     where
!
!     rho   = water density
!     k     = wave number
!     sigma = angular frequency
!     Cd    = drag coefficient
!     bv    = stem thickness
!     Nv    = vegetation density
!     Hrms  = rms wave height
!
!     and the coefficients
!
!     A = sinh^3 k*ah
!     B = 3*sinh k*ah
!     C = 3k*cosh^3 kh
!
!     with ah the vegetation height
!
!     The source term to be used in SWAN is based on the
!     corresponding dissipation rate and reads
!
!     Dtot = -epsv / rho / g
!
!     In the formulation, the mean average wavenumber according
!     to the WAM-formulation and the mean frequency will be employed
!
!     Now, the source term is:
!
!                      E
!     Sveg = Dtot *  ------ = factor * sqrt(Etot) * E
!                     Etot
!
!     and is linearized by means of the Picard iteration
!
!
!  4. Argument variables
!
!     DEP2        water depth
!     DISSC1      dissipation coefficient
!     ETOT        total energy per spatial gridpoint
!     IDCMIN      frequency dependent counter in directional space
!     IDCMAX      frequency dependent counter in directional space
!     IMATDA      coefficients of diagonal of matrix
!     ISSTOP      maximum counter of wave component in frequency
!                 space that is propagated
!     KMESPC      mean average wavenumber according to the WAM-formulation
!     NPLA2       number of plants per square meter (depth-averaged)
!     PLVEGT      array containing the vegetation source term for test-output
!     SMEBRK      mean frequency according to first order moment
!
      INTEGER ISSTOP, IDCMIN(MSC), IDCMAX(MSC)
      REAL    DEP2(MT) , IMATDA(MDC,MSC),DISSC1(MDC,MSC),PLVEGT(MDC,MSC,NPTST),NPLA2(MT)
      REAL    ETOT, SMEBRK, KMESPC,IMATRA(MDC,MSC)
!
!  6. Local variables
!
!     A     :     auxiliary variable
!     B     :     auxiliary variable
!     C     :     auxiliary variable
!     D     :     auxiliary variable
!     ID    :     counter of the spectral direction
!     IDDUM :     counter
!     IENT  :     number of entries
!     IK    :     counter
!     IL    :     counter
!     IS    :     counter of relative frequency band
!     KD    :     wavenumber times water depth
!     KVEGH :     wavenumber times plant height
!     LAYPRT:     part of layer below water level
!     SINHK :     sinh(kh)
!     SLAYH :     total sum of layer thicknesses
!     SLAYH1:     sum of layer thicknesses below water level
!     SLAYH2:     sum of layer thicknesses below water level
!     SVEG1 :     layer-independent dissipation factor
!     SVEG2 :     total sum of dissipation factor over layers
!     SVEGET:     source term containing dissipation due to vegetation
!                 to be stored in the array IMATDA
!
      INTEGER ID, IDDUM, IENT, IK, IL, IS, IG
      REAL    A, B, C, D, KD, KVEGH, LAYPRT, SINHK, SLAYH, SLAYH1, SLAYH2, SVEG1, SVEG2, SVEGET
!
!  9. Subroutines calling
!
!     SOURCE
!
! 12. Structure
!
!     Vegetation parameters are given per layer thickness, so for
!     each layer the contribution to wave damping is calculated
!
!     This routine checks in which layer the water level is present
!
!            -----------
!                        ILMAX
!            -----------
!              _ d       2
!            --|--------
!              |         1
!            --|--------
!
!     d     = waterdepth
!     ILMAX = number of layers in grid point
!
!     Subsequently, the vegetation parameters up to the layer where the
!     water level is in, are used to calculate dissipation for each layer
!
!     Thereafter, the contributions to disspation are summed up
!
!     With this summation the total dissipation due to vertical varying
!     vegetation is calculated
!
! 13. Source text
!

!     --- compute layer-independent vegetation dissipation factor

      KD    = KMESPC * DEP2(IG)
      IF ( KD.GT.10. ) RETURN
      C     = 3.*KMESPC*(COSH(KD))**3
      SVEG1 = SQRT(2./PI_W) * GRAV_W**2 * (KMESPC/SMEBRK)**3 * SQRT(ETOT)/ C

!      IF ( VARNPL ) SVEG1 = SVEG1 * NPLA2(KCGRD(1))
      SVEG1 = SVEG1 * SUM(VEGNSL(IG,1:ILMAX))

!     --- compute dissipation factor for each layer and summed up

      SLAYH = 0.
      DO IL = 1, ILMAX
         SLAYH = SLAYH + LAYH(IG,IL)
      ENDDO

      KVEGH = 0.
      C     = 0.
      D     = 0.
      SVEG2 = 0.

      IF ( DEP2(IG).GT.SLAYH ) THEN

         DO IL = 1, ILMAX
            KVEGH = KVEGH + KMESPC * LAYH(IG,IL)
            SINHK = SINH(KVEGH)
            A     = C
            B     = D
            C     = SINHK**3
            D     = 3.*SINHK
            A     = C - A
            B     = D - B
            SVEG2 = SVEG2 + VEGDRL(IG,IL)*VEGDIL(IG,IL)*VEGNSL(IG,IL)*(A + B)
         END DO

      ELSE IF ( DEP2(IG).LT.LAYH(IG,1) ) THEN

         SINHK = SINH(KD)
         A     = SINHK**3
         B     = 3.*SINHK
         SVEG2 = VEGDRL(IG,1)*VEGDIL(IG,1)*VEGNSL(IG,1)*(A + B)

      ELSE

         SLAYH1 = 0.
         SLAYH2 = 0.
         LAYPRT = 0.
         VGLOOP : DO IL = 1, ILMAX
            SLAYH1 = SLAYH1 + LAYH(IG,IL)
            IF (DEP2(IG).LE.SLAYH1) THEN
               DO IK = 1, IL-1
                  SLAYH2 = SLAYH2 + LAYH(IG,IK)
               END DO
               LAYPRT = DEP2(IG) - SLAYH2
               DO IK = 1, IL-1
                  KVEGH = KVEGH + KMESPC * LAYH(IG,IK)
                  SINHK = SINH(KVEGH)
                  A     = C
                  B     = D
                  C     = SINHK**3
                  D     = 3.*SINHK
                  A     = C - A
                  B     = D - B
                  SVEG2 = SVEG2+VEGDRL(IG,IK)*VEGDIL(IG,IK)*VEGNSL(IG,IK)*(A + B)
               END DO
               KVEGH = KVEGH + KMESPC * LAYPRT
               SINHK = SINH(KVEGH)
               A     = C
               B     = D
               C     = SINHK**3
               D     = 3.*SINHK
               A     = C - A
               B     = D - B
               SVEG2 = SVEG2 + VEGDRL(IG,IL)*VEGDIL(IG,IL)*VEGNSL(IG,IL)*(A + B)
               EXIT VGLOOP
            END IF
         END DO VGLOOP

      END IF

!     --- compute total dissipation

      SVEGET = SVEG1 * SVEG2
!
!     *** test output ***
!
      IF (TESTFL .AND. ITEST.GE.60) THEN
         WRITE (PRTEST, 110) IVEG, IG, DEP2(IG), SVEGET
 110     FORMAT (' SVEG :IVEG INDX DEP VEGFAC:', 2I5, 2E12.4)
      END IF

      DO IS = 1, ISSTOP
         DO IDDUM = IDCMIN(IS), IDCMAX(IS)
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1

!           *** store the results in the array IMATDA ***
!           *** if testfl store results in array for isoline plot ***

            IMATRA(ID,IS) = IMATRA(ID,IS) - SVEGET*AC2(ID,IS,IG)
            IMATDA(ID,IS) = IMATDA(ID,IS) + SVEGET
!            IF (TESTFL) PLVEGT(ID,IS,IPTST) = -1.* SVEGET
            DISSC1(ID,IS) = DISSC1(ID,IS) + SVEGET

         END DO
      END DO

      RETURN
    END SUBROUTINE SVEG
#   endif
!
!****************************************************************
!
   SUBROUTINE FRABRE ( HM, ETOT, QBLOC )                               
!
!****************************************************************
!
!     to compute the fraction of breaking waves in point ix,iy
!     of the computational grid
!
!  Method (updated...)
!
!      The fraction of breaking waves in a point ix,iy is given by
!      the implicit relation:
!
!        1 - Qb        ETOT
!        ------ = -8 * -----
!        ln Qb         HM**2
!
!        from which Qb can be found by solving the equation:
!
!                         ETOT
!        F = 1 - Qb + 8 * ----  * ln(Qb) = 0.
!                           2
!                         HM
!
!        The following appproximation is applied:
!
!                            2
!  (1)|   B = sqrt( 8 ETOT/HM ), i.e. B = Hrms/HM
!
!
!     |   Qo = 0.                                      B <= 0.5
!  (2)|                 2
!     |   Qo = ( 2B -1 )                         0.5 < B <= 1
!
!
!     applying the Newton-Raphson procedure (for 0.2<B<1.0):
!
!     |   Qb = 0.                                      B <= 0.2
!     |
!     |                                   2
!     |               2  Qo - exp((Qo-1)/B )
!  (3)|   Qb = Qo  - B   ------------------      0.2 < B <  1.0
!     |                   2               2
!     |                  B  - exp((Qo-1)/B )
!     |
!     |
!     |   Qb = 1.                                      B >= 1.0
!     |
!
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
!
   IMPLICIT NONE

   REAL    :: ETOT,  HM,  QBLOC
   INTEGER :: IENT
   REAL    :: B,  B2,  QO,  Z

   IF((HM > 0.) .AND. (ETOT >= 0.))THEN
     B = SQRT(8. * ETOT / (HM*HM) )
   ELSE
     B = 0.0
   END IF
!
   IF(B <= 0.5)THEN
     QO = 0.
   ELSE IF(B <= 1.0)THEN
     QO = (2.*B - 1.)**2
   END IF
!
   IF(B <= 0.2)THEN
     QBLOC = 0.0
   ELSE IF(B < 1.0)THEN
!
!    *** second iteration to find Qb ***
!
     B2 = B*B
     Z  = EXP((QO-1.)/B2)
     QBLOC = QO - B2 * (QO-Z)/(B2-Z)
   ELSE
     QBLOC = 1.0
   END IF
!
   RETURN
   END SUBROUTINE FRABRE
!
!****************************************************************
!
   SUBROUTINE SSURF (ETOT    ,HM      ,QB      ,SMEBRK  ,       &
                     IMATRA  ,IMATDA  ,IDCMIN  ,IDCMAX  ,       &
		     PLWBRK  ,ISSTOP  ,IG               )
!
!****************************************************************
!     Computation of the source term due to wave breaking.
!     White capping is not taken into account
!
!  Method
!
!     The source term for surf breaking is implemented following
!     the approach of Battjes/Janssen (1978) for the energy dissipation:
!
!             Alpha      -     2                  -   SMEBRK
!     Dtot =  ----  Qb * f * Hm              with f = ------
!              4                                      2 * Pi
!
!     Now the source term is:
!
!                      SIGMA * AC2(ID,IS,IX,IY)
!     Sbr =  - Dtot *  ------------------------  =
!                              Etot
!
!
!              Alpha * SMEBRK * Qb * Hm * Hm    SIGMA * AC2(ID,IS,IX,IY)
!         =  - ------------------------------ * -------------------------
!                       8 * Pi                            Etot
!
!
!         =  WS * SIGMA * AC2(ID,IS,IX,IY)   =  WS * E
!
!
!     with
!
!     Alpha = PSURF(1)                            ;
!
!                   SMEBRK Qb
!     WS    = Alpha ------ --                     ;
!                     Pi   BB
!                        2
!     BB    = 8 Etot / Hm  = - (1 - Qb) / ln (Qb) ;
!
!
!     The local maximum wave height Hm and mean frequency SMEBRK are computed
!     in subroutine SINTGRL.
!     The fraction of breaking waves Qb is calculated in the subroutine FRABRE
!
!     The new value for the dissipation is computed implicitly using
!     the last computed value for the action density Nold (at the spatial
!     gridpoint under consideration).
!
!     Sbr = WS * N
!
!         = Sbr_new + (d Sbr/d N) (Nnew - Nold)
!
!         = WS * Nnew + SbrD * (Nnew - Nold)
!
!         = (WS + SbrD)* Nnew - SbrD * Nold
!
!         = SURFA1 * Nnew - SURFA0 * Nold
!
!     In order to do this we need the derivative
!     of the source term Sbr to the action density N
!
!             d Sbr     d WS
!     SbrD =  -----  =  ---- * N + WS
!             d N       d N
!
!     Since BB and SMEBRK * N are proportional, we have
!
!     d Sbr     d WS                   SMEBRK  (d Qb/ d BB) *BB - Qb
!     -----  =  ---- * BB + WS = Alpha ------  --------------------- * BB + WS =
!     d N       d BB                     Pi           sqr(BB)
!
!
!           SMEBRK d Qb
!     Alpha ------ ----
!            Pi    d BB
!
!     With:
!
!     d Qb         1
!     ---- = -------------                 ;
!     d BB   (d BB / d Qb)
!
!                      2
!     d Qb           ln (Qb)
!     ---- = ---------------------------
!     d BB   ln (Qb) + (1 - Qb) (1 / Qb)
!
!            Qb (1 - Qb)
!          = ------------                  ;
!            BB (BB - Qb)
!
!     ------------------------------------------------------------
!     Get HM, QB and ETOT from the subroutine SINTGRL
!     For spectral direction IS and ID do,
!       get the mean energy frequency average over the full spectrum
!       If ETOT > 0 then
!         compute source term for energy dissipation SURFA0 and SURFA1
!       Else
!         source term for wave breaking is 0.
!       End if
!       ----------------------------------------------------------
!       Compute source terms for energy averaged frequency
!       Store results in the arrays IMATDA and IMATRA
!     ------------------------------------------------------------
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4 
!   USE ALL_VARS, ONLY : MT,AC2                                                       
   USE VARS_WAVE, ONLY : MT,AC2                                                       
!
   IMPLICIT NONE

   INTEGER :: ISSTOP,IDCMIN(MSC),IDCMAX(MSC)
   REAL    :: DISSC0(MDC,MSC),DISSC1(MDC,MSC),     &
              IMATDA(MDC,MSC),IMATRA(MDC,MSC),PLWBRK(MDC,MSC,NPTST)
   REAL    :: ETOT,HM,QB,SMEBRK                         
   INTEGER :: ID,IDDUM,IENT,IS,IG
   DOUBLE PRECISION BB,DIS0,SbrD,SURFA0,SURFA1,WS
!
!     ALFA = PSURF(1)   <default = 1.0>
!
   BB = 8. * DBLE(ETOT) / ( DBLE(HM)**2 )                              
!   SURFA0 = 0.
!   SURFA1 = 0.
!   IF(REAL(BB) > 0. .AND. REAL(ABS(BB - DBLE(QB))) > 0.)THEN  
!     IF(BB < 1.)THEN
!       WS  = ( DBLE(PSURF(1)) / DBLE(PI_W)) * DBLE(QB) * DBLE(SMEBRK) / BB 
!       SbrD = WS * (1. - DBLE(QB)) / (BB - DBLE(QB))                   
!     ELSE
!       WS  = ( DBLE(PSURF(1)) / DBLE(PI_W)) * DBLE(SMEBRK)               
!       SbrD = 0.
!     END IF
!     SURFA0 = SbrD
!     SURFA1 = WS + SbrD
!   ELSE
!     SURFA0 = 0.
!     SURFA1 = 0.
!   END IF
   IF(REAL(BB) > 0. .AND. REAL(ABS(BB - DBLE(QB))) > 0.)THEN
       WS  = ( DBLE(PSURF(1)) / DBLE(PI_W)) * DBLE(QB) * DBLE(SMEBRK) / BB
   ELSE
       WS = 0.
   ENDIF

!
!  *** store the results for surf wave breaking  ***
!  *** in the matrices IMATDA and IMATRA         ***
!
   DO IS = 1, ISSTOP
     DO IDDUM = IDCMIN(IS), IDCMAX(IS)
       ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
       IMATDA(ID,IS) = IMATDA(ID,IS) + REAL(WS)
       DIS0 = WS * DBLE(AC2(ID,IS,IG))
       IMATRA(ID,IS) = IMATRA(ID,IS) - REAL(DIS0)
     END DO
   END DO


   RETURN
   END SUBROUTINE SSURF

!
!****************************************************************
!
   SUBROUTINE SWCAP  (SPCDIR  ,SPCSIG  ,KWAVE   ,IDCMIN  ,      &
                      IDCMAX  ,ISSTOP  ,ETOT    ,IMATDA  ,      &
		      IMATRA  ,PLWCAP  ,CGO     ,UFRIC   ,      &
		      DEP2    ,DISSIP  ,DISIMP  ,IG      )              
!  (This subroutine has not been fully tested. Only tested for 
!   case IWCAP = 2)
!
!****************************************************************
!
!     Calculates the dissipation due to whitecapping
!
!  Method
!
!     Whitecapping dissipation is formulated as follows:
!
!     S_wc(sig,th) = - C_wc E(sig,th)
!
!     where the coefficient C_wc has four basic forms:
!
!     C_wc1 = C_K  sig~ (k/k~): According to Komen (generalised)
!     C_wc2 = C_BJ sig~ (k/k~): According to Battjes-Janssen (modified)
!     C_wc3 = C_LH            : According to Longuet Higgins (1969)
!     C_wc4 = C_CSM           : According to cumulative steepness
!                               method (Ris et al. (1999) and Donelan (1999))
!
!     In these formulations C_K is defined as (Komen; 1994 p. 145):
!
!                                       n1           n2
!     C_K = C1 [(1-delta) + delta (k/k~)  ] (S~/S_PM)
!
!     where C1, delta, n1 and n2 can be varied
!
!
!     C_BJ is defined as:
!
!                    2
!            alpha Hm Qb
!     C_BJ = -----------
!              8 Pi m0
!
!     where alpha can be varied.
!
!     for Hrms > Hm the formulation changes in a limit to (Hrms->Hm; Qb->1):
!
!            alpha
!     C_BJ = -----
!             Pi
!
!     and C_LH is defined as (Komen; 1994 p. 151-152):
!
!                            4   2                        2
!     C_LH = C3 Sqrt[(m0 sig0 )/g ] exp(A) sig0 (sig/sig0)
!
!     where
!                    2    2         4                2       2
!     A = -1/8 [1-eps ] [g /(m0 sig0 )]  with  [1-eps ] = [m2 ] / [m0 m4]
!
!     and C3 can be varied
!
!
!     Cumulative steepness method. Basic idea is that dissipation depends
!     on the steepness of the wave spectrum at and below a particular
!     frequency. Details can be found in "SWAN fysica plus" (2002),
!     report H3937/A832, implemented by Delf Hydraulics and Alkyon by
!     order of RIKZ/RWS as a part of the project HR-Ontwikkeling.
!                                                                         40.41
!     Recent modifications are described in:                              40.41
!     Hurdle, D.P., and G.Ph. van Vledder, 2004: Improved spectral wave   40.41
!     modelling of white-capping dissipation in swell sea systems. Proc.  40.41
!     Offshore Mechanics, Arctic Engineering Conference, June 2004,       40.41
!     Vancouver, Canada,  paper OMAE2004-51562.                           40.41
!                                                                         40.41
!     Default values: C_st = 4                                            40.41
!                     m    = 2                                            40.41
!                                                                         40.41
!     C_CSM = A*C_st S(sig,th)                                            40.41
!
!     where
!
!     C_st is a tuneable coefficient and
!
!     S(sig,th) = integrals from 0 to sig and from 0 to 2pi of
!
!                 k^2 |cos(th-th')|^m E(sig,th) dsig dth'
!
!     where coefficient m controls the directional dependence in case
!     of straining mechanism (if dominant then m = 2 else m > 10) and
!     th is the actual direction where straining mechanism takes
!     place.
!                                                                         40.41
!     A is a normalisation coefficient computed                           40.41
!                                                                         40.41
!              1     GAMMA(0.5*m+1.)                                      40.41
!       A = -------- ----------------                                     40.41
!           SQRT(PI) GAMMA(0.5*m+0.5)                                     40.41
!
!     In these equations the variables have the following meaning:
!
!     Hm   : Maximum wave height
!     Hrms : Root mean square of the wave heights
!     eps^2: Measure for the spectral bandwidth
!     m0   : Total wave energy density (=ETOT)
!     m2   : Second moment of the variance spectrum (=ETOT2)
!     m4   : Fourth moment of the variance spectrum (=ETOT4)
!     k    : Wave number (=KWAVE(IS,1))
!     k~   : Mean wave number
!     Qb   : Fraction of breaking waves
!     sig  : Frequency (=SPCSIG(IS))
!     sig0 : Average zero crossing frequency
!     sig~ : Mean frequency
!     S~   : Overall steepness (STP_OV)
!     S_PM : Overall steepness for a Pierson-Moskowitz spectrum
!     th   : direction theta (=SPCDIR(ID))
!
!     ------------------------------------------------------------
!     Calculate needed parameters
!     ------------------------------------------------------------
!     If IWCAP = 1, 2, or 5; Calculate C_K
!     If IWCAP = 4 or 5; Calcualte C_BJ
!     If IWCAP = 3; Calculate C_LH
!     If IWCAP = 6; Calculate cumulative steepness spectrum
!     IF IWCAP = 7; Calculate with Alves and Banner (2003)
!     ------------------------------------------------------------
!     For frequency dependent part of the spectrum
!       Calculate dissipation term due to whitecapping
!     ------------------------------------------------------------
!     For the whole frequency domain
!       Fill the matrices and PLWCAP-array
!     ------------------------------------------------------------
!     End of SWCAP
!     ------------------------------------------------------------
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   USE M_WCAP
!   USE ALL_VARS, ONLY : MT,AC2
   USE VARS_WAVE, ONLY : MT,AC2

   IMPLICIT NONE

   INTEGER, INTENT(IN)  :: ISSTOP, IDCMIN(MSC), IDCMAX(MSC)
   REAL, INTENT(IN)     :: DEP2(MT)
   REAL, INTENT(IN)     :: ETOT
   REAL, INTENT(IN)     :: KWAVE(MSC,MICMAX)                           
   REAL, INTENT(IN)     :: SPCDIR(MDC,6), SPCSIG(MSC)
   REAL, INTENT(OUT)    :: PLWCAP(MDC,MSC,NPTST)
   REAL, INTENT(IN OUT) :: IMATDA(MDC,MSC), IMATRA(MDC,MSC)
   REAL, INTENT(IN OUT) :: DISSIP(MDC,MSC), DISIMP(MDC,MSC)            
   REAL, INTENT(IN)     :: UFRIC                                       
   REAL, INTENT(IN)     :: CGO(MSC,MICMAX)                             

   INTEGER, SAVE     :: IENT = 0
   INTEGER           :: ID, IDDUM, IS, ID1, ID2, IF, IL, MXWCP,IG
   REAL              :: A, C_BJ, HM, HRMS, N1, N2
   REAL              :: QB_WC, SIG0, STP_OV, STP_PM
   REAL              :: DDIF, DELTA, DSTEEP, EBIN, XFAC          
   REAL              :: CPOW, CTOT, GAMMA                              
   REAL              :: BINSIZE                                      
   REAL, ALLOCATABLE :: C_K(:), C_LH(:), WCAP(:), WCIMPL(:)
   REAL, ALLOCATABLE :: CUMSTP(:,:), DSIGMA(:)                         
   REAL, ALLOCATABLE :: FCOS(:)                                        
   REAL              :: B, P                                           
   REAL              :: EF(MSC)                                        
   CHARACTER*20 NUMSTR, CHARS                                          
   CHARACTER*80 MSGSTR                                                 

   MXWCP = 7
   IF(IWCAP > MXWCP)THEN                                            
!
! Error message
!
     CHARS = NUMSTR(MXWCP+1,RNAN,'(I1)')
     CALL TXPBLA(CHARS,IF,IL)
     MSGSTR = 'Value for IWCAP should be less than '// CHARS(IF:IL)
     CALL MSGERR ( 4, MSGSTR )
     RETURN
   END IF
!
! Initialisation
!
   IF(ETOT    <= 0.) RETURN
   IF(ETOT2   <= 0.) RETURN
   IF(ETOT4   <= 0.) RETURN
   IF(ACTOT   <= 0.) RETURN
   IF(EDRKTOT <= 0.) RETURN
!
   ALLOCATE (C_K(MSC), C_LH(MSC), WCAP(MSC), WCIMPL(MSC))
   ALLOCATE (CUMSTP(0:MSC,MDC), DSIGMA(1:MSC))                         
   ALLOCATE (FCOS(1:MDC))                                              
   WCIMPL(1:MSC) = 0.
!
! Calculate coefficients
!
   IF((IWCAP == 1) .OR. (IWCAP == 2) .OR. (IWCAP == 5))THEN                 
!
! Calculate C_K
!
     STP_OV = KM_WAM * SQRT(ETOT)
     STP_PM = SQRT(PWCAP(2))
     N1     = PWCAP(11)
     N2     = 2. * PWCAP(9)
     C_K(:) = PWCAP(1) * (1. - PWCAP(10) +               &
              PWCAP(10) * (KWAVE(:,1) / KM_WAM)**N1) *   &
	      (STP_OV / STP_PM)**N2
!
   ENDIF
!
   IF((IWCAP == 4) .OR. (IWCAP == 5))THEN
!
! Calculate values for Hm and Qb
!
     HRMS   = SQRT(8. * ETOT)
     IF (IWCAP.EQ.4) HM = PWCAP(6) / KM01
     IF (IWCAP.EQ.5) HM = PWCAP(6) / (PWCAP(8) * KM_WAM)
     CALL FRABRE(HM, ETOT, QB_WC)
!
! Calculate C_BJ
!
     IF(HRMS >= HM)THEN
       C_BJ = PWCAP(7)  /  PI_W
     ELSE IF(HRMS > 0.)THEN
       C_BJ = (PWCAP(7) *  HM**2 * QB_WC) / (PI_W * HRMS**2)
     ELSE
       C_BJ = 0.
     END IF
   ENDIF
!
   IF(IWCAP == 3)THEN
!
! Calculate C_LH
!
     SIG0 = SQRT(ETOT2 / ETOT)
!
!       A = -(1./8.)*(ETOT2**2/(ETOT*ETOT4))*(GRAV_W**2/(ETOT*SIG0**4))
!       rewrite to prevent underflow
!
     A = -(1./8.) * GRAV_W**2 / ETOT4
     DO IS=1, ISSTOP
!          C_LH(IS) = PWCAP(5) * SQRT((ETOT * SIG0**4) / GRAV_W**2) *
!     &               EXP(A) * SIG0 * (SPCSIG(IS) / SIG0)**2
!          rewrite to prevent underflow:
!
       C_LH(IS) = PWCAP(5) * EXP(A) * SQRT(ETOT2) * SPCSIG(IS)**2 / GRAV_W
     END DO
   END IF
!
! In case of cumulative steepness method, calculate its spectrum          
!
   IF(IWCAP == 6)THEN                                              
     XFAC      = (SPCSIG(3)-SPCSIG(1))/(2.*SPCSIG(2))                 
     DSIGMA(:) = XFAC*SPCSIG(:)                                       
     DELTA     = SPCDIR(2,1)-SPCDIR(1,1)                              
!
!    --- precompute cos dependence                                    
!                                                                     
     DO ID1 = 1, MDC                                                  
       DDIF = REAL(ID1-1)*DELTA                                       
       FCOS(ID1) = (ABS(COS(DDIF)))**PWCAP(13)                        
     END DO                                                           
!
!    --- compute normalisation coefficient                            
!
     CPOW = PWCAP(13)                                                 
     IF(CPOW > 10)THEN                                              
       CTOT = SQRT(CPOW/(2.*PI_W))*(1.+0.25/CPOW)                       
     ELSE                                                             
       CTOT = 1./SQRT(PI_W)*GAMMA(0.5*CPOW+1.)/GAMMA(0.5*CPOW+0.5)      
     END IF                                                           
!
     CUMSTP = 0.                                                      
     DO IS = 1,ISSTOP                                                 
       BINSIZE = SPCSIG(IS)*DELTA*DSIGMA(IS)                         
       DO ID1 = 1,MDC                                                
         CUMSTP(IS,ID1) = CUMSTP(IS-1,ID1)                          
         DO ID2 = 1,MDC                                             
           EBIN= AC2(ID2,IS,IG)*BINSIZE                      
!           EBIN= AC2(ID2,IS,kcgrd(1))*BINSIZE                      
           DSTEEP = KWAVE(IS,1)**2*EBIN*FCOS(ABS(ID1-ID2)+1)       
           CUMSTP(IS,ID1) = CUMSTP(IS,ID1) + DSTEEP                
         END DO                                                     
       END DO                                                        
     END DO                                                           
     CUMSTP = PWCAP(12)*CUMSTP                                        
     CUMSTP = CUMSTP*CTOT                                             
   END IF                                                              
!
! Calculate dissipation according to Alves & Banner (2003)                
!
   IF(IWCAP == 7)THEN                                              
!                                                                         
!  Loop to calculate B(k)                                                 
!                                                                         
     DO IS = 1, ISSTOP                                                 
!                                                                         
!  Calculate E(f)                                                         
!                                                                         
       EF(IS) = 0.                                                     
       DO ID = 1,MDC                                                   
         EF(IS) = EF(IS) + AC2(ID,IS,IG)*SPCSIG(IS)*PI2_W*DDIR     
       END DO                                                          
!                                                                      
!  Calculate saturation spectrum B(k) from E(f)                          
!                                                                        
       B = (1./PI2_W) * CGO(IS,1) * KWAVE(IS,1)**3 * EF(IS)             
!                                                                        
!  Calculate exponent P of the relative saturation (B/Br)                
!                                                                        
       PWCAP(10)= 3. + TANH(25.76*(UFRIC*KWAVE(IS,1)/SPCSIG(IS)-0.1)) 
       P = 0.5*PWCAP(10)*(1. + TANH( 10.*( (B/PWCAP(12))**0.5 - 1.))) 
!                                                                        
!  Calculate WCAP(IS) from B(k) and P                                    
!                                                                        
       STP_OV = KM_WAM * SQRT(ETOT)                                   
       WCAP(IS) = PWCAP(1)*(B/PWCAP(12))**(P/2.) *                     &
                  STP_OV**PWCAP(9) * (KWAVE(IS,1)/KM_WAM)**PWCAP(11) * &
		  (GRAV_W**(0.5)*KWAVE(IS,1)**(0.5)/SPCSIG(IS))**(PWCAP(10)/2-1) * &
		  GRAV_W**(0.5)*KWAVE(IS,1)**(0.5)                                 
     END DO                                                           
   END IF                                                             
!
   IF(IWCAP < 6)THEN                                             
!
! Calculate the whitecapping source term WCAP(IS)
!
     DO IS=1, ISSTOP
       IF((IWCAP == 1) .OR. (IWCAP == 2) .OR.                         &
          ((IWCAP == 5) .AND. (C_BJ <= C_K(IS))))THEN
         WCAP(IS) = C_K(IS) * SIGM_10 * (KWAVE(IS,1) / KM_WAM)
       ELSE IF(IWCAP == 3)THEN
         WCAP(IS) = C_LH(IS)
       ELSE IF((IWCAP == 4) .OR. ((IWCAP == 5) .AND. (C_BJ >= C_K(IS))))THEN
         IF(IWCAP == 4) WCAP(IS) = C_BJ*SIGM01 *(KWAVE(IS,1)/KM01  )
         IF(IWCAP == 5) WCAP(IS) = C_BJ*SIGM_10*(KWAVE(IS,1)/KM_WAM)
!
! Calculate a term that is added to both sides of the equation to compensate
! for the strong non-linearity in the fraction of breaking waves Qb
!
         IF(HRMS < HM)THEN
           WCIMPL(IS)=WCAP(IS) * ((1.-QB_WC)/((HRMS**2/HM**2)-QB_WC))
           WCAP(IS)  =WCAP(IS) + WCIMPL(IS)
         END IF
       ELSE
         CALL MSGERR(2,'Whitecapping is inactive')
         WRITE (PRINTF,*) 'Occurs in gridpoint: ', IG
       END IF
     END DO

   END IF
!
! Fill the diagonal of the matrix and the PLWCAP-array
!
   IF(IWCAP /= 6)THEN                                              

     DO IS=1, ISSTOP
!
!      Only fill the values for the current sweep
!
       DO IDDUM = IDCMIN(IS), IDCMAX(IS)
         ID = MOD(IDDUM - 1 + MDC, MDC) + 1
         IMATDA(ID,IS) = IMATDA(ID,IS) + WCAP(IS)
         IMATRA(ID,IS) = IMATRA(ID,IS) - WCAP(IS) * AC2(ID,IS,IG)
	 
         DISSIP(ID,IS) = DISSIP(ID,IS) + WCAP(IS)                      
         IF (TESTFL) PLWCAP(ID,IS,IPTST) = -1.*(WCAP(IS)-WCIMPL(IS))
       END DO
     END DO

   ELSE

     DO IS=1, ISSTOP                                                   
!                                                                         
!      Only fill the values for the current sweep                      
!                                                                         
       DO IDDUM = IDCMIN(IS), IDCMAX(IS)                               
         ID = MOD(IDDUM - 1 + MDC, MDC) + 1                            
         IMATDA(ID,IS) = IMATDA(ID,IS) + CUMSTP(IS,ID)                 
         DISSIP(ID,IS) = DISSIP(ID,IS) + CUMSTP(IS,ID)                 
         IF (TESTFL) PLWCAP(ID,IS,IPTST) = -1.*CUMSTP(IS,ID)           
       END DO                                                          
     END DO                                                            

   END IF
!
! Add the implicit part to the right-hand side
!
   IF((IWCAP == 4) .OR. (IWCAP == 5))THEN
     DO IS=1, ISSTOP
!
!    Only fill the values for the current sweep
!
       DO IDDUM = IDCMIN(IS), IDCMAX(IS)
         ID = MOD(IDDUM - 1 + MDC, MDC) + 1
         IMATRA(ID,IS) = IMATRA(ID,IS) + WCIMPL(IS) * AC2(ID,IS,IG)
         DISIMP(ID,IS) = DISIMP(ID,IS) + WCIMPL(IS) * AC2(ID,IS,IG)      
       END DO
     END DO
   END IF
!
   DEALLOCATE (C_K, C_LH, WCAP, WCIMPL)
   DEALLOCATE (CUMSTP, DSIGMA, FCOS)                                   
!
   RETURN
   END SUBROUTINE SWCAP

!****************************************************************
!
!!$  SUBROUTINE BRKPAR(BRCOEF,ECOS,ESIN,AC2,SPCSIG,DEP2,RDX,RDY,IG)     
  SUBROUTINE BRKPAR(BRCOEF,ECOS,ESIN,AC2,SPCSIG,DEP2,IG)     
  
! (This subroutine has not been fully tested yet)    
!
!****************************************************************
!
      USE SWCOMM3                                                         
      USE SWCOMM4                                                         
      USE OCPCOMM4 
      USE MOD_PREC
      USE ALL_VARS, ONLY : MT,NTSN,NBSN,VX,VY,ART2 
#     if defined (SPHERICAL)
      USE MOD_SPHERICAL, ONLY : TPI,DEG2RAD
#     endif                                                            
!
      IMPLICIT NONE                                                       
!
!  2. Purpose
!
!     Determine the bottom slope in upwave direction and calculate
!     the slope dependent breaking parameter according to Nelson (1987)
!     Note that Nelson (1987) is used here since in Nelson (1994a,1994b)
!     an error is present in the equation.
!
!  3. Method
!
!     The breaker parameter is given by:
!
!     Hm / d =  0.55 + 0.88 exp ( -0.012 * cot beta)
!
!     with beta the angle the bed makes with the horizontal. This
!     above equation is only valid for positive slopes (Negative
!     slopes were not considered by Nelson. For very steep slopes
!     (>0.05 say) a very large breaker parameter is obtained (>>1).
!
!     To ensure wave breaking in laboratory cases (with very steep
!     slopes an upper limit of 0.81 (which corresponds to a bottom
!     slope of 0.01) is imposed on the model of Nelson.
!
!     For negative bottom slopes (not considered by Nelson) a value
!     op 0.73 is imposed (which is the average value in Table 2 of
!     Battjes and Janssen (1978).
!
      REAL, INTENT(OUT) :: BRCOEF    ! variable breaker coefficient       

      REAL, INTENT(IN)  :: SPCSIG(MSC)                                    
      REAL, INTENT(IN)  :: AC2(MDC,MSC,0:MT)  ! action densities         
      REAL, INTENT(IN)  :: ECOS(MDC), ESIN(MDC)  ! Cos and Sin of Theta   
      REAL, INTENT(IN)  :: DEP2(MT)         ! depths at grid points    

!!$      REAL, INTENT(IN)  :: RDX(10), RDY(10)
      INTEGER, INTENT(IN) :: IG                               

!     9. STRUCTURE
!
!     ------------------------------------------------------------
!     Calculate total energy per direction for all frequencies
!       determine action density per direction weighted with cos/sin
!       determine mean propagation direction of energy
!     ------------------------------------------------------------
!     calculate the depth derivative in the mean wave direction
!      according to dd/ds (see also subroutine SPROSD)
!     calculate the slope dependend breaking coefficient
!     ------------------------------------------------------------
!     End of NELSON
!     -------------------------------------------------------------
!
!************************************************************************
!
  INTEGER :: ID    ,IS      ! counters 
  INTEGER :: J,I1,I2
  REAL(SP) :: F1    
# if defined (SPHERICAL)
  REAL(DP) XTMP,XTMP1	    
# endif
!
!
  REAL  :: ETOTS,EEX,EEY,EAD,SIGMA1,COSDIR,SINDIR,DDDX,DDDY,DDDS,DETOT                                
!
!     *** determine the average wave direction ***
!
  EEX   = 0.
  EEY   = 0.
  ETOTS = 0.
  DO ID = 1, MDC
    EAD = 0.
    DO IS = 1, MSC
      SIGMA1 = SPCSIG(IS)                                            
      DETOT  = SIGMA1**2 * AC2(ID,IS,KCGRD(1))
      EAD    = EAD + DETOT
    END DO
    ETOTS = ETOTS + EAD
    EEX   = EEX + EAD * ECOS(ID)
    EEY   = EEY + EAD * ESIN(ID)
  END DO
!
  IF(ETOTS > 0.)THEN
    COSDIR = EEX / ETOTS
    SINDIR = EEY / ETOTS
  ELSE
    COSDIR = 1.
    SINDIR = 0.
  END IF
!
! *** Determine bottom slope in average wave propagation direction ***
!
!  DDDX =  RDX(1) * (DEP2(KCGRD(1)) - DEP2(KCGRD(2)))              &
!        + RDX(2) * (DEP2(KCGRD(1)) - DEP2(KCGRD(3)))
!  DDDY =  RDY(1) * (DEP2(KCGRD(1)) - DEP2(KCGRD(2)))              &
!        + RDY(2) * (DEP2(KCGRD(1)) - DEP2(KCGRD(3)))
  DDDX = 0.0_SP
  DDDY = 0.0_SP
  
  DO J=1,NTSN(IG)-1
    I1 = NBSN(IG,J)
    I2 = NBSN(IG,J+1)
    F1 = 0.50_SP*(DEP2(I1)+DEP2(I2)) 
#   if defined (SPHERICAL)
    DDDX=DDDX+F1*(VY(I1)-VY(I2))*TPI

    XTMP  = VX(I2)*TPI-VX(I1)*TPI  
    XTMP1 = VX(I2)-VX(I1)  
    IF(XTMP1 >  180.0_SP)THEN
      XTMP = -360.0_SP*TPI+XTMP     
    ELSE IF(XTMP1 < -180.0_SP)THEN
      XTMP =  360.0_SP*TPI+XTMP     
    END IF
    DDDY=DDDY+F1*XTMP*COS(DEG2RAD*VY(IG))
#   else
    DDDX = DDDX + F1*(VY(I1)-VY(I2))
    DDDY = DDDY + F1*(VX(I2)-VX(I1))
#   endif
  END DO
  DDDX = DDDX/ART2(IG)
  DDDY = DDDY/ART2(IG)

!
  DDDS = -1. * ( DDDX * COSDIR + DDDY * SINDIR )
!
! *** calculate breaking coefficient according to Nelson (1987) ***
!
  IF(DDDS >= 0.)THEN
    DDDS   = MAX ( 1.E-6 , DDDS)
    BRCOEF = PSURF(4) + PSURF(7) * EXP ( -PSURF(8) / DDDS )           
    BRCOEF = MIN ( PSURF(5) , BRCOEF )                                
  ELSE
    BRCOEF = PSURF(6)                                                 
  END IF
!
!  PSURF(2) = BRKVAR                              deleted             
!  PSURF(2) is no longer used to transmit br. coefficient             

  RETURN
  END SUBROUTINE BRKPAR
