!
!****************************************************************
!
  SUBROUTINE SWAPAR(IG,NICMAX,DEP2,KWAVE,CGO)
!
!****************************************************************
!
!     computes the wave parameters K and CGO in the nearby
!     points, depending of the sweep direction.
!     The nearby points are indicated with the index IC (see
!     FUNCTION ICODE(_,_)
!
!  Method
!
!     The wave number K(IS,iC) is computed with the dispersion relation:
!
!     S = GRAV K(IS,IC)tanh(K(IS,IC)DEP(IX,IY))
!
!     where S = is logarithmic distributed via LOGSIG
!
!     The group velocity CGO in the case without current is equal to
!
!                    1       K(IS,IC)DEP(IX,IY)          S
!     CGO(IS,IC) = ( - + --------------------------) -----------
!                    2   2 sinh 2K(IS,IC)DEP(IX,IY)  |k(IS,IC)|
!
!
  USE M_GENARR
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE OCPCOMM4
# if defined (EXPLICIT)
  USE MOD_ACTION_EX
# else         
  USE MOD_ACTION_IM
# endif    
  USE ALL_VARS 
   
  IMPLICIT NONE                                                     

  INTEGER :: IC,IS,ID,IG,NICMAX,INDX
  REAL    :: NN(1:MSC), ND(1:MSC)                                    
  REAL    :: DEP2(MT),KWAVE(MSC,NICMAX),CGO(MSC,NICMAX)
  REAL    :: DEPLOC

  DO IC = 1, NICMAX
    IF(IC == 1)THEN
      INDX  = IG
      DEPLOC = DEP2(INDX)
      IF(DEPLOC <= DEPMIN)THEN
!     *** depth is negative ***
        DO IS = 1, MSC
          KWAVE(IS,IC) = -1.                                           
          CGO(IS,IC)   = 0.                                            
        END DO	 
      ELSE
!     *** call KSCIP1 to compute KWAVE and CGO ***
        CALL KSCIP1(MSC,SPCSIG,DEPLOC,KWAVE(1,IC),CGO(1,IC),NN,ND)                                 
      ENDIF
    ELSE                                                    
      INDX  = NBVE(IG,IC-1)
      DEPLOC = DEP2(NV(INDX,1))+DEP2(NV(INDX,2))+DEP2(NV(INDX,3))
      DEPLOC = DEPLOC/3.0
      IF(DEPLOC <= DEPMIN)THEN
!     *** depth is negative ***
        DO IS = 1, MSC
          KWAVE(IS,IC) = -1.                                           
          CGO(IS,IC)   = 0.                                            
        END DO	 
      ELSE
!     *** call KSCIP1 to compute KWAVE and CGO ***
        CALL KSCIP1(MSC,SPCSIG,DEPLOC,KWAVE(1,IC),CGO(1,IC),NN,ND)
      ENDIF
    END IF  
  ENDDO                                                              

  RETURN
  END SUBROUTINE SWAPAR
 
 !
!
!****************************************************************
!
  SUBROUTINE SWAPAR1(I,IS,ID,DEP2,KWAVEL,CGOL)
!
!****************************************************************
!
!     computes the wave parameters K and CGO in the nearby
!     points, depending of the sweep direction.
!     The nearby points are indicated with the index IC (see
!     FUNCTION ICODE(_,_)
!
!  Method
!
!     The wave number K(IS,iC) is computed with the dispersion relation:
!
!     S = GRAV K(IS,IC)tanh(K(IS,IC)DEP(IX,IY))
!
!     where S = is logarithmic distributed via LOGSIG
!
!     The group velocity CGO in the case without current is equal to
!
!                    1       K(IS,IC)DEP(IX,IY)          S
!     CGO(IS,IC) = ( - + --------------------------) -----------
!                    2   2 sinh 2K(IS,IC)DEP(IX,IY)  |k(IS,IC)|
!
  USE M_GENARR
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE OCPCOMM4 
# if defined (EXPLICIT)
  USE MOD_ACTION_EX
# else        
  USE MOD_ACTION_IM
# endif    
  USE ALL_VARS, ONLY : NV,NTVE,MT
   
  IMPLICIT NONE                                                     
!
  INTEGER :: IC,IS,ID,I,INDX
  REAL    :: NN(1:MSC), ND(1:MSC)                                    
  REAL    :: DEP2(MT),KWAVEL,CGOL
  REAL    :: DEPLOC,SPCSIGL

  DEPLOC = (DEP2(NV(I,1))+DEP2(NV(I,2))+DEP2(NV(I,3)))/3.0
  IF(DEPLOC <= DEPMIN)THEN
!   *** depth is negative ***
    KWAVEL = -1.                                           
    CGOL   = 0.                                            
  ELSE
!   *** call KSCIP1 to compute KWAVE and CGO ***
    SPCSIGL = SPCSIG(IS)
    CALL KSCIP1(1,SPCSIGL,DEPLOC,KWAVEL,CGOL,NN,ND)                                 
  END IF

  RETURN
  END SUBROUTINE SWAPAR1
!    
!****************************************************************
!
   SUBROUTINE SPROXY (I1     ,IS     ,ID     ,CAXL   ,CAYL   ,  &
                      CG0L   ,ECOSL  ,ESINL  ,UX2L   ,UY2L   )
!
!****************************************************************
!
!        computes the propagation velocities of energy in X-, Y-
!        -space, i.e., CAX, CAY, in the presence or absence of
!        currents, for the action balance equation.
!
!        The propagation velocities are computed for the fully 360
!        degrees sector.
!
!     METHOD
!
!        The next equation are calculated:
!
!              @X     _
!        CAX = -- = n C cos (id) + Ux  = CGO cos(id) + Ux
!              @T
!
!              @Y     _
!        CAY = -- = n C sin(id)  + Uy  = CGO sin(id) + Uy
!              @T
!                                                         _
!
!       ******************************************************************
!       *  attention! in the action balance equation the term            *
!       *  dx                                                            *
!       *  -- = CGO + U = CX  with x, CGO, U and CX vectors              *
!       *  dt                                                            *
!       *  is in the literature the term dx/dt often indicated           *
!       *  with CX and CY in the action balance equation.                *
!       *  In this program we use:    CAX = CGO + U                      *
!       ******************************************************************
!
!   ------------------------------------------------------------
!   If depth is negative ( DEP(IX,IY) <= 0), then,
!     For every point in S and D-direction do,
!       Give propagation velocities default values :
!       CAX(ID,IS,IC)     = 0.   {propagation velocity of energy in X-dir.}
!       CAY(ID,IS,IC)     = 0.   {propagation velocity of energy in Y-dir.}
!     ---------------------------------------------------------
!   Else if current is on (ICUR > 0) then,
!     For every point in S and D-direction do,  {using the output of SWAPAR}
!       S = logaritmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S cos(D)
!       CAX = ( - + ------------------------) --------- + UX2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S sin(D)
!       CAY = ( - + ------------------------) -------- + UY2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!   Else if current is not on (ICUR = 0)
!     For every point in S and D-direction do
!       S = logarithmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S cos(D)
!       CAX = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S sin(D)
!       CAY = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!     ----------------------------------------------------------
!   End IF
!   ------------------------------------------------------------
!****************************************************************
   USE SWCOMM3
   USE SWCOMM4
   USE OCPCOMM4
   USE M_DIFFR
#  if defined (EXPLICIT)
   USE MOD_ACTION_EX
#  else      
   USE MOD_ACTION_IM
#  endif   

   IMPLICIT NONE

   INTEGER  IC,IS,ID,I1,NICMAX
   REAL     CAXL,CAYL,CG0L,ECOSL,ESINL,UX2L,UY2L

   CAXL = CG0L * ECOSL
   CAYL = CG0L * ESINL
!
!    --- adapt the velocities in case of diffraction
!
   IF(IDIFFR == 1 .AND. PDIFFR(3) /= 0.)THEN
!     CAXL = CAXL*DIFPARAM(I1)      
!     CAYL = CAYL*DIFPARAM(I1)      
   END IF
!
!    --- ambient currents added
!
   IF(ICUR == 1)THEN
     CAXL = CAXL + UX2L
     CAYL = CAYL + UY2L
   END IF

   RETURN
   END SUBROUTINE SPROXY
!
!

!****************************************************************
!
   SUBROUTINE SPROXY2 (CAXL   ,CAYL   ,  &
                      CG0L   ,ECOSL  ,ESINL  ,UX2L   ,UY2L   )
!
!****************************************************************
!
!        computes the propagation velocities of energy in X-, Y-
!        -space, i.e., CAX, CAY, in the presence or absence of
!        currents, for the action balance equation.
!
!        The propagation velocities are computed for the fully 360
!        degrees sector.
!
!     METHOD
!
!        The next equation are calculated:
!
!              @X     _
!        CAX = -- = n C cos (id) + Ux  = CGO cos(id) + Ux
!              @T
!
!              @Y     _
!        CAY = -- = n C sin(id)  + Uy  = CGO sin(id) + Uy
!              @T
!                                                         _
!
!       ******************************************************************
!       *  attention! in the action balance equation the term            *
!       *  dx                                                            *
!       *  -- = CGO + U = CX  with x, CGO, U and CX vectors              *
!       *  dt                                                            *
!       *  is in the literature the term dx/dt often indicated           *
!       *  with CX and CY in the action balance equation.                *
!       *  In this program we use:    CAX = CGO + U                      *
!       ******************************************************************
!
!   ------------------------------------------------------------
!   If depth is negative ( DEP(IX,IY) <= 0), then,
!     For every point in S and D-direction do,
!       Give propagation velocities default values :
!       CAX(ID,IS,IC)     = 0.   {propagation velocity of energy in X-dir.}
!       CAY(ID,IS,IC)     = 0.   {propagation velocity of energy in Y-dir.}
!     ---------------------------------------------------------
!   Else if current is on (ICUR > 0) then,
!     For every point in S and D-direction do,  {using the output of SWAPAR}
!       S = logaritmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S cos(D)
!       CAX = ( - + ------------------------) --------- + UX2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S sin(D)
!       CAY = ( - + ------------------------) -------- + UY2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!   Else if current is not on (ICUR = 0)
!     For every point in S and D-direction do
!       S = logarithmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S cos(D)
!       CAX = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S sin(D)
!       CAY = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!     ----------------------------------------------------------
!   End IF
!   ------------------------------------------------------------
!****************************************************************
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   USE M_DIFFR 
#  if defined (EXPLICIT)
   USE MOD_ACTION_EX
#  else      
   USE MOD_ACTION_IM
#  endif   
   
   IMPLICIT NONE                                              

   INTEGER  IC,IS,ID,I1,NICMAX
   REAL     CAXL(MDC,MSC),CAYL(MDC,MSC),CG0L(MDC,MSC),ECOSL(MDC),ESINL(MDC),UX2L,UY2L

   DO ID=1,MDC
     CAXL(ID,:) = CG0L(ID,:) * ECOSL(ID)
     CAYL(ID,:) = CG0L(ID,:) * ESINL(ID)
   END DO

!
!    --- adapt the velocities in case of diffraction
!
   IF(IDIFFR == 1 .AND. PDIFFR(3) /= 0.)THEN 
!     CAXL = CAXL*DIFPARAM(I1)      
!     CAYL = CAYL*DIFPARAM(I1)      
   END IF
!
!    --- ambient currents added
!
   IF(ICUR == 1)THEN 
     CAXL = CAXL + UX2L
     CAYL = CAYL + UY2L
   END IF

   RETURN
   END SUBROUTINE SPROXY2
!
!
!****************************************************************
!
   SUBROUTINE SPROXY3 (CAXL   ,CAYLA ,CAYLB,  &          !yzhang_w3
              CG0L   ,ECOSL  ,ESINL  ,UX2L   ,UY2L, DLTYETMPP, DLTXETMPP ,DLTXEA , DLTXEB)
!
!****************************************************************
!
!        computes the propagation velocities of energy in X-, Y-
!        -space, i.e., CAX, CAY, in the presence or absence of
!        currents, for the action balance equation.
!
!        The propagation velocities are computed for the fully 360
!        degrees sector.
!
!     METHOD
!
!        The next equation are calculated:
!
!              @X     _
!        CAX = -- = n C cos (id) + Ux  = CGO cos(id) + Ux
!              @T
!
!              @Y     _
!        CAY = -- = n C sin(id)  + Uy  = CGO sin(id) + Uy
!              @T
!                                                         _
!
!       ******************************************************************
!       *  attention! in the action balance equation the term            *
!       *  dx                                                            *
!       *  -- = CGO + U = CX  with x, CGO, U and CX vectors              *
!       *  dt                                                            *
!       *  is in the literature the term dx/dt often indicated           *
!       *  with CX and CY in the action balance equation.                *
!       *  In this program we use:    CAX = CGO + U                      *
!       ******************************************************************
!
!   ------------------------------------------------------------
!   If depth is negative ( DEP(IX,IY) <= 0), then,
!     For every point in S and D-direction do,
!       Give propagation velocities default values :
!       CAX(ID,IS,IC)     = 0.   {propagation velocity of energy in X-dir.}
!       CAY(ID,IS,IC)     = 0.   {propagation velocity of energy in Y-dir.}
!     ---------------------------------------------------------
!   Else if current is on (ICUR > 0) then,
!     For every point in S and D-direction do,  {using the output of SWAPAR}
!       S = logaritmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S cos(D)
!       CAX = ( - + ------------------------) --------- + UX2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S sin(D)
!       CAY = ( - + ------------------------) -------- + UY2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!   Else if current is not on (ICUR = 0)
!     For every point in S and D-direction do
!       S = logarithmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S cos(D)
!       CAX = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S sin(D)
!       CAY = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!     ----------------------------------------------------------
!   End IF
!   ------------------------------------------------------------
!****************************************************************
   USE SWCOMM3
   USE SWCOMM4
   USE OCPCOMM4
   USE M_DIFFR
   USE MOD_ACTION_EX

   IMPLICIT NONE
   INTEGER  IC,IS,ID,I1,NICMAX
   REAL   CAXL(MDC,MSC),CAYLA(MDC,MSC),CAYLB(MDC,MSC),CG0L(MDC,MSC),ECOSL(MDC),ESINL(MDC),UX2L,UY2L
   REAL     DLTXETMPP,DLTXEA,DLTXEB,DLTYETMPP

   DO ID=1,MDC
     CAXL(ID,:) = CG0L(ID,:) * ECOSL(ID) * DLTYETMPP
     CAYLA(ID,:) = CG0L(ID,:) * ESINL(ID) * DLTXEA
     CAYLB(ID,:) = CG0L(ID,:) * ESINL(ID) * DLTXEB
   END DO

!
!    --- adapt the velocities in case of diffraction
!
   IF(IDIFFR == 1 .AND. PDIFFR(3) /= 0.)THEN
!     CAXL = CAXL*DIFPARAM(I1)      
!     CAYL = CAYL*DIFPARAM(I1)      
   END IF
!
!    --- ambient currents added
!
   IF(ICUR == 1)THEN
     CAXL = CAXL + UX2L*DLTYETMPP
     CAYLA = CAYLA + UY2L*DLTXETMPP
     CAYLB = CAYLB + UY2L*DLTXETMPP
   END IF

   RETURN
   END SUBROUTINE SPROXY3
!
!
!****************************************************************

!****************************************************************
!
  SUBROUTINE SPROSD (KWAVE    ,CAS      ,CAD      ,          &
                     CGO      ,DEP2     ,DEP1     ,          &
		     ECOS     ,ESIN     ,UX2      ,          &
		     UY2      ,IDCMIN   ,IDCMAX   ,          &
		     COSCOS   ,SINSIN   ,SINCOS   ,          &
		     RDX      ,RDY      ,CAX      ,          &
		     CAY      ,IG       )
!
!****************************************************************
!
!     computes the propagation velocities of energy in S- and
!     D-space, i.e., CAS, CAD, in the presence or absence of
!     currents, for the action balance equation.
!
!  Method
!
!     The next equation are solved numerically
!
!           @S   @S   @D   _     @D   @D          _   @U
!     CAS = -- = -- [ -- + U . ( -- + --) ] - CGO K . --
!           @T   @D   @T         @X   @Y              @s
!
!           with:   @S       KS
!                   -- =  ---------
!                   @D    sinh(2KD)
!
!           @D      S      @D         @D           @Ux   @Uy
!     CAD = -- = ------- [ --sin(D) - --cos(D)] + [--- - ---] *
!           @T  sinh(2KD)  @X         @Y            @X   @Y
!
!                        @Uy               @Ux
!     * sin(D)cos(D) +   ---sin(D)sin(D) - ---cos(D)cos(D)
!                        @X                @Y
!   ------------------------------------------------------------
!   For current sweep and two adjacent sweeps do
!       determine interpolation factors RDXL and RDYL
!       determine depth and current gradients
!   ------------------------------------------------------------
!   For each frequency do
!       determine auxiliary quantities depending on sigma
!       For each direction in the sweep and two neighbouring
!           directions do
!           If IREFR=-1
!           Then compute reduction factor for contribution due
!                to depth gradient
!           ----------------------------------------------------
!           determine sweep in which this direction is located
!           using gradients of the proper sweep determine
!           Csigma (CAS) and Ctheta (CAD)
!   ------------------------------------------------------------
!   If ITFRE=0
!   Then make values of CAS=0
!   ------------------------------------------------------------
!   If IREFR=0
!   Then make values of CAD=0
!   ------------------------------------------------------------
!
!****************************************************************
!
  USE M_GENARR
  USE SWCOMM2                                           
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE TIMECOMM                                                        
  USE OCPCOMM4                                                        
  USE M_DIFFR   
# if defined (EXPLICIT)
  USE MOD_ACTION_EX
# else            
  USE MOD_ACTION_IM
# endif      
  USE ALL_VARS, ONLY : NTVE,NBVE,NV,VX,VY,XC,YC,ART1,MT,ISONB_W,NBSN,NTSN  
# if defined (SPHERICAL)
  USE CONTROL, ONLY : TPI,DEG2RAD                                                      
# endif      
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: IDCMIN(MSC), IDCMAX(MSC)
  REAL  :: CAS(MDC,MSC,MICMAX)                                        
  REAL  :: CAD(MDC,MSC,MICMAX)                                        
  REAL  :: CAX(MDC,MSC,MICMAX)                                        
  REAL  :: CAY(MDC,MSC,MICMAX)                                        
  REAL  :: CGO(MSC,MICMAX)                                            
  REAL  :: DEP2(MT) ,DEP1(MT) ,ECOS(MDC) ,ESIN(MDC) ,COSCOS(MDC) ,    &
	   SINSIN(MDC) ,SINCOS(MDC)
  REAL  :: KWAVE(MSC,MICMAX)                                          
  REAL  :: UX2(MT) ,UY2(MT) ,RDX(10) ,RDY(10)  
  INTEGER  IENT ,IS ,ID ,II ,SWPNGB ,IDDUM ,ID1 ,ID2 ,ISWEEP   
  INTEGER :: ISWP                                              
  INTEGER :: IC, KCGI                                             
  LOGICAL :: VALSWP                                            
  REAL    :: VLSINH ,KD1   ,COEF
  REAL    :: RDXL(2) ,RDYL(2) ,DET ,DX2 ,DY2 ,DX3 ,DY3
  REAL    :: DPDX   ,DPDY   ,DUXDX ,DUXDY ,DUYDX ,DUYDY
  REAL    :: CAST1    ,CAST2    ,CAST3(3) ,CAST4(3) ,               &
             CAST5    ,CAST6(3) ,CAST7(3) ,CAST8(3) ,CAST9(3) ,     &
	     CADT1    ,CADT2(3) ,CADT3(3) ,                         &
	     CADT4(3) ,CADT5(3) ,CADT6(3) ,CADT7(3)
  REAL    :: DLOC1, DLOC2, DLOC3

  INTEGER, PARAMETER :: SWP_ARRAY(1:3) = (/2,1,3/)

  REAL    :: SUMHDX,SUMHDY,SUMUXDX,SUMUXDY,SUMUYDX,SUMUYDY
  REAL    :: SUMUXHDY,SUMUYHDX
  INTEGER :: IG,I
  REAL    :: UX,UY
  REAL    :: X1,X2,X3,Y1,Y2,Y3,DX1,DY1
# if defined (SPHERICAL)
  REAL    :: DX_TMP
# endif      
      
  CAST1 = 0.
  CAST2 = 0.
  CAST5 = 0.
  CADT1 = 0.
  DLOC1 = DEP2(IG)

  CAS(:,:,1) = 0.
  CAD(:,:,1) = 0.
!
! *** coefficients for CAS -> function of IX and IY only ***
!
  IF(NSTATC == 0 .OR. .NOT.DYNDEP)THEN                            
!   *** stationary calculation ***
    CAST2 = 0.
  ELSE
!   nonstationary depth, CAST2 is @D/@t
    CAST2 = (DLOC1-DEP1(IG))*RDTIM
  END IF
      
  SUMHDX = 0.0
  SUMHDY = 0.0
  IF(ICUR /= 0)THEN
    SUMUXDX = 0.0
    SUMUXDY = 0.0
    SUMUYDX = 0.0
    SUMUYDY = 0.0
    SUMUXHDY = 0.0
    SUMUYHDX = 0.0
  END IF
      
  DO I = 1,NTVE(IG)
    UX = UX2(NV(NBVE(IG,I),1))+UX2(NV(NBVE(IG,I),2))+UX2(NV(NBVE(IG,I),3))
    UX = UX/3.0     
    UY = UY2(NV(NBVE(IG,I),1))+UY2(NV(NBVE(IG,I),2))+UY2(NV(NBVE(IG,I),3))
    UY = UY/3.0     
	
    DLOC2 = DEP2(NV(NBVE(IG,I),1)) +         &
            DEP2(NV(NBVE(IG,I),2)) +         &
            DEP2(NV(NBVE(IG,I),3)) 
    DLOC2 = DLOC2/3.0

    IF(IREFR == -1)THEN                                          
!     limitation of depths in neighbouring grid points
      DLOC2 = MIN(DLOC2, PNUMS(17)*DLOC1)
    ELSE  
      IF(ABS(DLOC2 - DLOC1) > 100.0)THEN
        DLOC2 = DLOC1
      ELSE	                                                         
!       no limitation                                                
        DLOC2 = DLOC2
      END IF                                           
    END IF                                                         
 
    X1 = VX(IG)
    Y1 = VY(IG)
	
    IF(NV(NBVE(IG,I),1) == IG)THEN
      X2 = VX(NV(NBVE(IG,I),2))	
      Y2 = VY(NV(NBVE(IG,I),2))	
      X3 = VX(NV(NBVE(IG,I),3))	
      Y3 = VY(NV(NBVE(IG,I),3))	
    ELSE IF(NV(NBVE(IG,I),2) == IG)THEN
      X2 = VX(NV(NBVE(IG,I),3))	
      Y2 = VY(NV(NBVE(IG,I),3))	
      X3 = VX(NV(NBVE(IG,I),1))	
      Y3 = VY(NV(NBVE(IG,I),1))	
    ELSE  
      X2 = VX(NV(NBVE(IG,I),1))	
      Y2 = VY(NV(NBVE(IG,I),1))	
      X3 = VX(NV(NBVE(IG,I),2))	
      Y3 = VY(NV(NBVE(IG,I),2))	
    END IF

# if defined (SPHERICAL)
    if ((abs(X1-X2)).GT.300.0)then
      DX1 = 0.5*(X1+X2+360.0)
      if (DX1.GE.360)DX1=DX1-360.0
    else
    DX1 = 0.5*(X1+X2)
    endif
    DY1 = 0.5*(Y1+Y2)
# else
    DX1 = 0.5*(X1+X2)
    DY1 = 0.5*(Y1+Y2)
# endif

# if defined (SPHERICAL)
    if  (Y1.GE.89.99)then
      DX1=X2
    elseif (Y2.GE.89.99)then
      DX1=X1
    endif
# endif

    DX1 = XC(NBVE(IG,I))-DX1
    DY1 = YC(NBVE(IG,I))-DY1

# if defined (SPHERICAL)
    DX_TMP = DX1
    DX1 = DX_TMP*TPI
    IF(DX_TMP >  180.0_SP)THEN
      DX1 = -360.0_SP*TPI+DX1
    ELSE IF(DX_TMP < -180.0_SP)THEN
      DX1 =  360.0_SP*TPI+DX1
    END IF 
    DX1 = DX1*COS(DEG2RAD*VY(IG))
    DY1 = DY1*TPI 
# endif    
	

# if defined (SPHERICAL)
    if ((abs(X1-X3)).GT.300.0)then
      DX2 = 0.5*(X1+X3+360.0)
      if (DX2.GE.360)DX2=DX2-360.0
    else
    DX2 = 0.5*(X1+X3)
    endif
    DY2 = 0.5*(Y1+Y3)
# else
    DX2 = 0.5*(X1+X3)
    DY2 = 0.5*(Y1+Y3)
# endif

# if defined (SPHERICAL)
    if  (Y1.GE.89.99)then
      DX2=X3
    elseif (Y3.GE.89.99)then
      DX2=X1
    endif
# endif
    DX2 = DX2-XC(NBVE(IG,I))
    DY2 = DY2-YC(NBVE(IG,I))
# if defined (SPHERICAL)
    DX_TMP = DX2
    DX2 = DX_TMP*TPI
    IF(DX_TMP >  180.0_SP)THEN
      DX2 = -360.0_SP*TPI+DX2
    ELSE IF(DX_TMP < -180.0_SP)THEN
      DX2 =  360.0_SP*TPI+DX2
    END IF 
    DX2 = DX2*COS(DEG2RAD*VY(IG))
    DY2 = DY2*TPI 
# endif    
	
    DX = DX1+DX2
    DY = DY1+DY2
	
    SUMHDX = SUMHDX - DLOC2*DX
    SUMHDY = SUMHDY - DLOC2*DY
	
    IF(ICUR /= 0)THEN
      SUMUXDX = SUMUXDX - UX*DX
      SUMUXDY = SUMUXDY - UX*DY
      SUMUYDX = SUMUYDX - UY*DX
      SUMUYDY = SUMUYDY - UY*DY
      SUMUXHDY = SUMUXHDY - UX*DLOC2*DY
      SUMUYHDX = SUMUYHDX - UY*DLOC2*DX
    END IF
  END DO	    

  IF(ISONB_W(IG) == 1)THEN
    DX = VX(NBSN(IG,2))-VX(NBSN(IG,NTSN(IG)-1))
# if defined (SPHERICAL)
    DX_TMP = DX
    DX = DX_TMP*TPI
    IF(DX_TMP >  180.0_SP)THEN
      DX = -360.0_SP*TPI+DX
    ELSE IF(DX_TMP < -180.0_SP)THEN
      DX =  360.0_SP*TPI+DX
    END IF 
    DX = DX*COS(DEG2RAD*VY(IG))
# endif    
    DX = 0.5*DX
    DY = VY(NBSN(IG,2))-VY(NBSN(IG,NTSN(IG)-1))
# if defined (SPHERICAL)
    DY = DY*TPI 
# endif    
    DY = 0.5*DY
	
    SUMHDX = SUMHDX - DEP2(IG)*DX
    SUMHDY = SUMHDY - DEP2(IG)*DY

    IF(ICUR /= 0)THEN
      SUMUXDX = SUMUXDX - UX2(IG)*DX
      SUMUXDY = SUMUXDY - UX2(IG)*DY
      SUMUYDX = SUMUYDX - UY2(IG)*DX
      SUMUYDY = SUMUYDY - UY2(IG)*DY
      SUMUXHDY = SUMUXHDY - UX2(IG)*DEP2(IG)*DY
      SUMUYHDX = SUMUYHDX - UY2(IG)*DEP2(IG)*DX
    END IF
  END IF	

  DO IS =1,MSC
    KD1 = KWAVE(IS,1)*DLOC1
    IF(KD1 > 30.0)KD1 = 30.
    VLSINH = SINH(2.*KD1)
    COEF   = SPCSIG(IS)/VLSINH
!
!   *** coefficients for CAS -> function of IS only ***
!
    CAST1 = KWAVE(IS,1)*COEF
    CAST5 = CGO(IS,1)*KWAVE(IS,1)
!
!   *** coefficients for CAD -> function of IS only ***
!
    CADT1 =  COEF
!
!   loop over spectral directions
!
    DO IDDUM = IDCMIN(IS)-1, IDCMAX(IS)+1  
      ID = MOD(IDDUM-1+MDC, MDC)+1
      IF(ICUR == 0)THEN
        CAS(ID,IS,1) = CAST1*CAST2
        CAD(ID,IS,1) = CADT1*(ESIN(ID)*SUMHDY+ECOS(ID)*SUMHDX)
        CAD(ID,IS,1) = CAD(ID,IS,1)/ART1(IG)
	    
!       --- adapt the velocity in case of diffraction                
        IF(IDIFFR == 1)THEN                                        
!          CAD(ID,IS,1) = DIFPARAM(IG)*CAD(ID,IS,1)          &
!                       - DIFPARDX(IG)*CGO(IS,1)*ESIN(ID)    &
!	                + DIFPARDY(IG)*CGO(IS,1)*ECOS(ID)          
          PRINT*,'NOT FINISH YET. SEE SPROSD 001'
          STOP
        END IF
      ELSE
        IF(IDIFFR == 0)THEN                                        
          CAS(ID,IS,1)= CAST1*(CAST2*ART1(IG)+                            &
	 	        SUMUXHDY-SUMUYHDX-DLOC1*SUMUXDY+DLOC1*SUMUYDX)-   &
		        CAST5 *                                           &
		       (COSCOS(ID)*SUMUXDY-SINCOS(ID)*(SUMUXDX-SUMUYDY)-  &
		        SINSIN(ID)*SUMUYDX)
          CAS(ID,IS,1)=CAS(ID,IS,1)/ART1(IG)		     

          CAD(ID,IS,1) = CADT1*(ESIN(ID)*SUMHDY+ECOS(ID)*SUMHDX)      +   &
		         SINCOS(ID)*(SUMUXDY+SUMUYDX) +                   &
		         SINSIN(ID)*SUMUYDY+COSCOS(ID)*SUMUXDX
          CAD(ID,IS,1) = CAD(ID,IS,1)/ART1(IG)		   
        ELSE IF(IDIFFR == 1)THEN                                    
!          CAS(ID,IS,1)= CAST1 * (CAST2 + CAST3(ISWEEP) + CAST4(ISWEEP)) - &
!		         DIFPARAM(IG)*CAST5 *                              &
!		        (COSCOS(ID) * CAST6(ISWEEP) +                      &
!		         SINCOS(ID) * (CAST7(ISWEEP) + CAST8(ISWEEP)) +    &
!		         SINSIN(ID) * CAST9(ISWEEP) )                  

!          CAD(ID,IS,1) = DIFPARAM(IG)*CADT1 * (ESIN(ID) * CADT2(ISWEEP) - &
!			  ECOS(ID) * CADT3(ISWEEP)) -                      &
!	    	          DIFPARDX(IG)*CGO(IS,1)*ESIN(ID) +                &
!			  DIFPARDY(IG)*CGO(IS,1)*ECOS(ID) +                &
!		          SINCOS(ID) * (CADT4(ISWEEP) - CADT5(ISWEEP)) +   &
!		          SINSIN(ID) *  CADT6(ISWEEP) -                    &
!		          COSCOS(ID) *  CADT7(ISWEEP)                         
          PRINT*,'NOT FINISH YET. SEE SPROSD 002'
	  STOP
        END IF                                                        
      ENDIF
    END DO    !IDDUM
  END DO      !IS	  

!
! *** for most cases CAS and CAD will be activated. Therefore ***
! *** for IREFR is set 0 (no refraction) or ITFRE = 0 (no     ***
! *** frequency shift) we have put the IF statement outside   ***
! *** the internal loop above                                 ***
!
  IF(ITFRE == 0)THEN
    CAS(:,:,1) = 0.0
  ENDIF
!
  IF(IREFR == 0)THEN
    CAD(:,:,1) = 0.0
  ENDIF
!
  RETURN
  END SUBROUTINE SPROSD
!****************************************************************
!
  SUBROUTINE DSPHER (CAD, CAX, CAY, IG, ECOS, ESIN)       
!
!****************************************************************
!
!     computes the propagation velocities of energy in Theta-
!     space, i.e., CAD, due to use of spherical coordinates
!
!  Method
!
!     References:
!     W. E. Rogers, J. M. Kaihatu, H. A. H. Petit, N. Booij and L. H. Holthuijsen,
!     "Multiple-scale Propagation in a Third-Generation Wind Wave Model"
!     in preparation
!
!             Cg Cos(theta) Tan(latitude)
!     CAD = - ---------------------------
!                    Rearth2
!
!     The group velocity CG in the direction of the wave propagation
!     in case with a current is equal to:
!
!                     1      K(IS,IC)DEP(IX,IY)        S
!     CG(ID,IS,IC)= ( - + -----------------------) --------- +
!                     2  sinh 2K(IS,IC)DEP(IX,IY)  |k(IS,IC)|
!
!                     + (UX2(IX,IY)cos(D) + UY2(IX,IY)sin(D))
!
!     which is equivalent with CAX*cos(D) + CAY*sin(D)
!
!****************************************************************
!
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE OCPCOMM4  
  USE ALL_VARS, ONLY : VY                                                      
!
  IMPLICIT NONE
!
  REAL  :: CAD(MDC,MSC,MICMAX)                                        
  REAL  :: CAX(MDC,MSC,MICMAX)                                        
  REAL  :: CAY(MDC,MSC,MICMAX)                                        
  REAL  :: ECOS(MDC)
  REAL  :: ESIN(MDC)
  INTEGER :: IS, ID, IG
  REAL     TANLAT, CTTMP
!
!************************************************************************
!
!
! *** TANLAT is Tan of Latitude
!
  TANLAT = TAN(DEGRAD*(VY(IG)+YOFFS))
!
  DO ID = 1, MDC
    CTTMP = ECOS(ID) * TANLAT / REARTH2
    DO IS = 1, MSC
      CAD(ID,IS,1) = CAD(ID,IS,1) -                &
	        (CAX(ID,IS,1)*ECOS(ID) + CAY(ID,IS,1)*ESIN(ID)) * CTTMP  
    END DO
  END DO

  RETURN
  END SUBROUTINE DSPHER
 
!
!*******************************************************************
!
      SUBROUTINE ADDDIS (DISSXY     ,LEAKXY     ,                  &
                         AC2        ,ANYBIN     ,                  &
			 DISC0      ,DISC1      ,                  &
			 LEAKC1     ,SPCSIG     )               
!    (This subroutine has not been tested yet)
!
!*******************************************************************
!
      USE SWCOMM3 
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
!     30.72: IJsbrand Haagsma
!
!  1. Updates
!
!     20.53, Aug. 95: New subroutine
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!
!  2. Purpose
!
!     Adds dissipation and leak
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                              
!
!     IX          Counter of gridpoints in x-direction
!     IY          Counter of gridpoints in y-direction
!     MXC         Maximum counter of gridppoints in x-direction
!     MYC         Maximum counter of gridppoints in y-direction
!     MSC         Maximum counter of relative frequency
!     MDC         Maximum counter of directional distribution
!
!     REALS:
!     ---------
!
!     SP          Dummy variable
!     TEMP        Dummy variable
!
!     one and more dimensional arrays:
!     ---------------------------------
!     AC2       4D    Action density as function of D,S,X,Y and T
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
!     SWOMPU
!
! 11. Remarks
!
!     DISSXY and LEAKXY are dissipation and leak integrated over the
!     spectrum for each point in the computational grid
!     DISSC0 and DISSC1 give the dissipation distributed over the
!     spectral space in one point of the computational grid
!
! 12. Structure
!
!     -------------------------------------------------------------
!     -------------------------------------------------------------
!
! 13. Source text
!
      REAL     DISSXY(MT)    ,LEAKXY(MT)      ,          &
               DISC0(MDC,MSC)   ,DISC1(MDC,MSC)     ,          &
	       LEAKC1(MDC,MSC)  ,AC2(MDC,MSC,0:MT)            
!
      LOGICAL  ANYBIN(MDC,MSC)
      INTEGER, SAVE :: IENT=0
      CALL STRACE (IENT, 'ADDDIS')
!
      DO 100 ISC = 1, MSC
        DSDD = DDIR * FRINTF * SPCSIG(ISC)**2
        DO 90 IDC = 1, MDC
          IF (ANYBIN(IDC,ISC)) THEN
            DISSXY(KCGRD(1)) = DISSXY(KCGRD(1)) + DSDD*(DISC0(IDC,ISC) +  &
	                       DISC1(IDC,ISC) * AC2(IDC,ISC,KCGRD(1)))
            LEAKXY(KCGRD(1)) = LEAKXY(KCGRD(1)) + DSDD *                  &
	                       LEAKC1(IDC,ISC) * AC2(IDC,ISC,KCGRD(1))
          ENDIF
  90    CONTINUE
 100  CONTINUE
      RETURN
      END SUBROUTINE ADDDIS
!
!****************************************************************
!
   SUBROUTINE SPREDT (AC2     ,CAX     ,CAY     ,IDCMIN  ,IDCMAX  ,    &
	              ISSTOP  ,RDX     ,RDY     )     
!  (This subroutine has not been tested yet)
!
!****************************************************************
!
!        to estimate the action density depending of the sweep
!        direction during the first iteration of a stationary
!        computation. The reason for this is that AC2 is zero
!        at first iteration and no initialisation is given in
!        case of stationarity (NSTATC=0). Action density should
!        be nonzero because of the computation of the source
!        terms. The estimate is based on solving the equation
!
!            dN       dN
!        CAX -- + CAY -- = 0
!            dx       dy
!
!        in an explicit manner.
!
!     METHOD
!
!
!          [RDX1*CAX + RDY1*CAY]*N(i-1,j) + [RDX2*CAX + RDY2*CAY]*N(i,j-1)
! N(i,j) = ---------------------------------------------------------------
!                      (RDX1+RDX2) * CAX  +  (RDY1+RDY2) * CAY
!
!   ------------------------------------------------------------
!   For every sweep direction do,
!     For every point in S and D direction in sweep direction do,
!       predict values for action density at new point from values
!       of neighbour gridpoints taking into account spectral propagation
!       direction (with currents !!) and the boundary conditions.
!       --------------------------------------------------------
!       If wave action AC2 is negative, then
!         Give wave action initial value 1.E-10
!     ---------------------------------------------------------
!****************************************************************
!
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4 
   USE ALL_VARS, ONLY : MT                                                       

   IMPLICIT NONE
   
   INTEGER :: IS ,ID ,IDDUM ,ISSTOP                                       
   REAL    :: FAC_A ,FAC_B
   REAL    :: AC2(MDC,MSC,0:MT)                                         
   REAL    :: CAX(MDC,MSC,MICMAX)                                        
   REAL    :: CAY(MDC,MSC,MICMAX)                                        
   REAL    :: RDX(10),  RDY(10)                                          
   INTEGER :: IDCMIN(MSC) ,IDCMAX(MSC)
   REAL    :: CDEN ,CNUM ,WEIG1, WEIG2
!
   DO IS = 1, ISSTOP
     DO IDDUM = IDCMIN(IS), IDCMAX(IS)
       ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!
!      *** Computation of weighting coefs WEIG1 AND WEIG2 ***
!
       CDEN = RDX(1) * CAX(ID,IS,1) + RDY(1) * CAY(ID,IS,1)
       CNUM = (RDX(1) + RDX(2)) * CAX(ID,IS,1)              &
	    + (RDY(1) + RDY(2)) * CAY(ID,IS,1)
       WEIG1 = CDEN/CNUM
       WEIG2 = 1. - WEIG1
!
       FAC_A = WEIG1 * AC2(ID,IS,KCGRD(2))                         
       FAC_B = WEIG2 * AC2(ID,IS,KCGRD(3))                         
!
       IF (ACUPDA) AC2(ID,IS,KCGRD(1)) = MAX ( 0. , (FAC_A + FAC_B))      
!
     END DO   !IDDUM	    
   END DO     !IS  
!
   RETURN
   END SUBROUTINE SPREDT
!
!******************************************************************
!
  SUBROUTINE SWPSEL(IDCMIN    ,IDCMAX    ,CAX       ,             &
                    CAY       ,ISCMIN    ,ISCMAX    ,             &
		    IDTOT     ,ISTOT     ,IDDLOW    ,             &
		    IDDTOP    ,ISSTOP    ,DEP2      ,             &
		    UX2       ,UY2       ,SPCDIR    )
!
!******************************************************************
!
!     compute the frequency dependent counters in directional space
!     in a situation with a current and without a current.
!     The counters are only computed for the gridpoint
!     considered. This means IC = 1 
!
!******************************************************************
  USE SWCOMM1                                                         
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE OCPCOMM4                                                        
  USE ALL_VARS, ONLY : MT                                                       

  IMPLICIT NONE
!
  REAL    :: SPCDIR(MDC,6)                                               
  INTEGER :: IS ,ID ,IDSUM ,IDCLOW ,IDCHGH ,IDTOT ,ISTOT ,    & 
             IDDLOW ,IDDTOP ,ISSLOW ,ISSTOP ,IENT ,IDDUM ,ISCLOW ,    &
	     ISCHGH ,IX ,IY ,IC         
  REAL    :: CAXMID ,CAYMID ,GROUP ,UABS ,THDIR    
  INTEGER :: IDCMIN(MSC) ,IDCMAX(MSC) ,ISCMIN(MDC) ,ISCMAX(MDC) ,SECTOR(MSC)
  REAL    :: CAX(MDC,MSC,MICMAX) ,CAY(MDC,MSC,MICMAX) ,DEP2(MT) ,     &
             UX2(MT) ,UY2(MT) ,RDX(10) ,RDY(10)    
  LOGICAL :: LOWEST ,LOWBIN ,HGHBIN
!
! *** initialize arrays in frequency direction ***
!
  DO ID = 1, MDC
    ISCMIN(ID) = 1
    ISCMAX(ID) = 1
  END DO  
!
! *** initialize array's in theta direction ***
!
  DO IS = 1, MSC
    IDCMIN(IS) = 1
    IDCMAX(IS) = MDC
    SECTOR(IS) = 1
  END DO  
!
!     *** calculate minimum and maximum counters in frequency ***
!     *** space if a current is present: ISCMIN and ISCMAX    ***
!
  IDDLOW =  9999
  IDDTOP = -9999
  DO IS = 1 , MSC
    IF(SECTOR(IS) > 0)THEN
      IDDLOW = MIN ( IDDLOW , IDCMIN(IS) )
      IDDTOP = MAX ( IDDTOP , IDCMAX(IS) )
    END IF
  END DO
!
!     *** Determine counters ***
!
  DO IDDUM = IDDLOW, IDDTOP
    ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
    LOWEST = .TRUE.
    DO IS = 1, MSC
      IF(LOWEST)THEN
        ISCLOW = IS
        LOWEST = .FALSE.
      END IF
      ISCHGH = IS
    END DO
!
!       *** set the minimum and maximum counters in arrays ***
!
    IF(.NOT.LOWEST)THEN
      ISCMIN(ID) = ISCLOW
      ISCMAX(ID) = ISCHGH
    ELSE
    END IF
!
  END DO   !IDDUM
!
!     *** calculate the maximum number of counters in both ***
!     *** directional space and frequency space            ***
!
  IF(IDDLOW /= 9999)THEN
    IDTOT = ( IDDTOP - IDDLOW ) + 1
    IF(ICUR == 1)THEN
      IF(IDTOT < 3)THEN
        IDDTOP = IDDTOP + 1
        IF(IDTOT == 1) IDDLOW = IDDLOW - 1
        IDTOT = 3
      END IF
    END IF
  ELSE
    IDTOT = 0
  END IF
!
! *** set variables ***
!
  IDTOT  =     1
  ISTOT  =     1
  ISSLOW =  9999                                                     
  ISSTOP = -9999                                                     

  DO IS = 1, MSC
    IDCLOW  = 0
    IDCHGH  = 0
    IDSUM   = 0
    DO ID = 1, MDC
      IDSUM = IDSUM + 1
      ISSLOW = MIN(IS,ISSLOW)                                   
      ISSTOP = MAX(IS,ISSTOP)                                   
    END DO
  END DO  
!
  IF(ISSLOW /= 9999)THEN
    ISSLOW = 1
!   minimal value of ISSTOP is 4 (or MSC if MSC<4)                    
    IF(ICUR > 0) ISSTOP = MAX(MIN(4,MSC),ISSTOP)                    
    ISTOT = ( ISSTOP - ISSLOW ) + 1
  ELSE
    ISTOT = 0
  END IF
!
!     *** check if IDTOT is less then MDC ***
!
  IF(IDTOT > MDC)THEN
    IDDLOW = 1
    IDDTOP = MDC
    IDTOT  = MDC
  END IF
!
!     *** check if the lowest frequency is not blocked !    ***
!     *** this can occur in real cases if the depth is very ***
!     *** small and the current velocity is large           ***
!     *** the propagation velocity Cg = sqrt (gd) < U       ***
!
  IF(ICUR == 1 .AND. FULCIR .AND.                      &
     ISSLOW /= 1 .AND. ISSLOW /= 9999)THEN                          
!        CALL MSGERR (2,'The lowest freqency is blocked')                  
7002 FORMAT (I4, 1X, I4, 1X, I2)                                       
    IC = 1
    GROUP = SQRT ( GRAV_W * DEP2(KCGRD(IC)) )
  END IF

  RETURN
  END SUBROUTINE SWPSEL
