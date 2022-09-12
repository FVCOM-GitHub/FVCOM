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

!==============================================================================|
# if defined (HEATING_CALCULATED)
   SUBROUTINE COARE26Z(UR,ZU,TA,ZT,RH,ZQ,PA,TS,DLW,DSW,TAU,HSB,HLB,LAT,USR,DTER,Cdn_10) ! Siqi Li, 2021-01-27
!   SUBROUTINE COARE26Z(UR,ZU,TA,ZT,RH,ZQ,PA,TS,DLW,DSW,TAU,HSB,HLB,LAT,USR,DTER) !mdr include USR, DTER
!  A=coare26(u,zu,Ta,zt,rh,zq,Pa,Ts,dlw,dsw)
! Simplified non-vectorized version of coare2.6 code
! with cool skin option retained but warm layer and 
! surface wave options removed, and rain set to zero. 
! Assumes input are single scalars and that u is the
! magitude of the difference between wind and surface 
! current vectors, if latter available. Output:
! A = [tau qsen qlat Cd Ch Ce Cdn_10 Chn_10 Cen_10].

   USE MOD_PREC
   USE MOD_HEATFLUX, ONLY : HEATING_FRESHWATER
   
   IMPLICIT NONE
    
   REAL(SP) :: UR,ZU,TA,ZT,RH,ZQ,PA,TS,DLW,DSW,TAU,HSB,HLB,LAT,LE,L,JCOOL,L10
   REAL(SP) :: U,US,T,P,RL,RS,QSAT26S,QS,Q,ZI,BETA,VON,FDG,TDK,GRVS,GRAV,RGAS
   REAL(SP) :: CPA,CPV,RHOA,VISA,AL,BE,CPW,RHOW,VISW,TCW,BIGC,WETC,RNS,RNL,DU
   REAL(SP) :: DT,DQ,UG,DTER,DQER,UT,U10,USR,ZO10,CD10,CH10,CT10,ZOT10,CD,CT,CC
   REAL(SP) :: RIBCU,RIBU,NITS,ZETU,PSIU_26S,PSIT_26S,TSR,QSR,TKT,CHARN,ZET
   REAL(SP) :: ZO,ZOQ,ZOT,BF,QOUT,DELS,QCOL,ALQ,XLAMX,CH,CE,CDN_10,CHN_10,CEN_10,RR
   
   INTEGER  :: I
! Set jcool=0 if Ts is surface, =1 if Ts is bulk.
! rcb checked 6/9/04
! set jcool=1 if Ts is bulk, 0 if Ts is true skin jcool=1;
   JCOOL=1.
! rename variables from fairall et al coare3 code


! wind speed (m/s) at height zr (m)
   U=UR
!  surface current speed in the wind direction(m/s)
   US=0.*UR
!  water temperature (deg C)
   TS=TS
!  BULK AIR TEMPERATURE (C) AT HEIGHT ZT(m)
   T=TA
!  RELATIVE HUMIDITY (%) AT HEIGHT zq(M)
   RH=RH
!  SURFACE PRESSURE (mb)
!   P=PA
   P=PA/100.  !ejw 8/15/2011 from pa to mb in this subroutine only
!  DOWNWARD LONGWAVE RADIATION (W/m2)
   RL=DLW
!  DOWNWARD SHORTWAVE RADIATION (W/m2)
   RS=DSW
!  CONVERT RH TO SPECIFIC HUMIDITY (G/KG)
   CALL QSAT26(TS,P,QSAT26S)

   IF(.NOT. HEATING_FRESHWATER)THEN
     QS=0.98_SP*QSAT26S/1000._SP
   ELSE 
     QS=QSAT26S/1000._SP !MDR, don't apply 0.98 for freshwater
   END IF  

!  SPECIFIC HUMIDITY OF AIR (G/KG)  
   CALL QSAT26(T,P,QSAT26S)
   Q=(0.01_SP*RH)*QSAT26S/1000._SP

!   SET RAIN TO ZERO RAIN=0*U
!   SET RAIN RATE (MM/HR) - KEEP AS OPTION

!  ***********SET LOCAL CONSTANTS *********
!    PBL HEIGHT (M)
   ZI=600._SP
!     LATITUDE (DEG,N=+)- GEORGES BANK
!        LAT=42.
! ************SET CONSTANTS **************
   BETA=1.2_SP
   VON=0.4_SP
   FDG=1.00_SP
   TDK=273.16_SP
   CALL GRV(LAT,GRVS)
   GRAV=GRVS
! ************AIR CONSTANTS **************
   RGAS=287.1_SP
   LE=(2.501_SP-0.00237_SP*TS)*1000000._SP
   CPA=1004.67_SP
   CPV=CPA*(1._SP+0.84_SP*Q)
   RHOA=P*100._SP/(RGAS*(T+TDK)*(1+0.61_SP*Q))
   VISA=1.326_SP*0.00001_SP*(1+6.542_SP*0.001_SP*T+  &
        8.301_SP*0.000001_SP*T*T-4.84_SP*0.000000001_SP*T*T*T)
! ***********COOL SKIN CONSTANTS******************
   AL=2.1_SP*0.00001_SP*((TS+3.2_SP)**0.79_SP)
   IF(.NOT. HEATING_FRESHWATER)THEN
     BE=0.026_SP
   ELSE
!!MDR salinity expansion coefficient, BE, to zero for freshwater
! confirmed by email with CW Fairall 4-24-2013
     BE=0.0_SP 
   END IF  
   CPW=4000._SP
   IF(.NOT. HEATING_FRESHWATER)THEN
     RHOW=1022._SP
   ELSE
     RHOW=1000._SP !MDR, freshwater density
   END IF  

   VISW=0.000001_SP
   TCW=0.6_SP
   BIGC=16._SP*GRAV*CPW*((RHOW*VISW)**3)/(TCW*TCW*RHOA*RHOA)
   WETC=0.622_SP*LE*QS/(RGAS*((TS+TDK)**2))
! ***************COMPUTE AUX STUFF *********
   IF(.NOT. HEATING_FRESHWATER)THEN
     RNS=RS*0.945_SP
   ELSE
     RNS=RS !Mark Rowe, remove albedo here because albedo is included in forcings
   END IF  

   RNL=0.97_SP*(5.67_SP*0.00000001_SP*((TS-0.3_SP*JCOOL+TDK)**4)-RL)


! **************BEGIN BULK LOOP ***********

! **************FIRST GUESS ***************
   DU=U-US
   DT=TS-T-0.0098_SP*ZT
   DQ=QS-Q
   TA=T+TDK
   UG=0.5_SP
   DTER=0.3_SP
   DQER=WETC*DTER
   UT=SQRT(DU*DU+UG*UG) 
   U10=UT*LOG(10._SP/1e-4)/LOG(ZU/1e-4)
   USR=0.035_SP*U10

   ZO10=0.011_SP*USR*USR/GRAV+0.11_SP*VISA/USR
   CD10=(VON/LOG(10._SP/ZO10))**2
   CH10=0.00115_SP
   CT10=Ch10/SQRT(CD10)
   ZOT10=10._SP/EXP(VON/CT10)
   CD=(VON/LOG(ZU/ZO10))**2
   CT=VON/LOG(ZT/ZOT10)
   CC=VON*CT/CD
   RIBCU=-ZU/(ZI*0.004_SP*(BETA**3))
   RIBU=-GRAV*ZU*((DT-DTER*JCOOL)+.61_SP*TA*DQ)/(TA*(UT**2))
!       same as edson
!MDR   NITS=6._SP   
   NITS=3._SP !MDR set NITS to 3, after F90 version of COARE 3.0
   IF(RIBU < 0)THEN
     ZETU=CC*RIBU/(1._SP+RIBU/RIBCU)
   ELSE
!MDR     ZETU=CC*RIBU/(1._SP+27._SP/(9*RIBU*CC))
     ZETU=CC*RIBU*(1.0+27.0/9.0*RIBU/CC) !MDR, corrected from F90 version
   END IF

   L10=ZU/ZETU
   IF(ZETU > 50)THEN
     NITS=1
   END IF
   CALL PSIU_26(ZU/L10,PSIU_26S)
   USR=UT*VON/(LOG(ZU/ZO10)-PSIU_26S)
   CALL PSIT_26(ZT/L10,PSIT_26S)
   TSR=-(DT-DTER*JCOOL)*VON*FDG/(LOG(ZT/ZOT10)-PSIT_26S)
   CALL PSIT_26(ZQ/L10,PSIT_26S)
   
   QSR=-(DQ-WETC*DTER*JCOOL)*VON*FDG/(LOG(ZQ/ZOT10)-PSIT_26S)
   TKT=.001_SP
   CHARN=0.011_SP
   
   IF(UT > 10)THEN 
     CHARN=0.011_SP+(UT-10)/(18-10)*(0.018_SP-0.011_SP)
   END IF
   
   IF(UT > 18)THEN
     CHARN=0.018_SP
   END IF

!*************** bulk loop ************
    do I=1,nits !MDR uncomment to do iterations, was hardwired to 1 iteration
!MDR   DO I=1,1
     ZET=VON*GRAV*ZU/TA*(TSR*(1._SP+0.61_SP*Q)+0.61_SP*TA*QSR)/   &
         (USR*USR)/(1._SP+0.61_SP*Q)
     ZO=CHARN*USR*USR/GRAV+0.11_SP*VISA/USR
     RR=ZO*USR/VISA
     L=ZU/ZET
#    if !defined(DOUBLE_PRECISION)
     ZOQ=AMIN1(1.15e-4_SP,5.5e-5_SP/(RR**0.6_SP)) 
#    else
     ZOQ=DMIN1(1.15e-4_SP,5.5e-5_SP/(RR**0.6_SP)) 
#    endif
     ZOT=ZOQ
     CALL PSIU_26(ZU/L,PSIU_26S)
     USR=UT*VON/(LOG(ZU/ZO)-PSIU_26S)
     CALL PSIT_26(ZT/L,PSIT_26S)
     TSR=-(DT-DTER*JCOOL)*VON*FDG/(LOG(ZT/ZOT)-PSIT_26S)
     CALL PSIT_26(ZQ/L,PSIT_26S)
     QSR=-(DQ-WETC*DTER*JCOOL)*VON*FDG/(LOG(ZQ/ZOQ)-PSIT_26S) 
     BF=-GRAV/TA*USR*(TSR+.61_SP*TA*QSR)               

     IF(BF > 0)THEN
       UG=BETA*((BF*ZI)**0.333_SP)
     ELSE
       UG=0.2_SP
     END IF
     
     UT=SQRT(DU*DU+UG*UG) 
     RNL=0.97_SP*(5.67_SP*0.00000001_SP*((TS-DTER*JCOOL+TDK)**4)-RL)
     HSB=-RHOA*CPA*USR*TSR
     HLB=-RHOA*LE*USR*QSR
     QOUT=RNL+HSB+HLB
 
!     DELS=RNS*(.065_SP+11*TKT-6.6_SP*0.00001_SP/(TKT*(1-EXP(-TKT/8.0_SP*0.0001_SP))))
!Mark Rowe 4-24-2013 the line above has misplaced parentheses in two places, and was giving 
!DELS = -Inf  By comparison to COARE3.0 Matlab version, below is correct
     DELS=RNS*(.065_SP+11.0_SP*TKT-6.6_SP*0.00001_SP/TKT*(1.0_SP-EXP(-TKT/(8.0_SP*0.0001_SP))))
     QCOL=QOUT-DELS
     ALQ=AL*QCOL+BE*HLB*CPW/LE

     IF(ALQ > 0)THEN
       XLAMX=6._SP/(1._SP+(BIGC*ALQ/USR**4)**0.75_SP)**0.333_SP
       TKT=XLAMX*VISW/(SQRT(RHOA/RHOW)*USR)
     ELSE
       XLAMX=6.0_SP
#    if !defined(DOUBLE_PRECISION)
       TKT=AMIN1(.01_SP,XLAMX*VISW/(SQRT(RHOA/RHOW)*USR))
#    else
       TKT=DMIN1(.01_SP,XLAMX*VISW/(SQRT(RHOA/RHOW)*USR))
#    endif
     END IF
!           cool skin

     DTER=QCOL*TKT/TCW
     DQER=WETC*DTER

   END DO
!  of  end bulk iter loop

!****** compute fluxes ******************************************
!             wind stress
   TAU=RHOA*USR*USR*DU/UT
!              sensible heat flux
   HSB=RHOA*CPA*USR*TSR
!              latent heat flux
   HLB=RHOA*LE*USR*QSR


!****** compute transfer coeffs relative to ut @ meas. ht ********
#  if !defined(DOUBLE_PRECISION)
   CD=TAU/RHOA/UT/AMAX1(.1_SP,DU)
#  else
   CD=TAU/RHOA/UT/DMAX1(.1_SP,DU)
#  endif
   CH=-USR*TSR/UT/(DT-DTER*JCOOL)
   CE=-USR*QSR/(DQ-DQER*JCOOL)/UT

!****** compute 10-m neutral coeff relative to ut ****************
   CDN_10=VON*VON/LOG(10._SP/ZO)/LOG(10._SP/ZO)
   CHN_10=VON*VON*FDG/LOG(10._SP/ZO)/LOG(10._SP/ZOT)
   CEN_10=VON*VON*FDG/LOG(10._SP/ZO)/LOG(10._SP/ZOQ)
   
!******** rain heat flux (save to use if desired) *************
! dwat=2.11e-5*((t+tdk)/tdk)^1.94; %! water vapour diffusivity
! dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa); %!heat diffusivity
! alfac= 1/(1+(wetc*Le*dwat)/(cpa*dtmp)); %! wet bulb factor
! RF= rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600;
!**************************************************************
!----------------------------------------------------------
   RETURN
   END SUBROUTINE COARE26Z
   

!-----------------------------------------------------------------------|
   SUBROUTINE PSIT_26(ZET,PSI)
! computes temperature structure function
   USE MOD_PREC
   
   IMPLICIT NONE
   REAL(SP) :: ZET,PSI,X,PSIK,PSIC,F,C
   
   IF(ZET < 0.0_SP)THEN
     X=(1._SP-15._SP*ZET)**.5_SP
     PSIK=2._SP*LOG((1._SP+X)/2._SP)
     X=(1._SP-34.15_SP*ZET)**.3333_SP
     PSIC=1.5_SP*LOG((1._SP+X+X*X)/3._SP)-SQRT(3._SP)    &
          *ATAN((1._SP+2._SP*X)/SQRT(3._SP))+4._SP*ATAN(1._SP)/SQRT(3._SP)
     F=ZET*ZET/(1._SP+ZET*ZET)
     PSI=(1._SP-F)*PSIK+F*PSIC
   ELSE
#    if !defined(DOUBLE_PRECISION)
     C=AMIN1(50._SP,.35_SP*ZET)
#    else
     C=DMIN1(50._SP,.35_SP*ZET)
#    endif
     PSI=-((1._SP+2._SP/3._SP*ZET)**1.5_SP+.6667_SP     &
         *(ZET-14.28_SP)/EXP(C)+8.525_SP)
   END IF
   
   RETURN
   END SUBROUTINE PSIT_26
  

!----------------------------------------------------------
   SUBROUTINE PSIU_26(ZET,PSI)
! computes velocity structure function
   USE MOD_PREC
   
   IMPLICIT NONE
   REAL(SP) :: ZET,PSI,X,PSIK,PSIC,F,C

   IF(ZET < 0.0_SP)THEN
     X=(1._SP-15._SP*ZET)**.25_SP
     PSIK=2._SP*LOG((1._SP+X)/2._SP)+LOG((1._SP+X*X)/2._SP)   &
          -2._SP*ATAN(X)+2._SP*ATAN(1._SP)
     X=(1._SP-10.15_SP*ZET)**.3333_SP
     PSIC=1.5_SP*LOG((1._SP+X+X*X)/3._SP)-SQRT(3._SP)   &
          *ATAN((1._SP+2._SP*X)/SQRT(3._SP))+4._SP*ATAN(1._SP)/SQRT(3._SP)
     F=ZET*ZET/(1._SP+ZET*ZET)
     PSI=(1._SP-F)*PSIK+F*PSIC
   ELSE
#    if !defined(DOUBLE_PRECISION)
     C=AMIN1(50._SP,.35_SP*ZET)
#    else
     C=DMIN1(50._SP,.35_SP*ZET)
#    endif
     PSI=-((1._SP+1.0_SP*ZET)**1.0_SP+.667_SP*(ZET-14.28_SP)/EXP(C)+8.525_SP)
   END IF
   
   RETURN
   END SUBROUTINE PSIU_26


!-----------------------------------------------------------
   SUBROUTINE QSAT26(T,P,QS)
! computes saturation specific humidity
   USE MOD_PREC
   
   IMPLICIT NONE
   REAL(SP) :: T,P,QS,ES
   
   ES=6.112_SP*EXP(17.502_SP*T/(T+241.0_SP))     &
      *(1.0007_SP+3.46_SP*0.000001_SP*P)
   QS=ES*622/(P-.378_SP*ES)
   
   RETURN
   END SUBROUTINE QSAT26


!------function g=grv(lat)----------------------------------------------|
   SUBROUTINE GRV(LAT,G)
   USE MOD_PREC
   
   IMPLICIT NONE
   REAL(SP) :: PI,GAMMA,C1,C2,C3,C4,PHI,X,G,LAT

   PI=3.1415926_SP
! computes g given lat in deg
   GAMMA=9.7803267715_SP
   C1=0.0052790414_SP
   C2=0.0000232718_SP
   C3=0.0000001262_SP
   C4=0.0000000007_SP
   PHI=LAT*PI/180._SP
   X=SIN(PHI)
   G=GAMMA*(1._SP+C1*X**2+C2*X**4+C3*X**6+C4*X**8)
   
   RETURN
   END SUBROUTINE GRV

# endif 
