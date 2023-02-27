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

!------------------------------------------------------------------------------
# if defined (HEATING_CALCULATED)
  subroutine coare40vn(u,zu,t,zt,rh,zq,PA,ts,Rl,Rs,tau,hsb,hlb,lat,zi,rain,cp,sigH,fmiss,USR,DTER,Cdn_10) ! Siqi Li, 2021-01-27
!subroutine coare40vn(u,zu,t,zt,rh,zq,P,ts,Rl,Rs,tau,hsb,hlb,lat,zi,rain,cp,sigH,fmiss,USR,DTER)
! No-vectorized verion - Revised from the vectorized of coare40vn.m
! Zhongxiang Wu
! 4/8/2016
! Calibrated by Siqi Li on July, 2017
        
! Vectorized version of COARE 3 code (Fairall et al, 2003) with 
! modification based on the CLIMODE, MBL and CBLAST experiments 
! (Edson et al., JPO, 43, 2013). 
!
! The current version of the code includes the wind-speed, wave-age and
! sea-state dependent parameterizations of the Charnock variabile as
! described in Edson et al. (2013).  The parameterization is chosen by the 
! inputed values of
!    cp = phase speed of dominant waves (m/s)  
!  sigH = significant wave height (m)
!
! An important component of this code is whether the inputed ts 
! represents the skin temperature of a near surface temperature.  
! How this variable is treated is determined by the jcool parameter:
! set jcool=1 if Ts is bulk ocean temperature (default),
!     jcool=0 if Ts is true ocean skin temperature. 
!********************************************************************
!
! The code assumes u,t,rh,ts are vectors; 
! sensor heights zu,zt,zl, latitude lat, and PBL height zi are constants;
! air pressure P and radiation Rs,Rl may be vectors or constants. 
! Default values are assigned for P,Rs,Rl,lat,and zi if these data are not 
! available.  Input NaNs to indicate no data. Defaults should be set to 
! representative regional values if possible.
!
! Inputs:  
!
!     u = relative wind speed (m/s) at heigth zu 
!    zu = height of wind speed measurement (m)
!     t = bulk air temperature (degC) at height zt
!    zt = height of temperature measurement (m)
!    rh = relative humidity (%) at height zq
!    zq = height of relative humidity measurement (m)
!---> Siqi Li, 2021-01-27
!     P = surface air pressure (mb) (default = 1015)   (removed)
!     P = surface air pressure (Pa) (default = 101500)
!<--- Siqi Li, 2021-01-27
!    ts = water temperature (degC) see jcool below
!    Rs = downward shortwave radiation (W/m^2) (default = 150) 
!    Rl = downward longwave radiation (W/m^2) (default = 370)
!   lat = latitude (default = +45 N)
!    zi = PBL height (m) (default = 600m)
!  rain = rain rate (mm/hr) - not required for turbulent flux estimates,
!         set to NaN if unavailable.
!    cp = phase speed of dominant waves (m/s)  
!  sigH = significant wave height (m)
!

! The user controls the output.  This is currently set as:
! 
! A=[usr tau hsb hlb hbb hsbb wbar  tsr qsr zot zoq Cd Ch Ce  L 
!      1   2   3   4   5    6    7    8   9  10  11 12 13 14 15
!    zet dter dqer tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N
!     16   17   18  19  20  21  22   23   24  25 26   27 28  29   30
! Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10];
!     31     32     33 34 35   36  37  38   39        
!  where
!
!   usr = friction velocity that includes gustiness (m/s)
!   tau = wind stress (N/m^2)
!   hsb = sensible heat flux into ocean (W/m^2)
!   hlb = latent heat flux into ocean (W/m^2)
!   hbb = buoyany flux into ocean (W/m^2)
!   hsbb = "sonic" buoyancy flux measured directly by sonic anemometer 
!   wbar = the vertical velocity required to "Webb Correct" direct
!          covariance fluxes
!   tsr = temperature scaling parameter (K)
!   qsr = specific humidity scaling parameter (g/Kg)
!   zot = thermal roughness length (m)
!   zoq = moisture roughness length (m)
!   Cd = wind stress transfer (drag) coefficient at height zu   
!   Ch = sensible heat transfer coefficient (Stanton number) at height zu   
!   Ce = latent heat transfer coefficient (Dalton number) at height zu
!    L = Obukhov length scale (m) 
!  zet = Monin-Obukhov stability parameter zu/L 
! dter = cool-skin temperature depression (degC)
! dqer = cool-skin humidity depression (degC)
!  tkt = cool-skin thickness (m)
!  Urf = wind speed at reference height (user can select height below)
!  Trf = temperature at reference height
!  Qrf = specific humidity at reference height
! RHrf = relative humidity at reference height
! UrfN = neutral value of wind speed at reference height
!  Rnl = Upwelling IR radiation computed by COARE
!   Le = latent heat of vaporization
! rhoa = density of air
!   UN = neutral value of wind speed at zu
!  U10 = wind speed adjusted to 10 m
! UN10 = neutral value of wind speed at 10m
!Cdn_10 = neutral value of drag coefficient at 10m    
!Chn_10 = neutral value of Stanton number at 10m    
!Cen_10 = neutral value of Dalton number at 10m
!    RF = sensible heat flux due to rain (W/m^2)
!    Qs = surface value of specific humidity (g/kg) without cool skin
!         correction
!  Evap = evaporation rate from surface in mm/hr
!   T10 = temperarure at 10 m 
!   Q10 = specific humidity at 10 m
!  RH10 = relative humidity at 10 m
!

! Notes: 1) u is the relative wind speed, i.e., the magnitude of the
!           difference between the wind (at zu) and ocean surface current 
!           vectors.
!        2) Set jcool=0 in code if ts is true surface skin temperature,
!           otherwise ts is assumed the bulk temperature and jcool=1.
!        3) Set P=NaN to assign default value if no air pressure data 
!           available. 
!        4) Set Rs=NaN, Rl=NaN if no radiation data available.  This assigns 
!           default values to Rs, Rl so that cool skin option can be applied. 
!        5) Set lat=NaN and/or zi=NaN to assign default values if latitude
!           and/or PBL height not given. 
!        6) The code to compute the heat flux caused by precipitation is 
!           included if rain data is available (default is no rain).
!        7) Code updates the cool-skin temperature depression dter and thickness
!           tkt during iteration loop for consistency.
!        8) Number of iterations set to nits = 6.

! Reference:
!
!  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
!  Bulk parameterization of air sea fluxes: updates and verification for the 
!  COARE algorithm, J. Climate, 16, 571-590.

!  ----------------------------------------------------------
    use mod_prec
    USE MOD_HEATFLUX,ONLY : HEATING_FRESHWATER
    implicit none
	
    real(SP) :: u,zu,t,zt,rh,zq,P,ts,Rs,Rl,tau,hsb,hlb,lat,zi,rain,cp,sigH,fmiss
    real(SP) :: PA    ! Siqi Li, 2021-01-27
    real(SP) :: Pdef,Rsdef,Rldef,Latdef,Zidef
    real(SP) :: us,Qs,Q,Pv
    real(SP) :: zref,Beta,von,fdg,tdk,grav
    real(SP) :: Rgas,Le,cpa,cpv,rhoa,rhodry,visa,lapse
    real(SP) :: Al,be,cpw,rhow,visw,tcw,bigc,wetc
    real(SP) :: Rns,Rnl
    real(SP) :: rovcp,du,dt,dq,ta,tv,ug,dter,dqer,ut,u10,usr,zo10
    real(SP) :: Cd10,Ch10,Ct10,zot10,Cd,Ct,CC,Ribcu,Ribu,zetu
    real(SP) :: L10,gf,tsr,qsr,tkt
    real(SP) :: charn,charnC,umax,a1,a2,A,B,charnW,Ad,Bd,zoS,charnS
    real(SP) :: zet,L,zo,rr,zoq,zot,cdhf,cqhf,cthf,tvsr,tssr,Bf
    real(SP) :: qout,dels,qcol,alq,xlamx,sst,usr50,tsr50,qsr50,L50,zet50
    real(SP) :: dter50,dqer50,tkt50,u10N,hbb,hsbb,wbar,Evap,Ch,Ce
    real(SP) :: Cdn_10,Chn_10,Cen_10,zrf_u,zrf_t,zrf_q
    real(SP) :: psi,psi10,psirf,psiT,psi10T,psirfT,psirfQ
    real(SP) :: S,S10,Urf,UN,UrfN,UN2,U10N2,UrfN2
    real(SP) :: dwat,dtmp,dqs_dt,alfac,RF,T10,Trf,TN,T10N,TrfN,TN2,T10N2,TrfN2
    real(SP) :: SSQ,Q10,Qrf,QN,Q10N,QrfN,QN2,Q10N2,QrfN2,RHrf,RH10
    real(SP) :: qsat26sea,qsat26lake,grvf,psiu_26f,psit_26f,RHcalc
    integer :: i,nits,jcool
    integer :: waveage,seastate
    data Pdef,Rsdef,Rldef,Latdef,Zidef/ &
         101300,150,  370,  45,    600/             ! Siqi Li, 2021-01-27
!         1030,150,  370,  45,    600/

        
! convert input to column vectors
!u=u(:);t=t(:);rh=rh(:);P=P(:);ts=ts(:);
!Rs=Rs(:);Rl=Rl(:);lat=lat(:);zi=zi(:);
!zu=zu(:);zt=zt(:);zq=zq(:);
!rain=rain(:);
!N=length(u);

    jcool = 1
! set local variables to default values if input is NaN
    if(abs(p-fmiss)   < 0.00001_SP) p  =pdef
    if(abs(Rs-fmiss)  < 0.00001_SP) Rs =Rsdef
    if(abs(Rl-fmiss)  < 0.00001_SP) Rl =Rldef
    if(abs(Lat-fmiss) < 0.00001_SP) Lat=Latdef
    if(abs(Zi-fmiss)  < 0.00001_SP) Zi =Zidef

    P = PA/100.0_SP ! Siqi Li, 2021-01-27

    waveage=0
    seastate=0
!---> Siqi Li, 2021-01-27
    if(abs(cp) > 1.0E-6_SP) then
      waveage=1
      if(abs(sigh) > 1.0E-6_SP) seastate=1
    endif
!    if(abs(cp) > 0.0) then
!      waveage=1
!      if(abs(sigh) > 0.0) seastate=1
!    endif
!<--- Siqi Li, 2021-01-27

!************************************************************************
! Check on which parameterization you are using.  You can remove the 
! pause once you are familiar with how the code selects the 
! appropriate parameterization based on the input variables.
!************************************************************************
!---> Siqi Li, 2021-01-27
!    if (waveage.eq.1 .and. seastate.eq.1)  &
!       print *,' Use the seastate dependent parameterization.'
!    if (waveage.eq.1 .and. seastate.eq.0)  &
!       print *, ' Use the waveage dependent parameterization.'
!<--- Siqi Li, 2021-01-27

! input variable u is assumed relative wind speed (magnitude of difference
! between wind and surface current vectors). to follow orginal Fairall code, set
! surface current speed us=0. if us data are available, construct u prior to
! using this code.
    us = 0.0_SP
! convert rh to specific humidity
    IF(.NOT. HEATING_FRESHWATER)THEN
      Qs = qsat26sea(ts,P)/1000.0_SP    ! surface water specific humidity (g/kg) ! Siqi Li, 2021-01-27: I think the unit is kg/kg
    ELSE 
      Qs = qsat26lake(ts,P)/1000.0_SP    ! surface water specific humidity (g/kg) ! Siqi Li, 2021-01-27: I think the unit is kg/kg
    END IF  

    call qsat26air(t,P,rh,Q,Pv)   ! specific humidity of air (g/kg)
    Q=Q/1000.0_SP                      ! specific humidity of air (kg/kg)

!***********  set constants **********************************************
    zref=10.0_SP
    Beta = 1.2_SP
    von  = 0.4_SP
    fdg  = 1.0_SP ! Turbulent Prandtl number
    tdk  = 273.16_SP
    grav = grvf(lat)

!***********  air constants **********************************************
    Rgas = 287.05_SP
    Le   = (2.501_SP-.00237_SP*ts)*1e6_SP
    cpa  = 1004.67_SP
    cpv  = cpa*(1.0_SP+0.84_SP*Q)
    rhoa = P*100.0_SP/(Rgas*(t+tdk)*(1.0_SP+0.61_SP*Q))
    rhodry = (P-Pv)*100.0_SP/(Rgas*(t+tdk))
    visa = 1.326e-5_SP*(1.0_SP+6.542e-3_SP*t+8.301e-6_SP*t**2-4.84e-9_SP*t**3)
    lapse=grav/cpa

!***********  cool skin constants  ***************************************
    Al   = 2.1e-5_SP*(ts+3.2_SP)**0.79_SP
    IF(.NOT. HEATING_FRESHWATER)THEN
      be   = 0.026_SP
    ELSE
!!MDR salinity expansion coefficient, BE, to zero for freshwater
! confirmed by email with CW Fairall 4-24-2013
      be   = 0.0_SP
    END IF  

    cpw  = 4000.SP

    IF(.NOT. HEATING_FRESHWATER)THEN
        rhow = 1022._SP
    ELSE
        rhow = 1000._SP !MDR, freshwater density
    END IF  

    visw = 1e-6_SP
    tcw  = 0.6_SP
    bigc = 16._SP*grav*cpw*(rhow*visw)**3/(tcw**2*rhoa**2)
    wetc = 0.622_SP*Le*Qs/(Rgas*(ts+tdk)**2)

!***********  net radiation fluxes ***************************************
    IF(.NOT. HEATING_FRESHWATER)THEN
      Rns = 0.945_SP*Rs ! albedo correction
    ELSE
      Rns = Rs !Mark Rowe, remove albedo here because albedo is included in forcings
    END IF  

    Rnl = 0.97_SP*(5.67e-8_SP*(ts-0.3_SP*jcool+tdk)**4-Rl) ! initial value

!****************  begin bulk loop ********************************************

!***********  first guess ************************************************
    rovcp=Rgas/cpa
    lapse=grav/cpa

    du = u-us
    dt = ts-t-lapse*zt
    dq = Qs-Q     
    ta = t+tdk
    tv = ta*(1.0_SP+0.61_SP*Q)
    ug = 0.5_SP
    dter  = 0.3_SP
    dqer = dter*wetc 
    ut    = sqrt(du**2+ug**2)
    u10   = ut*log(10.0_SP/1e-4_SP)/log(zu/1e-4_SP)
    usr   = 0.035_SP*u10
    zo10  = 0.011_SP*usr**2/grav + 0.11_SP*visa/usr
    Cd10  = (von/log(10.0_SP/zo10))**2
    Ch10  = 0.00115_SP
    Ct10  = Ch10/sqrt(Cd10)
    zot10 = 10._SP/exp(von/Ct10)
    Cd    = (von/log(zu/zo10))**2
    Ct    = von/log(zt/zot10)
    CC    = von*Ct/Cd
    Ribcu = -zu/(zi*0.004_SP*Beta**3)
    Ribu  = -grav*zu/ta*((dt-dter*jcool)+.61_SP*ta*dq)/ut**2 !(ta*(tu**2)) in 26z
    if(Ribu.lt.0.0_SP) then
      zetu=CC*Ribu/(1.0_SP+Ribu/Ribcu)
    else
      zetu = CC*Ribu*(1.0_SP+27.0_SP/9.0_SP*Ribu/CC)
!      zetu=CC*Ribu/(1.+27/(9*Ribu*cc))   ! in 26z          
    endif
        
!!! Late ---- what for this        
!        k50=find(zetu>50) ! stable with very thin M-O length relative to zu
        
    L10 = zu/zetu
    gf=ut/du
    usr = ut*von/(log(zu/zo10)-psiu_26f(zu/L10))
        
!---> Changed by Siqi Li
!        tsr = -(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_26(zt/L10))
    tsr = -(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_26f(zt/L10))

!        qsr = -(dq-dqer*jcool)*von*fdg/(log(zq/zot10)-psit_26(zq/L10))
    qsr = -(dq-dqer*jcool)*von*fdg/(log(zq/zot10)-psit_26f(zq/L10))
!<--- Changed by Siqi Li
    tkt = 0.001_SP
        
!**********************************************************
!  The Charnock variable for COARE 3.0
!**********************************************************
    if(ut.gt.10.0_SP) then 
      charn=0.011_SP+(ut-10.0_SP)/(18._SP-10._SP)*(0.018_SP-0.011_SP)
    elseif(ut.gt.18._SP) then
      charn=0.018_SP
    else
      charn = 0.011_SP
    endif
!**********************************************************
!  The following gives the new formulation for the
!  Charnock variable in COARE 3.5
!**********************************************************
    charnC=0.11_SP
    umax=19._SP
    a1=0.0017_SP
    a2=-0.0050_SP
    if(u10.gt.umax) then  !wind-speed dependent coefficients
      charnC=a1*umax+a2
    else
      charnC=a1*u10+a2  
    endif
        
    A=0.114_SP  !wave-age dependent coefficients
    B=0.622_SP
    charnW=A*(usr/cp)**B

    Ad=0.091_SP  !Sea-state/wave-age dependent coefficients
    Bd=2.0_SP
    zoS=sigH*Ad*(usr/cp)**Bd
    charnS=zoS*grav/usr/usr

    nits=10 ! number of iterations
            ! Note: nits=1 in 26z_v1.0
!**************  bulk loop **************************************************

    do i=1,nits
      zet=von*grav*zu/tv*(tsr +0.61_SP*ta*qsr)/(usr**2)
      if (waveage.eq.1) then
         if (seastate.eq.1) then
            charn=charnS
         else
            charn=charnW
         endif
      else
         charn=charnC
      endif
      L=zu/zet
      zo=charn*usr*usr/grav+0.11_SP*visa/usr ! surface roughness

!---> Siqi Li, 2021-01-17
! Add the upper and lower limits for zo based on Davis et al. 2008.
! Only when wave is considered.
      if (waveage == 1) then
        zo = MAX(zo, 0.125E-6_SP)
        zo = MIN(zo, 2.850E-3_SP)
      end if
!<--- Siqi Li, 2021-01-17

      rr=zo*usr/visa
!      zoq=min(1.6e-4_SP,6e-3_SP/rr**1.6_SP) !These thermal roughness lengths give
#     if !defined (DOUBLE_PRECISION)
      zoq=amin1(1.6e-4_SP,6.e-3_SP/rr**1.6_SP) !These thermal roughness lengths give
#     else	
      zoq=dmin1(1.6e-4_SP,6.e-3_SP/rr**1.6_SP) !These thermal roughness lengths give
#     endif	
      zot=zoq                   !Stanton and Dalton numbers for COARE 4.0
      cdhf=von/(log(zu/zo)-psiu_26f(zu/L))
!---> Changed by Siqi Li
!      cqhf=von*fdg/(log(zq/zoq)-psit_26(zq/L))
      cqhf=von*fdg/(log(zq/zoq)-psit_26f(zq/L))
!      cthf=von*fdg/(log(zt/zot)-psit_26(zt/L))
      cthf=von*fdg/(log(zt/zot)-psit_26f(zt/L))
!<--- Changed by Siqi Li
      usr=ut*cdhf
      qsr=-(dq-dqer*jcool)*cqhf
      tsr=-(dt-dter*jcool)*cthf
      tvsr=tsr+0.61_SP*ta*qsr
      tssr=tsr+0.51_SP*ta*qsr
      Bf=-grav/tv*usr*tvsr
      if(Bf.gt.0.0_SP) then
!         ug=max(0.2,Beta*(Bf*zi)**0.333)
#     if !defined (DOUBLE_PRECISION)
        ug=amax1(0.2_SP,Beta*(Bf*zi)**0.333_SP)
#     else
        ug=dmax1(0.2_SP,Beta*(Bf*zi)**0.333_SP)
#     endif	 
      else
         ug=0.2_SP
      endif
      ut=sqrt(du**2+ug**2)
      gf=ut/du
      hsb=-rhoa*cpa*usr*tsr
      hlb=-rhoa*Le*usr*qsr
      qout=Rnl+hsb+hlb
      dels=Rns*(0.065_SP+11._SP*tkt-6.6e-5_SP/tkt*(1-exp(-tkt/8.0e-4_SP)))
      qcol=qout-dels
      alq=Al*qcol+be*hlb*cpw/Le
      xlamx=6.0_SP
!      tkt=min(0.01, xlamx*visw/(sqrt(rhoa/rhow)*usr))
#     if !defined (DOUBLE_PRECISION)
      tkt=amin1(0.01_SP, xlamx*visw/(sqrt(rhoa/rhow)*usr))
#     else      
      tkt=dmin1(0.01_SP, xlamx*visw/(sqrt(rhoa/rhow)*usr))
#     endif
      
      if(alq.gt.0.0_SP) then
         xlamx=6._SP/(1._SP+(bigc*alq/usr**4)**0.75_SP)**0.333_SP
         tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)
      endif
      dter=qcol*tkt/tcw
      sst=ts-dter  

      IF(.NOT. HEATING_FRESHWATER)THEN
         dqer=Qs-qsat26sea(sst,P)/1000._SP !wetc.*dter
      ELSE
         dqer=Qs-qsat26lake(sst,P)/1000._SP !wetc.*dter
      END IF	 

      Rnl=0.97_SP*(5.67e-8_SP*(ts-dter*jcool+tdk)**4-Rl) ! update dter
      if (i.eq.1) then ! save first iteration solution for case of zetu>50
         if(zetu.gt.50._SP) then ! stable with very thin M-O length relative to zu
            usr50=usr
            tsr50=tsr
            qsr50=qsr
            L50=L
            zet50=zet
            dter50=dter
            dqer50=dqer
            tkt50=tkt
         endif
      end if
      u10N = usr/von/gf*log(10._SP/zo)
      if (waveage.eq.1) then
         if (seastate.eq.1) then
            zoS=sigH*Ad*(usr/cp)**Bd-0.11_SP*visa/usr
            charnS=zoS*grav/usr/usr
         else
            charnW=A*(usr/cp)**B
         end if
      else
         charnC=a1*u10N+a2
         if(u10N.gt.umax)then
            charnC=a1*umax+a2
         endif
      end if
    end do

    ! insert first iteration solution for case with zetu>50
    if(zetu.gt.50._SP) then
       usr=usr50
       tsr=tsr50
       qsr=qsr50
       L=L50
       zet=zet50
       dter=dter50
       dqer=dqer50
       tkt=tkt50
    endif
!****************  compute fluxes  ********************************************
    tau=rhoa*usr*usr/gf      ! wind stress
    hsb=-rhoa*cpa*usr*tsr     ! sensible heat flux
    hlb=-rhoa*Le*usr*qsr      ! latent heat flux
    hbb=-rhoa*cpa*usr*tvsr    ! buoyancy flux
    hsbb=-rhoa*cpa*usr*tssr   ! sonic heat flux
    wbar=1.61_SP*hlb/Le/(1._SP+1.61_SP*Q)/rhoa+hsb/rhoa/cpa/ta !Useful for Webb Correction
    !hlwebb=rhoa*wbar*Q*Le
    Evap=1000._SP*hlb/Le/1000._SP*3600._SP   !mm/hour

!*****  compute transfer coeffs relative to ut @ meas. ht  ********************
    Cd=tau/rhoa/ut/max(.1_SP,du)
    Ch=-usr*tsr/ut/(dt-dter*jcool)
    Ce=-usr*qsr/(dq-dqer*jcool)/ut

!***  compute 10-m neutral coeff relative to ut (output if needed) ************
!---> Siqi Li, 2021-01-27 : there is no need to time 1000.
    Cdn_10=von**2./log(10._SP/zo)**2  
    Chn_10=von**2*fdg/log(10._SP/zo)/log(10._SP/zot)
    Cen_10=von**2*fdg/log(10._SP/zo)/log(10._SP/zoq)
!    Cdn_10=1000*von**2./log(10./zo)**2 
!    Chn_10=1000*von**2*fdg/log(10./zo)/log(10./zot)
!    Cen_10=1000*von**2*fdg/log(10./zo)/log(10./zoq)
!<--- Siqi Li, 2021-01-27

!***  compute 10-m neutral coeff relative to ut (output if needed) ************
!  Find the stability functions
!********************************
    zrf_u=2._SP             !User defined reference heights to
    zrf_t=2._SP             !compute values at zrf.
    zrf_q=2._SP
    psi=psiu_26f(zu/L)
    psi10=psiu_26f(10._SP/L)
    psirf=psiu_26f(zrf_u/L)
    psiT=psit_26f(zt/L)
    psi10T=psit_26f(10._SP/L)
    psirfT=psit_26f(zrf_t/L)
    psirfQ=psit_26f(zrf_q/L)
    gf=ut/du

!*********************************************************
!  Determine the wind speeds relative to ocean surface
!  Note that usr is the friction velocity that includes 
!  gustiness usr = sqrt(Cd) S, which is equation (18) in
!  Fairall et al. (1996)
!*********************************************************
    S = ut
    U = du
    S10 = S + usr/von*(log(10._SP/zu)-psi10+psi)
    U10 = S10/gf
    ! or U10 = U + usr/von/gf*(log(10/zu)-psi10+psi)
    Urf = U + usr/von/gf*(log(zrf_u/zu)-psirf+psi)
    UN = U + psi*usr/von/gf
    U10N = U10 + psi10*usr/von/gf
    UrfN = Urf + psirf*usr/von/gf
        
    UN2 = usr/von/gf*log(zu/zo)
    U10N2 = usr/von/gf*log(10._SP/zo)
    UrfN2  = usr/von/gf*log(zrf_u/zo)
        
!******** rain heat flux (save to use if desired) *****************************
    if(rain.gt.0._SP) then
       dwat=2.11e-5_SP*((t+tdk)/tdk)**1.94_SP !! water vapour diffusivity
       dtmp=(1._SP + 3.309e-3_SP*t - 1.44e-6_SP*t*t)*0.02411_SP/(rhoa*cpa) !! heat diffusivity
       dqs_dt=Q*Le/(Rgas*(t+tdk)**2) !! Clausius-Clapeyron
       alfac= 1._SP/(1._SP+0.622_SP*(dqs_dt*Le*dwat)/(cpa*dtmp)) !! wet bulb factor
       RF= rain*alfac*cpw*((ts-t-dter*jcool)+ &
            (Qs-Q-dqer*jcool)*Le/cpa)/3600._SP
    else
       RF=0._SP
    end if

    lapse=grav/cpa
    SST=ts-dter*jcool

    T = t
    T10 = T + tsr/von*(log(10._SP/zt)-psi10T+psiT) + lapse*(zt-10._SP)
    Trf = T + tsr/von*(log(zrf_t/zt)-psirfT+psiT) + lapse*(zt-zrf_t)
    TN = T + psiT*tsr/von
    T10N = T10 + psi10T*tsr/von
    TrfN = Trf + psirfT*tsr/von

    TN2 = SST + tsr/von*log(zt/zot)-lapse*zt
    T10N2 = SST + tsr/von*log(10._SP/zot)-lapse*10._SP
    TrfN2 = SST + tsr/von*log(zrf_t/zot)-lapse*zrf_t

    IF(.NOT. HEATING_FRESHWATER)THEN
      dqer=(Qs-qsat26sea(SST,P)/1000._SP)*jcool !wetc*dter*jcool
    ELSE
      dqer=(Qs-qsat26lake(SST,P)/1000._SP)*jcool !wetc*dter*jcool
    END IF

    SSQ=Qs-dqer
    SSQ=SSQ*1000._SP
    Q=Q*1000._SP
    qsr=qsr*1000._SP
    Q10 = Q + qsr/von*(log(10._SP/zq)-psi10T+psiT)
    Qrf = Q + qsr/von*(log(zrf_q/zq)-psirfQ+psiT)
    QN = Q + psiT*qsr/von/sqrt(gf)
    Q10N = Q10 + psi10T*qsr/von
    QrfN = Qrf + psirfQ*qsr/von
        
    QN2 = SSQ + qsr/von*log(zq/zoq)
    Q10N2 = SSQ + qsr/von*log(10._SP/zoq)
    QrfN2 = SSQ + qsr/von*log(zrf_q/zoq)
    RHrf=RHcalc(Trf,P,Qrf/1000._SP)
    RH10=RHcalc(T10,P,Q10/1000._SP)

!--->Siqi Li
    hsb=-hsb
    hlb=-hlb
!<---Siqi Li

!****************  output  ****************************************************

!A=[usr tau hsb hlb hbb hsbb wbar  tsr qsr zot zoq Cd Ch Ce  L zet dter dqer tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 gf]
!   1   2   3   4   5   6    7      8   9  10  11  12 13 14 15 16   17   18  19  20  21  22  23   24   25 26  27  28  29  30    31     32     33   34 35  36  37  38   39  40]
      end subroutine coare40vn
!------------------------------------------------------------------------------
!      function psi=psit_26(zet)
      function psit_26f(zet)        
    use mod_prec
! computes temperature structure function
        implicit none
	real(SP) :: zet
	real(SP) :: dzet,x,psik,psic,f
	real(SP) :: psit_26f
	
!        dzet=min(50.,0.35*zet) ! stable
#       if !defined (DOUBLE_PRECISION)
        dzet=amin1(50._SP,0.35_SP*zet) ! stable
#       else
        dzet=dmin1(50._SP,0.35_SP*zet) ! stable
#       endif	

        psit_26f=-((1._SP+0.6667_SP*zet)**1.5+0.6667_SP*(zet-14.28_SP)*exp(-dzet)+8.525_SP)
        
!        k=find(zet<0) ! unstable
        if(zet.lt.0.0_SP) then ! unstable
           x=(1._SP-16._SP*zet)**0.5
           psik=2._SP*log((1._SP+x)/2._SP)
           x=(1._SP-34.15_SP*zet)**0.3333
           psic=1.5_SP*log((1._SP+x+x**2)/3._SP)-sqrt(3._SP) &
                *atan((1._SP+2._SP*x)/sqrt(3._SP))+4._SP*atan(1._SP)/sqrt(3._SP)
           f=zet**2./(1._SP+zet**2)
           psit_26f=(1._SP-f)*psik+f*psic
        endif
      end function psit_26f
!------------------------------------------------------------------------------
!      function psi=psiu_26(zet)
      function psiu_26f(zet)
    use mod_prec
! computes velocity structure function
        implicit none
	real(SP) :: zet
	real(SP) :: a,b,c,d,dzet,x,psik,psic,f
	real(SP) :: psiu_26f
	
!        dzet=min(50.,0.35*zet) ! stable
#       if !defined (DOUBLE_PRECISION)
        dzet=amin1(50._SP,0.35_SP*zet) ! stable
#       else	
        dzet=dmin1(50._SP,0.35_SP*zet) ! stable
#       endif
	
        a=0.7_SP
        b=3._SP/4._SP
        c=5._SP
        d=0.35_SP
        psiu_26f=-(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)
!        k=find(zet<0) ! unstable
        if(zet.lt.0._SP) then ! unstable
           x=(1._SP-16._SP*zet)**0.25
           psik=2.0_SP*log((1.0_SP+x)/2.0_SP)+log((1.0_SP+x*x)/2.0_SP)-2.0_SP*atan(x)+2.0_SP*atan(1.0_SP)
           x=(1._SP-10.15_SP*zet)**0.3333
           psic=1.5_SP*log((1._SP+x+x**2)/3._SP)-sqrt(3._SP) &
                *atan((1._SP+2._SP*x)/sqrt(3._SP))+4.0_SP*atan(1._SP)/sqrt(3._SP)
           f=zet**2./(1._SP+zet**2)
           psiu_26f=(1._SP-f)*psik+f*psic
        endif
      end function psiu_26f
!------------------------------------------------------------------------------
!      function psi=psiu_40(zet)
      function psiu_40(zet)				
    use mod_prec
! computes velocity structure function
        implicit none
	real(SP) :: zet
	real(SP) :: a,b,c,d,dzet,x,psik,psic,f
	real(SP) :: psiu_40
	
!        dzet=min(50.,0.35*zet) ! stable
#       if !defined (DOUBLE_PRECISION)
        dzet=amin1(50._SP,0.35_SP*zet) ! stable
#       else	
        dzet=dmin1(50._SP,0.35_SP*zet) ! stable
#       endif
        a=1._SP
        b=3._SP/4._SP
        c=5._SP
        d=0.35_SP
        psiu_40=-(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)
        !k=find(zet<0) ! unstable
        if(zet.lt.0._SP) then
           x=(1._SP-18._SP*zet)**0.25
           psik=2._SP*log((1._SP+x)/2._SP)+log((1._SP+x*x)/2._SP)-2._SP*atan(x)+2._SP*atan(1._SP)
           x=(1._SP-10._SP*zet)**0.3333
           psic=1.5_SP*log((1._SP+x+x**2)/3._SP)-sqrt(3.0_SP)  &
                *atan((1._SP+2._SP*x)/sqrt(3._SP))+4.0_SP*atan(1._SP)/sqrt(3._SP)
           f=zet**2./(1._SP+zet**2)
           psiu_40=(1._SP-f)*psik+f*psic
        endif
      end function psiu_40
!------------------------------------------------------------------------------
!      function exx=bucksat(T,P)
      function bucksat(T,P)
    use mod_prec
! computes saturation vapor pressure [mb]
! given T [degC] and P [mb]
        implicit none
	real(SP) :: T,P
	real(SP) :: bucksat
	
        bucksat=6.1121_SP*exp(17.502_SP*T/(T+240.97_SP))*(1.0007_SP+3.46e-6_SP*P)
      end function bucksat
!------------------------------------------------------------------------------
      !function qs=qsat26sea(T,P)
      function qsat26sea(T,P)
    use mod_prec
! computes surface saturation specific humidity [g/kg]
! given T [degC] and P [mb]
        implicit none
	real(SP) :: T,P
	real(SP) :: ex,es,bucksat
	real(SP) :: qsat26sea
	
        ex=bucksat(T,P)
        es=0.98_SP*ex ! reduction at sea surface
        qsat26sea=622._SP*es/(P-0.378_SP*es)
      end function qsat26sea
!------------------------------------------------------------------------------
      !function qs=qsat26lake(T,P)
      function qsat26lake(T,P)
    use mod_prec
! computes surface saturation specific humidity [g/kg]
! given T [degC] and P [mb]
        implicit none
	real(SP) :: T,P
	real(SP) :: ex,es,bucksat
	real(SP) :: qsat26lake
	
        ex=bucksat(T,P)
!JQI        es=0.98*ex ! reduction at sea surface
        es=1.0_SP*ex ! reduction at sea surface !MDR, don't apply 0.98 for freshwater
        qsat26lake=622._SP*es/(P-0.378_SP*es)
      end function qsat26lake
!------------------------------------------------------------------------------
      subroutine qsat26air(T,P,rh,q,em)
    use mod_prec  
! computes saturation specific humidity [g/kg]
! given T [degC] and P [mb]
        implicit none
	real(SP) :: T,P,rh
	real(SP) :: es,bucksat
	real(SP) :: q,em
	
        es=bucksat(T,P)
        em=0.01_SP*rh*es
        q=622._SP*em/(P-0.378_SP*em)
      end subroutine qsat26air
!------------------------------------------------------------------------------
!      function g=grv(lat)
      function grvf(lat)
    use mod_prec
! computes g [m/sec**2] given lat in deg
        implicit none
        real(SP) :: lat
	real(SP) :: c1,c2,c3,c4,gamma,phi,x
	real(SP) :: grvf
        real(SP), parameter :: pi=3.1415926_SP
        
        gamma=9.7803267715_SP
        c1=0.0052790414_SP
        c2=0.0000232718_SP
        c3=0.0000001262_SP
        c4=0.0000000007_SP
        phi=lat*pi/180._SP
        x=sin(phi)
        grvf=gamma*(1._SP+c1*x**2+c2*x**4+c3*x**6+c4*x**8)
      end function grvf

!------------------------------------------------------------------------------
!      function RHrf=RHcalc(T,P,Q)
      function RHcalc(T,P,Q)
    use mod_prec
! computes relative humidity given T,P, & Q
        implicit none
        real(SP) :: T,P,Q
        real(SP) :: es,em
        real(SP) :: RHcalc

        es=6.1121_SP*exp(17.502_SP*T/(T+240.97_SP))*(1.0007_SP+3.46e-6_SP*P)
        em=Q*P/(0.378_SP*Q+0.622_SP)
        RHcalc=100._SP*em/es
      end function RHcalc

# endif
! ---------------------------------------------- 
