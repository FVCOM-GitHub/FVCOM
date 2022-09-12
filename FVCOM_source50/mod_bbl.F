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

!=======================================================================
! FVCOM Bottom Boundary Layer Module 
!
! Copyright:    2009(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Contact:      G. Cowles 
!               School for Marine Science and Technology, Umass-Dartmouth
!
! Based on the Bottom Boundary Layer with Wave-Current Interaction as implemented
!     in ROMS by J. Warner (USGS)
!
! Comments:  Bottom Boundary Layer Module 
!
! Current FVCOM dependency
!
!   init_sed.F:  - user defined sediment model initial conditions
!   mod_ncdio.F: - netcdf output includes concentration/bottom/bed fields
!   fvcom.F:     - main calls sediment setup 
!   internal.F:  - calls sediment advance
!
! History

!
!  Later

!
!=======================================================================
Module Mod_BBL

# if defined (SEDIMENT) && (WAVE_CURRENT_INTERACTION)

Use Mod_Prec 
Use Mod_Types
Use ALL_VARS
Use Lims
# if defined (ORIG_SED)
Use Mod_Sed
# elif defined (CSTMS_SED)
Use Mod_Sed_CSTMS
# endif
USE TIMECOMM
USE SWCOMM3
USE VARS_WAVE
USE MOD_WAVE_CURRENT_INTERACTION

implicit none

public

!=======================================================================
!                                                                      !
!  Ubot         Wind-induced, bed wave orbital U-velocity (m/s) at     !
!                 RHO-points.                                          !
!  Ur           Bottom U-momentum above bed (m/s) at RHO-points.       !
!  Vbot         Wind-induced, bed wave orbital V-velocity (m/s) at     !
!                 RHO-points.                                          !
!  Vr           Bottom V-momentum above bed (m/s) at RHO-points.       !
!  bustrc       Kinematic bottom stress (m2/s2) due currents in the    !
!                 XI-direction at RHO-points.                          !
!  bustrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the XI-direction at horizontal RHO-points.     !
!  bustrcwmax   Kinematic bottom stress (m2/s2) due to maximum wind    !
!                 and currents in the XI-direction at RHO-points.      !
!  bvstrc       Kinematic bottom stress (m2/s2) due currents in the    !
!                 ETA-direction at RHO-points.                         !
!  bvstrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the ETA-direction at horizontal RHO-points.    !
!  bvstrcwmax   Kinematic bottom stress (m2/s2) due to maximum wind    !
!                 and currents in the ETA-direction RHO-points.        !
!                                                                      !
!=======================================================================

  integer , allocatable :: Iconv(:)


 character(len=120), parameter :: bbl_model_version = "BBL Module with Sediment and Wave"
 character(len=120) :: bblfile

!    vonKar        von Karman constant.
  real(sp), parameter  :: vonKar  = 0.41            ! non-dimensional
  real(sp), parameter  :: Cdb_max = 0.5
  real(sp), parameter  :: Cdb_min = 0.000001
!-----------------------------------------------------------------------
!  Closure parameters associated with Styles and Glenn (1999) bottom
!  currents and waves boundary layer.
!-----------------------------------------------------------------------
!
!    sg_Cdmax      Upper limit on bottom darg coefficient.
!    sg_alpha      Free parameter indicating the constant stress
!                    region of the wave boundary layer.
!    sg_g          Acceleration of gravity (m/s2).
!    sg_kappa      Von Karman constant.
!    sg_mp         Nondimensional closure constant.
!    sg_n          Maximum number of iterations for bisection method.
!    sg_nu         Kinematic viscosity of seawater (m2/s).
!    sg_pi         Ratio of circumference to diameter.
!    sg_tol        Convergence criterion.
!    sg_ustarcdef  Default bottom stress (m/s).
!    sg_z100       Depth (m), 100 cm above bottom.
!    sg_z1p        Nondimensional closure constant.
!    sg_zrmin      Minimum allowed height (m) of current above bed.
!                    Otherwise, logarithmic interpolation is used.
!    sg_znotcdef   Default apparent hydraulic roughness (m).
!    sg_znotdef    Default hydraulic roughness (m).
!
  integer, parameter   :: sg_n = 20

  real(sp), parameter  :: sg_pi = pi

  real(sp), parameter  :: sg_Cdmax = 0.01         ! non-dimensional
  real(sp), parameter  :: sg_alpha = 1.0          ! non-dimensional
  real(sp), parameter  :: sg_g = 9.81             ! (m/s2)
  real(sp), parameter  :: sg_kappa = 0.41         ! non-dimensional
  real(sp), parameter  :: sg_nu = 0.00000119      ! (m2/s)
  real(sp), parameter  :: sg_tol = 0.0001         ! non-dimensional
  real(sp), parameter  :: sg_ustarcdef = 0.01     ! (m/s)
  real(sp), parameter  :: sg_z100 = 1.0           ! (m)
  real(sp), parameter  :: sg_z1min = 0.20         ! (m)
  real(sp), parameter  :: sg_z1p   = 1.0          ! non-dimensional
  real(sp), parameter  :: sg_znotcdef = 0.01      ! (m)
  real(sp)             :: sg_znotdef              ! (m)

  complex  :: sg_mp



contains


subroutine init_bbl
!=======================================================================
!   initialize the Bottom Boundary Layer  Module                       !
!=======================================================================
  Use Input_Util

  implicit none

  character(len=80)    :: fname
  integer linenum,i,k1
  real(sp)           :: ftemp
  character(len=120) :: stemp
  logical             :: fexist

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: init_bbl" 

  fname = SEDIMENT_MODEL_FILE
  !ensure sediment parameter file exists
  if(fname == "none" .or. fname == "NONE" .or. fname == "None")then
    bblfile = "./"//trim(input_dir)//"/"//trim(casename)//'_sediment.inp'
  else
    bblfile = "./"//trim(input_dir)//"/"//trim(fname)
  endif  
  
  

  inquire(file=trim(bblfile),exist=fexist)
  if(.not.fexist)then
    write(*,*)'BBL parameter file: ',trim(sedfile),' does not exist'
    write(*,*)'stopping'
    call pstop
  end if

  !read in MB_BBL parameters and options 
  Call Get_Val(MB_BBL_USE,bblfile,'MB_BBL_USE',line=linenum,echo=.false.)  
  Call Get_Val(MB_CALC_ZNOT,bblfile,'MB_CALC_ZNOT',line=linenum,echo=.false.)  
  Call Get_Val(MB_CALC_UB,bblfile,'MB_CALC_UB',line=linenum,echo=.false.) 
  Call Get_Val(MB_Z0BIO,bblfile,'MB_Z0BIO',line=linenum,echo=.false.) 
  Call Get_Val(MB_Z0BL,bblfile,'MB_Z0BL',line=linenum,echo=.false.) 
  Call Get_Val(MB_Z0RIP,bblfile,'MB_Z0RIP',line=linenum,echo=.false.) 

  !read in SG_BBL parameters and options
  Call Get_Val(SG_BBL_USE,bblfile,'SG_BBL_USE',line=linenum,echo=.false.) 
  Call Get_Val(SG_CALC_ZNOT,bblfile,'SG_CALC_ZNOT',line=linenum,echo=.false.) 
  Call Get_Val(SG_CALC_UB,bblfile,'SG_CALC_UB',line=linenum,echo=.false.) 
  Call Get_Val(SG_LOGINT,bblfile,'SG_LOGINT',line=linenum,echo=.false.) 

  !read in SSW_BBL parameters and options
  Call Get_Val(SSW_BBL_USE,bblfile,'SSW_BBL_USE',line=linenum,echo=.false.) 
  Call Get_Val(SSW_CALC_ZNOT,bblfile,'SSW_CALC_ZNOT',line=linenum,echo=.false.) 
  Call Get_Val(SSW_LOGINT,bblfile,'SSW_LOGINT',line=linenum,echo=.false.) 
  Call Get_Val(SSW_CALC_UB,bblfile,'SSW_CALC_UB',line=linenum,echo=.false.) 
  Call Get_Val(SSW_FORM_DRAG_COR,bblfile,'SSW_FORM_DRAG_COR',line=linenum,echo=.false.) 
  Call Get_Val(SSW_ZOBIO,bblfile,'SSW_ZOBIO',line=linenum,echo=.false.)
  Call Get_Val(SSW_ZOBL,bblfile,'SSW_ZOBL ',line=linenum,echo=.false.)
  Call Get_Val(SSW_ZORIP,bblfile,'SSW_ZORIP',line=linenum,echo=.false.)
  Call Get_Val(SGWC,bblfile,'SGWC',line=linenum,echo=.false.)
  Call Get_Val(M94WC,bblfile,'M94WC',line=linenum,echo=.false.)
  Call Get_Val(GM82_RIPRUF,bblfile,'GM82_RIPRUF',line=linenum,echo=.false.)
  Call Get_Val(N92_RIPRUF,bblfile,'N92_RIPRUF',line=linenum,echo=.false.)
  Call Get_Val(R88_RIPRUF,bblfile,'R88_RIPRUF',line=linenum,echo=.false.)
  
  if(MB_BBL_USE)then
      if(msr)write(ipt,*)'!-----------------------------------------'
      if(msr)write(ipt,*)'!     BOTTOM BOUNDARY LAYER SETUP         '
      if(msr) write(ipt,*)'! mb_bbl        : active'
    if(MB_CALC_ZNOT)then
      if(msr)write(ipt,*)'! MB_CALC_ZNOT  : active'
    else
      if(msr)write(ipt,*)'! MB_CALC_ZNOT  : inactive'
    end if
    if(MB_CALC_UB)then
      if(msr)write(ipt,*)'! MB_CALC_UB    : active'
    else
      if(msr)write(ipt,*)'! MB_CALC_UB    : inactive'
    end if
    if(MB_Z0BIO)then
      if(msr)write(ipt,*)'! MB_Z0BIO      : active'
    else
      if(msr)write(ipt,*)'! MB_Z0BIO      : inactive'
    end if
    if(MB_Z0BL)then
      if(msr)write(ipt,*)'! MB_Z0BL       : active'
    else
      if(msr)write(ipt,*)'! MB_Z0BL       : inactive'
    end if
    if(MB_Z0RIP)then
      if(msr)write(ipt,*)'! MB_Z0RIP      : active'
    else
      if(msr)write(ipt,*)'! MB_Z0RIP      : inactive'
    end if
      if(msr)write(ipt,*)'!-----------------------------------------'
  end if

  if(SG_BBL_USE)then
      if(msr)write(ipt,*)'!-----------------------------------------'
      if(msr)write(ipt,*)'!     BOTTOM BOUNDARY LAYER SETUP         '
      if(msr)write(ipt,*)'! sg_bbl         : active'
    if(SG_CALC_ZNOT)then
      if(msr)write(ipt,*)'! SG_CALC_ZNOT   : active'
    else
      if(msr)write(ipt,*)'! SG_CALC_ZNOT   : inactive'
    end if
    if(SG_CALC_UB)then
      if(msr)write(ipt,*)'! SG_CALC_UB     : active'
    else
      if(msr)write(ipt,*)'! SG_CALC_UB     : inactive'
    end if
    if(SG_LOGINT)then
      if(msr)write(ipt,*)'! SG_LOGINT      : active'
    else
      if(msr)write(ipt,*)'! SG_LOGINT      : inactive'
    end if
      if(msr)write(ipt,*)'!-----------------------------------------'
  end if

  if(SSW_BBL_USE)then
      if(msr)write(ipt,*)'!-----------------------------------------'
      if(msr)write(ipt,*)'!     BOTTOM BOUNDARY LAYER SETUP         '
      if(msr)write(ipt,*)'! ssw_bbl        : active'
    if(SSW_CALC_ZNOT)then
      if(msr)write(ipt,*)'! SSW_CALC_ZNOT          : active'
    else
      if(msr)write(ipt,*)'! SSW_CALC_ZNOT          : inactive'
    end if
    if(SSW_LOGINT)then
      if(msr)write(ipt,*)'! SSW_LOGINT             : active'
    else
      if(msr)write(ipt,*)'! SSW_LOGINT             : inactive'
    end if
    if(SSW_CALC_UB)then
      if(msr)write(ipt,*)'! SSW_CALC_UB            : active'
    else
      if(msr)write(ipt,*)'! SSW_CALC_UB            : inactive'
    end if
    if(SSW_FORM_DRAG_COR)then
      if(msr)write(ipt,*)'! SSW_FORM_DRAG_COR      : active'
    else
      if(msr)write(ipt,*)'! SSW_FORM_DRAG_COR      : inactive'
    end if
    if(SSW_ZOBIO)then
      if(msr)write(ipt,*)'! SSW_ZOBIO              : active'
    else
      if(msr)write(ipt,*)'! SSW_ZOBIO              : inactive'
    end if
    if(SSW_ZOBL)then
      if(msr)write(ipt,*)'! SSW_ZOBL               : active'
    else
      if(msr)write(ipt,*)'! SSW_ZOBL               : inactive'
    end if
    if(SSW_ZORIP)then
      if(msr)write(ipt,*)'! SSW_ZORIP              : active'
    else
      if(msr)write(ipt,*)'! SSW_ZORIP              : inactive'
    end if
    if(SGWC)then
      if(msr)write(ipt,*)'! SGWC                   : active'
    else
      if(msr)write(ipt,*)'! SGWC                   : inactive'
    end if
    if(M94WC)then
      if(msr)write(ipt,*)'! M94WC                  : active'
    else
      if(msr)write(ipt,*)'! M94WC                  : inactive'
    end if
    if(GM82_RIPRUF)then
      if(msr)write(ipt,*)'! GM82_RIPRUF            : active'
    else
      if(msr)write(ipt,*)'! GM82_RIPRUF            : inactive'
    end if
    if(N92_RIPRUF)then
      if(msr)write(ipt,*)'! N92_RIPRUF             : active'
    else
      if(msr)write(ipt,*)'! N92_RIPRUF             : inactive'
    end if
    if(R88_RIPRUF)then
      if(msr)write(ipt,*)'! R88_RIPRUF             : active'
    else
      if(msr)write(ipt,*)'! R88_RIPRUF             : inactive'
    end if
      if(msr)write(ipt,*)'!-----------------------------------------'
  end if

  allocate( Iconv     (0: MT)  )
  allocate( Ubot      (0: MT)  )
  allocate( Vbot      (0: MT)  )
  allocate( Ur        (0: MT)  )
  allocate( Vr        (0: MT)  )
  allocate( bustrc    (0: MT)  )
  allocate( bvstrc    (0: MT)  )
  allocate( bustrw    (0: MT)  )
  allocate( bvstrw    (0: MT)  )
  allocate( bustrcwmax(0: MT)  )
  allocate( bvstrcwmax(0: MT)  )
  allocate( bustr     (0: MT)  )
  allocate( bvstr     (0: MT)  )
  allocate( taucwmax  (0: MT)  )
  !allocate( Dwave     (0: MT)  ) now allocated in mod_wave_current

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: init_bbl"

end subroutine init_bbl


subroutine sg_bbl
!=======================================================================
!                                                                      !
!  This routine compute bottom stresses based on the wave-current      !
!  algorithm and the ripple geometry and movable bed roughness.        ! 
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Styles, R. and S.M. glenn,  2000: Modeling stratified wave and      !
!    current bottom boundary layers in the continental shelf, JGR,     !
!    105, 24119-24139.                                                 !
!                                                                      !
!=======================================================================
integer :: Iter, i, j, k

real(SP), parameter :: eps = 1.0E-10

real(SP) :: Fwave_bot, Kb, Kbh, KboKb0, Kb0, Kdelta, Ustr
real(SP) :: anglec, anglew
real(SP) :: cff, cff1, cff2, cff3, og, fac, fac1, fac2
real(SP) :: sg_ab, sg_abokb, sg_a1, sg_b1, sg_chi, sg_c1, sg_dd
real(SP) :: sg_epsilon, sg_eta, sg_fofa, sg_fofb, sg_fofc, sg_fwm
real(SP) :: sg_kbs, sg_lambda, sg_mu, sg_phicw, sg_ro, sg_row
real(SP) :: sg_shdnrm, sg_shld, sg_shldcr, sg_scf, sg_ss, sg_star
real(SP) :: sg_ub, sg_ubokur, sg_ubouc, sg_ubouwm, sg_ur
real(SP) :: sg_ustarc, sg_ustarcw, sg_ustarwm, sg_znot, sg_znotp
real(SP) :: sg_zr, sg_zrozn, sg_z1, sg_z1ozn, sg_z2, twopi, zd1, zd2

real(sp), dimension(1:kbm1) :: Un,Vn

real(SP), dimension(0:MT) :: Ab
real(SP), dimension(0:MT) :: Tauc
real(SP), dimension(0:MT) :: Tauw
!real(SP), dimension(0:MT) :: Taucwmax
real(SP), dimension(0:MT) :: Ur_sg
real(SP), dimension(0:MT) :: Vr_sg
real(SP), dimension(0:MT) :: Ub
real(SP), dimension(0:MT) :: Ucur
real(SP), dimension(0:MT) :: Umag
real(SP), dimension(0:MT) :: Vcur
real(SP), dimension(0:MT) :: Zr
real(SP), dimension(0:MT) :: phic
real(SP), dimension(0:MT) :: phicw
real(SP), dimension(0:MT) :: rheight
real(SP), dimension(0:MT) :: rlength
real(SP), dimension(0:MT) :: u100
real(SP), dimension(0:MT) :: znot
real(SP), dimension(0:MT) :: znotc

integer  :: cnt,jj,cc

logical :: iterate

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: sg_bbl"

!
!-----------------------------------------------------------------------
!  Initalize to default values.
!-----------------------------------------------------------------------
!
 sg_mp=CMPLX(SQRT(1.0/(2.0*sg_z1p)),SQRT(1.0/(2.0*sg_z1p)))
 sg_znotdef=BOTTOM_ROUGHNESS_LENGTHSCALE

 do i=1,m
     Tauc(i)=0.0
     Tauw(i)=0.0
     Taucwmax(i)=0.0
     u100(i)=0.0
     rheight(i)=0.0
     rlength(i)=0.0
     znot(i)=sg_znotdef
     znotc(i)=0.0
  end do
!
!-----------------------------------------------------------------------
!  Set currents above bed.
!-----------------------------------------------------------------------
!
  do i=1,m
     Zr(i)=dt(i)*(zz(i,kbm1)-z(i,kb))
     Ur_sg(i) = 0.0
     Vr_sg(i) = 0.0
     Un       = 0.0
     Vn       = 0.0
     cnt = 0 
     do jj=1,ntve(i)
        cnt =cnt + 1
        cc = nbve(i,jj)
        Ur_sg(i) = Ur_sg(i) + u(cc,kbm1)
        Vr_sg(i) = Vr_sg(i) + v(cc,kbm1)
        Un(:)    = Un(:)    + u(cc,:)
        Vn(:)    = Vn(:)    + v(cc,:)
     end do
     Ur_sg(i) = Ur_sg(i) / cnt
     Vr_sg(i) = Vr_sg(i) / cnt
     Un(:)    = Un(:)    / cnt
     Vn(:)    = Vn(:)    / cnt
     
     if(SG_LOGINT)then
!
!  If current height is less than z1ur, interpolate logarithmically
!  to z1ur.
!
        IF (Zr(i).lt.sg_z1min) THEN
          DO k=kbm1-1,1,-1
            zd1=dt(i)*(zz(i,k+1)-z(i,kb)) !!!??? z_r(i,j,k-1)-z_w(i,j,0)
            zd2=dt(i)*(zz(i,k)  -z(i,kb)) !!!??? z_r(i,j,k  )-z_w(i,j,0)
            IF ((zd1.lt.sg_z1min).and.(sg_z1min.lt.zd2)) THEN
              fac=1.0/LOG(zd2/zd1)
              fac1=fac*LOG(zd2/sg_z1min)
              fac2=fac*LOG(sg_z1min/zd1)
              Ur_sg(i)=fac1*Un(k+1)+fac2*Un(k)
              Vr_sg(i)=fac1*Vn(k+1)+fac2*Vn(k)
              Zr(i)=sg_z1min
            END IF
          END DO
        END IF
     end if

  end do
!
!-----------------------------------------------------------------------
!  Compute bed wave orbital velocity (m/s) and excursion amplitude
!  (m) from wind-induced waves.  Use linear wave theory dispersion
!  relation for wave number.
!-----------------------------------------------------------------------
!
   twopi=2.0*pi
   og=1.0/grav

   do i=1,m
!
!  Compute first guess for wavenumber, Kb0.  Use deep water (Kb0*h>1)
!  and shallow water (Kb0*H<1) approximations.
!
          Fwave_bot=twopi/MAX(Pwave_bot(i),0.05)
!
!  Compute bed wave orbital velocity and excursion amplitude.
!
     if(SG_CALC_UB)then
          Kb0=Fwave_bot*Fwave_bot*og
          IF (Kb0*h(i).ge.1.0) THEN
            Kb=Kb0
          ELSE
            Kb=Fwave_bot/SQRT(grav*h(i))
          END IF
!
!  Compute bottom wave number via Newton-Raphson method.
!
          ITERATE=.TRUE.
          DO Iter=1,sg_n
            IF (ITERATE) THEN
              Kbh=Kb*h(i)
              KboKb0=Kb/Kb0
              Kdelta=(1.0-KboKb0*TANH(Kbh))/(1.0+Kbh*(KboKb0-1.0/KboKb0))
              ITERATE=ABS(Kb*Kdelta) .ge. sg_tol
              Kb=Kb*(1.0+Kdelta)
            END IF
          END DO
          Ab(i)=0.5*HSC1(i)/SINH(Kb*h(i))+eps
          Ub(i)=Fwave_bot*Ab(i)+eps
     else
          Ub(i)=ABS(Ub_swan(i))+eps
          Ab(i)=Ub(i)/Fwave_bot+eps
     endif
!
!  Compute bottom current magnitude at RHO-points.
!
          Ucur(i)=Ur_sg(i)
          Vcur(i)=Vr_sg(i)
          Umag(i)=SQRT(Ucur(i)*Ucur(i)+Vcur(i)*Vcur(i))+eps
!
!  Compute angle between currents and waves (radians)
!
          IF (Ucur(i).eq.0.0) THEN
            phic(i)=0.5*pi*SIGN(1.0,Vcur(i))
          ELSE
            phic(i)=ATAN2(Vcur(i),Ucur(i))
          ENDIF
          !phicw(i)=Dwave(i)-phic(i)  
          phicw(i)=1.5_sp*pi-Dwave(i)-phic(i) !-angler(i,j)
    end do
!
!-----------------------------------------------------------------------
!  Set default logarithmic profile.
!-----------------------------------------------------------------------
!
    do i=1,m
          IF (Umag(i).gt.0.0) THEN
!!          Ustr=MIN(sg_ustarcdef,Umag(i,j)*vonKar/                     &
!!   &                            LOG(Zr(i,j)/sg_znotdef))
!!          Tauc(i,j)=Ustr*Ustr
            cff1=vonKar/LOG(Zr(i)/sg_znotdef)
            cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
            Tauc(i)=cff2*Umag(i)*Umag(i)
          END IF
    end do

!
!-----------------------------------------------------------------------
!  Wave-current interaction case.
!-----------------------------------------------------------------------
!
   do i=1,m
          sg_dd=bottom(i,isd50)
          sg_ss=bottom(i,idens)/(rho(i,kbm1)+1000.0)
          sg_ab=Ab(i)
          sg_ub=Ub(i)
          sg_phicw=phicw(i)
          sg_ur=Umag(i)
          sg_zr=Zr(i)
!
!  Compute hydraulic roughness "Znot" (m), ripple height "eta" (m),
!  and ripple length "lambda" (m).
!
     if(SG_CALC_ZNOT)then
          sg_star=sg_dd/(4.0*sg_nu)*SQRT((sg_ss-1.0)*sg_g*sg_dd)
!
!  Compute critical shield parameter based on grain diameter.
!  (sg_scf is a correction factor).
!
          sg_scf=1.0
          IF (sg_star.le.1.5) THEN
            sg_shldcr=sg_scf*0.0932*sg_star**(-0.707)
          ELSE IF ((1.5.lt.sg_star).and.(sg_star.lt.4.0)) THEN
            sg_shldcr=sg_scf*0.0848*sg_star**(-0.473)
          ELSE IF ((4.0.le.sg_star).and.(sg_star.lt.10.0)) THEN
            sg_shldcr=sg_scf*0.0680*sg_star**(-0.314)
          ELSE IF ((10.0.le.sg_star).and.(sg_star.lt.34.0)) THEN
            sg_shldcr=sg_scf*0.033
          ELSE IF ((34.0.le.sg_star).and.(sg_star.lt.270.0)) THEN
            sg_shldcr=sg_scf*0.0134*sg_star**(0.255)
          ELSE
            sg_shldcr=sg_scf*0.056
          END IF
!
!  Calculate skin friction shear stress based on Ole Madsen (1994)
!  empirical formula. Check initiation of sediment motion criteria,
!  to see if we compute sg_znot based on the wave-formed ripples.
!  If the skin friction calculation indicates that sediment is NOT
!  in motion, the ripple model is invalid and take the default value,
!  sg_znotdef.
!
          sg_abokb=sg_ab/sg_dd
          IF (sg_abokb.le.100.0) THEN
            sg_fwm=EXP(7.02*sg_abokb**(-0.078)-8.82)
          ELSE
            sg_fwm=EXP(5.61*sg_abokb**(-0.109)-7.30)
          END IF
          sg_ustarwm=SQRT(0.5*sg_fwm)*sg_ub
          sg_shdnrm=(sg_ss-1.0)*sg_dd*sg_g
          sg_shld=sg_ustarwm*sg_ustarwm/sg_shdnrm
          IF ((sg_shld/sg_shldcr).le.1.0) THEN
            sg_znot=sg_znotdef
            sg_eta=0.0
            sg_lambda=0.0
          ELSE
!
!  Calculate ripple height and length and bottom roughness
!
            sg_chi=4.0*sg_nu*sg_ub*sg_ub/                            &
     &             (sg_dd*((sg_ss-1.0)*sg_g*sg_dd)**1.5)
            IF (sg_chi.le.2.0) THEN
              sg_eta=sg_ab*0.30*sg_chi**(-0.39)
              sg_lambda=sg_ab*1.96*sg_chi**(-0.28)
            ELSE
              sg_eta=sg_ab*0.45*sg_chi**(-0.99)
              sg_lambda=sg_ab*2.71*sg_chi**(-0.75)
            END IF
            sg_kbs=sg_ab*0.0655*                                     &
     &             (sg_ub*sg_ub/((sg_ss-1.0)*sg_g*sg_ab))**1.4
            sg_znot=(sg_dd+2.3*sg_eta+sg_kbs)/30.0
          END IF
     else
          sg_znot=sg_znotdef
          sg_chi=4.0*sg_nu*sg_ub*sg_ub/                              &
     &           (sg_dd*((sg_ss-1.0)*sg_g*sg_dd)**1.5)
          IF (sg_chi.le.2.0) THEN
            sg_eta=sg_ab*0.32*sg_chi**(-0.34)
            sg_lambda=sg_ab*2.04*sg_chi**(-0.23)
          ELSE
            sg_eta=sg_ab*0.52*sg_chi**(-1.01)
            sg_lambda=sg_ab*2.7*sg_chi**(-0.78)
          END IF
     endif
          znot(i)=sg_znot
          rheight(i)=sg_eta
          rlength(i)=sg_lambda
!
!  Compute only when nonzero currents and waves.
!
          sg_zrozn=sg_zr/sg_znot
          IF ((sg_ur.gt.0.0).and.(sg_ub.gt.0.0).and.              &
     &        (sg_zrozn.gt.1.0)) THEN
!
! Compute bottom stress based on ripple roughness.
!
            sg_ubokur=sg_ub/(sg_kappa*sg_ur)
            sg_row=sg_ab/sg_znot
            sg_a1=1.0E-6
            CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur,     &
     &                       sg_a1, sg_mu, sg_epsilon, sg_ro, sg_fofa)
            sg_abokb=sg_ab/(30.0*sg_znot)
            IF (sg_abokb.le.100.0) THEN
              sg_fwm=EXP(-8.82+7.02*sg_abokb**(-0.078))
            ELSE
              sg_fwm=EXP(-7.30+5.61*sg_abokb**(-0.109))
            END IF
            sg_ubouwm=SQRT(2.0/sg_fwm)
!
!  Determine the maximum ratio of wave over combined shear stresses,
!  sg_ubouwm (ub/ustarwm).
!
            CALL sg_purewave (sg_row, sg_ubouwm, sg_znotp, sg_ro)
!
!  Set initial guess of the ratio of wave over shear stress, sg_c1
!  (ub/ustarc).
!
            sg_b1=sg_ubouwm
            sg_fofb=-sg_fofa
            sg_c1=0.5*(sg_a1+sg_b1)
            CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur,     &
     &                       sg_c1, sg_mu, sg_epsilon, sg_ro, sg_fofc)
!
!  Solve PDE via bi-section method.
!
            ITERATE=.TRUE.
            DO Iter=1,sg_n
              IF (ITERATE) THEN
                IF ((sg_fofb*sg_fofc).lt.0.0) THEN
                  sg_a1=sg_c1
                ELSE
                  sg_b1=sg_c1
                END IF
                sg_c1=0.5*(sg_a1+sg_b1)
                CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur, &
     &                           sg_c1, sg_mu, sg_epsilon, sg_ro,       &
     &                           sg_fofc)
                ITERATE=(sg_b1-sg_c1) .ge. sg_tol
                IF (ITERATE) Iconv(i)=Iter
              END IF
            END DO
            sg_ubouc=sg_c1
!
!  Compute bottom shear stress magnitude (m/s).
!
            sg_ustarcw=sg_ub/sg_ubouc
            sg_ustarwm=sg_mu*sg_ustarcw
!!          sg_ustarc=MIN(sg_ustarcdef,sg_epsilon*sg_ustarcw)
!!          sg_ustarc=MIN(SQRT(Tauc(i,j)),sg_epsilon*sg_ustarcw)
            sg_ustarc=MAX(SQRT(Tauc(i)),sg_epsilon*sg_ustarcw)
            Tauc(i)=sg_ustarc*sg_ustarc
            Tauw(i)=sg_ustarwm*sg_ustarwm
            Taucwmax(i)=SQRT((Tauc(i)+                            &
     &                          Tauw(i)*COS(phicw(i)))**2+          &
     &                         (Tauw(i)*SIN(phicw(i)))**2)
!
!  Compute apparent hydraulic roughness (m).
!
            IF (sg_epsilon.gt.0.0) THEN
              sg_z1=sg_alpha*sg_kappa*sg_ab/sg_ubouc
              sg_z2=sg_z1/sg_epsilon
              sg_z1ozn=sg_z1/sg_znot
              znotc(i)=sg_z2*EXP(-(1.0-sg_epsilon+sg_epsilon*LOG(sg_z1ozn)))
!
!  Compute mean (m/s) current at 100 cm above the bottom.
!
              IF (sg_z100.gt.sg_z2) THEN
                u100(i)=sg_ustarc*                                 &
     &                    (LOG(sg_z100/sg_z2)+1.0-sg_epsilon+        &
     &                    sg_epsilon*LOG(sg_z1ozn))/sg_kappa
              ELSE IF ((sg_z100.le.sg_z2).and.(sg_zr.gt.sg_z1)) THEN
                u100(i)=sg_ustarc*sg_epsilon*                       &
     &                    (sg_z100/sg_z1-1.0+LOG(sg_z1ozn))/sg_kappa
              ELSE
                u100(i)=sg_ustarc*sg_epsilon*                       &
     &                    LOG(sg_z100/sg_znot)/sg_kappa
              END IF
            END IF
          END IF
  end do

!
!-----------------------------------------------------------------------
!  Compute kinematic bottom stress components due current and wind-
!  induced waves.
!-----------------------------------------------------------------------
!
     do i=1,m
          anglec=Ur_sg(i)/Umag(i)
          bustr(i)=Tauc(i)*anglec
          anglec=Vr_sg(i)/Umag(i)
          bvstr(i)=Tauc(i)*anglec
     end do
     
      do i=1,m
          anglec=Ucur(i)/Umag(i)
          anglew=COS(Dwave(i))
          bustrc(i)=Tauc(i)*anglec
          bustrw(i)=Tauw(i)*anglew
          bustrcwmax(i)=Taucwmax(i)*anglew
          Ubot(i)=Ub(i)*anglew
          Ur(i)=Ucur(i)
!
          anglec=Vcur(i)/Umag(i)
          anglew=SIN(Dwave(i))
          bvstrc(i)=Tauc(i)*anglec
          bvstrw(i)=Tauw(i)*anglew
          bvstrcwmax(i)=Taucwmax(i)*anglew
          Vbot(i)=Ub(i)*anglew
          Vr(i)=Vcur(i)
!
          bottom(i,ibwav)=Ab(i)
          bottom(i,irhgt)=rheight(i)
          bottom(i,irlen)=rlength(i)
          bottom(i,izdef)=znot(i)
          bottom(i,izapp)=znotc(i)

# if defined (WET_DRY)
          if(iswetn(i)/=1)then
             bustr(i)       = 0.0_sp
             bvstr(i)       = 0.0_sp
             bustrc(i)      = 0.0_sp
             bvstrc(i)      = 0.0_sp
             bustrw(i)      = 0.0_sp
             bvstrw(i)      = 0.0_sp
             bustrcwmax(i)  = 0.0_sp
             bvstrcwmax(i)  = 0.0_sp
             Ubot(i)        = 0.0_sp
             Vbot(i)        = 0.0_sp
             Ur(i)          = 0.0_sp
             Vr(i)          = 0.0_sp
             taucwmax(i)    = 0.0_sp
          end if
# endif	  
     end do


!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.

# if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustr     , bvstr      )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustrc    , bvstrc     )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustrw    , bvstrw     )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustrcwmax, bvstrcwmax )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Ubot      , Vbot       )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Ur        , Vr         )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,TauCWmax               ) 
# if defined (ORIG_SED)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bottom                 )
# elif defined (CSTMS_SED)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,sedbed%bottom                 )
# endif
# endif

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: sg_bbl"

    return

end subroutine sg_bbl




subroutine mb_bbl
!=======================================================================
!                                                                      !
!  This subroutine computes bottom stresses for combined waves and     !
!  currents using the parametric approximation by Soulsby 1997:        !
!                                                                      !
!        tauCW=tauC*[1+1.2*(tauW/(tauC+tauW))^3.2]                     !
!                                                                      !
!  and                                                                 !
!                                                                      !
!        tauCWmax=SQRT([tauCW+tauW*COS(phiCW)]^2+[tauW*SIN(phiCW)]^2)  !
!                                                                      !
!  where                                                               !
!                                                                      !
!  tauCW      Combined wave-averaged stress (in current direction).    !
!  tauC       Stress due to currents, if waves would be absent.        ! 
!  tauW       Amplitude of stress due to waves without currents.       !
!  tauCWmax   Maximum combined wave-averaged stress.                   !
!  phiCW      Angle between current and waves.                         !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Dyer 1986, Coastal and Estuarine Sediment Dynamics, Wiley,        !
!      342 pp.                                                         !
!                                                                      !
!    Harris and Wiberg 2001, Comp. and Geosci. 27, 675-690.            !
!                                                                      !
!    Li and Amos 2001, Comp. and Geosci. 27, 619-645.                  !
!                                                                      !
!    Soulsby 1997, Dynamics of Marine Sands, Telford  Publ., 249 pp.   !
!                                                                      !
!    Soulsby 1995, Bed shear-stresses due to combined waves and        !
!      currents, in: Stive et al: Advances in Coastal Morphodynamics,  !
!      Wiley, 4.20-4.23                                                !
!                                                                      !
!    Wiberg and Harris 1994, J. Geophys. Res. 99(C4), 775-789.         !
!                                                                      !
!=======================================================================
  integer :: i, ised, j
  integer :: jj,cnt,cc

  real(sp), parameter :: K1 = 0.6666666666
  real(sp), parameter :: K2 = 0.3555555555
  real(sp), parameter :: K3 = 0.1608465608
  real(sp), parameter :: K4 = 0.0632098765
  real(sp), parameter :: K5 = 0.0217540484
  real(sp), parameter :: K6 = 0.0065407983

  real(sp), parameter :: eps = 1.0E-10

  real(sp), parameter :: scf1 = 0.5 * 1.39
  real(sp), parameter :: scf2 = 0.52
  real(sp), parameter :: scf3 = 2.0 - scf2
  real(sp), parameter :: scf4 = 1.2
  real(sp), parameter :: scf5 = 3.2

  real(sp) :: twopi = 2.0 * pi

  real(sp) :: RHbiomax = 0.006    ! maximum biogenic ripple height
  real(sp) :: RHmin = 0.001       ! arbitrary minimum ripple height
  real(sp) :: RLmin = 0.01        ! arbitrary minimum ripple length

  real(sp) :: Ab, Fwave_bot, Kbh, Kbh2, Kdh
  real(sp) :: angleC, angleW, phiC, phiCW
  real(sp) :: cff, cff1, cff2, d50, viscosity, wset
  real(sp) :: RHbio, RHbiofac, RHmax, RLbio
  real(sp) :: rhoWater, rhoSed
  real(sp) :: rhgt, rlen, thetw
  real(sp) :: Znot, ZnotC, Znot_bl, Znot_rip
  real(sp) :: tau_bf, tau_ex,  tau_en, tau_up, tau_w, tau_wb
  real(sp) :: tau_c, tau_cb, tau_cs, tau_cw, tau_cwb, tau_cws

  real(sp), dimension(0:mt) :: Ub
  real(sp), dimension(0:mt) :: Ucur
  real(sp), dimension(0:mt) :: Vcur
  real(sp), dimension(0:mt) :: Zr
  real(sp), dimension(0:mt) :: Ur_mb
  real(sp), dimension(0:mt) :: Vr_mb
  real(sp), dimension(0:mt) :: Umag
  real(sp), dimension(0:mt) :: tauC
  real(sp), dimension(0:mt) :: tauW
  real(sp), dimension(0:mt) :: tauCW
!  real(sp), dimension(0:mt) :: tauCWmax


  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: mb_bbl"

!
!-----------------------------------------------------------------------
!  Initalize stresses due to currents and waves.
!-----------------------------------------------------------------------
!
  RHbiofac=1.0/EXP(4.11)
  do i=1,m
     tauC(i)=0.0
     tauCW(i)=0.0
     tauCWmax(i)=0.0
  end do
!
!-----------------------------------------------------------------------
!  Set currents above bed.
!-----------------------------------------------------------------------
!
  do i=1,m
     Zr(i)=dt(i)*(zz(i,kbm1)-z(i,kb))
     Ur_mb(i) = 0.0
     Vr_mb(i) = 0.0
     cnt = 0 
     do jj=1,ntve(i)
        cnt =cnt + 1
        cc = nbve(i,jj)
        Ur_mb(i) = Ur_mb(i) + u(cc,kbm1)
        Vr_mb(i) = Vr_mb(i) + v(cc,kbm1)
     end do
     Ur_mb(i) = Ur_mb(i) / cnt
     Vr_mb(i) = Vr_mb(i) / cnt
  end do
!
!=======================================================================
!  Compute bottom stresses.
!=======================================================================
!
  do i=1,m
          RHbio=0.0
          RLbio=0.0
          rlen=bottom(i,irlen)
          rhgt=bottom(i,irhgt)
          rhoWater=rho(i,kbm1)+1000.0
          viscosity=0.0013/rhoWater         ! kinematic viscosity
!
!-----------------------------------------------------------------------
!  Compute bed wave orbital velocity (m/s) and excursion amplitude (m)
!  from wind-induced waves.  Use Dean and Dalrymple (1991) 6th-degree
!  polynomial to approximate wave number on shoaling water.
!-----------------------------------------------------------------------
!
          Fwave_bot=twopi/MAX(Pwave_bot(i),0.05)
          if(MB_CALC_UB)then
             Kdh=h(i)*Fwave_bot*Fwave_bot/grav
             Kbh2=Kdh*Kdh+                                              &
     &         Kdh/(1.0+Kdh*(K1+Kdh*(K2+Kdh*(K3+Kdh*(K4+             &
     &              Kdh*(K5+K6*Kdh))))))
             Kbh=SQRT(Kbh2)
!
!  Compute bed wave orbital velocity and excursion amplitude.
!
             Ab=0.5*HSC1(i)/SINH(Kbh)+eps
             Ub(i)=Fwave_bot*Ab
          else
             Ub(i)=Ub_swan(i)
             Ab=Ub(i)/Fwave_bot+eps
          endif
!
!  Compute bottom current magnitude at RHO-points.
!
          Ucur(i)=Ur_mb(i)
          Vcur(i)=Vr_mb(i)
          Umag(i)=SQRT(Ucur(i)**2+Vcur(i)**2)+eps
!
!  Compute angle between currents and waves (radians)
!
          IF (Ucur(i).ne.0.0) THEN
            phiC=ATAN2(Vcur(i),Ucur(i))
          ELSE
            phiC=0.5*pi*SIGN(1.0,Vcur(i))
          END IF
          phiCW=Dwave(i)-phiC
!
!-----------------------------------------------------------------------
!  Determine skin roughness from sediment size. Set default logarithmic
!  profile consistent with current-only case.
!  Establish local median grain size for all calculations here and
!  determine local values of critical stresses.
!  Since most parameterizations have been derived ignoring multiple
!  grain sizes, we apply this single d50 also in the case of mixed beds. 
!-----------------------------------------------------------------------
!
          d50=bottom(i,isd50)              ! (m)
          tau_cb=bottom(i,itauc)/1025.     ! (m2/s2)
          wset=bottom(i,iwset)             ! (m/s)
          rhoSed=bottom(i,idens)/rhoWater  ! (nondimensional)
!
! Critical stress (m2/s2) for transition to sheet flow
! (Li and Amos, 2001, eq. 11).
!
          tau_up=0.172*(rhoSed-1.0)*grav*(d50**0.624)
!
! Threshold stress (m2/s2) for break off (Grant and Madsen,1982).
!
          tau_bf=0.79*(viscosity**(-0.6))*                        &
     &           (((rhoSed-1.0)*grav)**0.3)*(d50**0.9)*tau_cb
!
!  Set Znot for currents as maximun of user value or grain roughness.
!
           ZnotC=d50/12.0
           Znot=MAX(BOTTOM_ROUGHNESS_LENGTHSCALE,ZnotC) 
!
!-----------------------------------------------------------------------
!  Set tauC (m2/s2) acc. to current-only case (skin friction) [m/s]
!-----------------------------------------------------------------------
!
          cff1=vonKar/LOG(Zr(i)/Znot)
          cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
          tauC(i)=cff2*Umag(i)*Umag(i)

          cff1=vonKar/LOG(Zr(i)/ZnotC)
          tau_cs=cff1*cff1*Umag(i)*Umag(i)
!
!-----------------------------------------------------------------------
!  If significant waves (Ub > 0.01 m/s), wave-current interaction case
!  according to Soulsby (1995). Otherwise, tauCWmax = tauC for sediment
!  purposes.
!-----------------------------------------------------------------------
!
          IF (Ub(i).gt.0.01) THEN  
!            
!  Determine skin stresses (m2/s2) for pure waves and combined flow
!  using Soulsby (1995) approximation of the wave friction factor,
!  fw=2*scf1*(Znot/Ab)**scf2 so tauW=0.5fw*Ub**2.
!
            tau_w=scf1*((ZnotC*Fwave_bot)**scf2)*(Ub(i)**scf3)
!      
!  Wave-averaged, combined wave-current stress.(Eqn. 69, Soulsby, 1997).
!
            tau_cw=tau_cs*                                              &
     &             (1.0+scf4*((tau_w/(tau_w+tau_cs))**scf5))
!      
!  Maximum of combined wave-current skin stress (m2/s2) component for
!  sediment.(Eqn. 70, Soulsby, 1997).
!
            tau_cws=SQRT((tau_cw+tau_w*COS(phiCW))**2+                  &
     &                   (tau_w*SIN(phiCW))**2)
!!          tauw(i,j)=tau_cws
            tauCWmax(i)=tau_cws
            tauW(i)=tau_w
!
!  Set combined stress for Znot.
!
            tau_w=scf1*((Znot*Fwave_bot)**scf2)*(Ub(i)**scf3)
            tau_cw=tauC(i)*                                           &
     &             (1.0+scf4*((tau_w/(tau_w+tauC(i)))**scf5))
                           
            if(MB_Z0BL)then
!
!-----------------------------------------------------------------------
!  Compute bedload roughness for ripple predictor and sediment purposes.
!  At high transport stages, friction depends on thickness of bedload
!  layer. Shear stress due to combined grain and bedload roughness is
!  used to predict ripples and onset of suspension (Li and Amos, 2001).
!-----------------------------------------------------------------------
!
               if(MB_CALC_ZNOT)then
                  tau_ex=MAX((tau_cws-tau_cb),0.0)
                  cff=(1.0/((rhoSed-1.0)*grav*d50))
                  Znot_bl=17.4*d50*(cff*tau_ex)**0.75
                  ZnotC=ZnotC+Znot_bl
               endif
!     
!-----------------------------------------------------------------------
!  Compute stresses (m2/s2)for sediment purposes, using grain and
!  bedload roughness.
!-----------------------------------------------------------------------
!
              cff1=vonKar/LOG(Zr(i)/ZnotC)
              tau_c=cff1*cff1*Umag(i)*Umag(i)
              tau_wb=scf1*((ZnotC*Fwave_bot)**scf2)*(Ub(i)**scf3)
              tau_cw=tau_c*(1.0+scf4*((tau_wb/(tau_wb+tau_c))**scf5))
!            
!  Maximum of combined wave-current stress (m2/s2) component for
!  sediment purposes.
!
              tau_cwb=SQRT((tau_cw+tau_wb*COS(phiCW))**2+                  &
     &                   (tau_wb*SIN(phiCW))**2)
              tauCWmax(i)=tau_cwb
              tauW(i)=tau_wb
            endif

            if(MB_Z0RIP)then
!
!-----------------------------------------------------------------------
!  Determine bedform roughness ripple height (m) and ripple length (m)
!  for sandy beds. Use structure according to Li and Amos (2001).
!-----------------------------------------------------------------------
!
!  Check median grain diameter
!
              IF (d50.ge.0.000063) THEN
!      
!  Enhanced skin stress if pre-exisiting ripples (Nielsen, 1986).
!
                RHmax=0.8*rlen/pi
                rhgt=MAX(MIN(RHmax,rhgt),RHmin)
                tau_en=MAX(tau_cws,tau_cws*(rlen/(rlen-pi*rhgt))**2) 
!
                IF ((tau_cws.lt.tau_cb).and.(tau_en.ge.tau_cb)) THEN 
                   rhgt=(19.6*(SQRT(tau_cws/tau_cb))+20.9)*d50
                   rlen=rhgt/0.12        ! local transport
                ELSE 
                  IF ((tau_cws.ge.tau_cb).and.(tau_cwb.lt.tau_bf)) THEN 
                    rhgt=(22.15*(SQRT(tau_cwb/tau_cb))+6.38)*d50
                    rlen=rhgt/0.12       ! bed load in eq. range
                  ELSE 
                    IF ((tau_cwb.ge.tau_bf).and.(tau_cwb.lt.tau_up)) THEN 
                      rlen=535.0*d50     ! break off regime
                      rhgt=0.15*rlen*(SQRT(tau_up)-SQRT(tau_cwb))/     &
     &                                (SQRT(tau_up)-SQRT(tau_bf )) 
                    ELSE IF (tau_cwb.ge.tau_up) THEN 
                      rlen=0.0          ! sheet flow, plane bed
                      rhgt=0.0
                    ELSE 
                      rlen=bottom(i,irlen)
                      rhgt=bottom(i,irhgt)
                    END IF
                  END IF
                END IF
              END IF
            endif

            if(MB_Z0BIO)then
!
!-----------------------------------------------------------------------
!  Determine (biogenic) bedform roughness ripple height (m) and ripple
!  length (m) for silty beds, using Harris and Wiberg (2001).
!-----------------------------------------------------------------------
!
!  Use 10 cm default biogenic ripple length, RLbio (Wheatcroft 1994).
!
!  NOTE:  For mixed beds take average of silty and sandy bedform
!         roughnesses weighted according to silt and sand fractions.
!
               IF (d50.lt.0.000063) THEN
                 RLbio=0.1
                 thetw=tau_cws*(1.0/((rhoSed-1.0)*grav*d50))
                 RHbio=(thetw**(-1.67))*RLbio*RHbiofac
                 rhgt=MIN(RHbio,RHbiomax)
                 rlen=RLbio
               END IF
            endif

            if(MB_Z0RIP.or. MB_Z0BIO)then
!
!  Ripple roughness using Grant and Madsen (1982) roughness length.
!
                if(MB_CALC_ZNOT)then
                   Znot_rip=0.92*rhgt*rhgt/(MAX(rlen,RLmin))
                   ZnotC=ZnotC+Znot_rip
                endif
!
!-----------------------------------------------------------------------
! Compute bottom stress (m2/s2) components based on total roughnes.
!-----------------------------------------------------------------------
!
                cff1=vonKar/LOG(Zr(i)/ZnotC)
                tau_c=cff1*cff1*Umag(i)*Umag(i)
                tau_w=scf1*((ZnotC*Fwave_bot)**scf2)*(Ub(i)**scf3)
                tau_cw=tau_c*                                               &
     &               (1.0+scf4*((tau_w/(tau_w+tau_c))**scf5))
!!          tau_cwb=SQRT((tau_cw+tau_w*COS(phiCW))**2+                  &
!!   &                   (tau_w*SIN(phiCW))**2)
!!          tauCWmax(i,j)=tau_cwb
            endif
!
!-----------------------------------------------------------------------
!  Compute effective bottom shear velocity (m/s) relevant for flow and
!  eddy-diffusivities/viscosity.
!-----------------------------------------------------------------------
!
!!          tauC(i,j)=tau_cw
            tauCW(i)=tau_cw
            tauW(i)=tau_w
          ELSE IF (Ub(i).le.0.01) THEN   
!
!  If current-only, tauCWmax=tauC(skin) for use in sediment model; tauC
!  is determined using roughness due to current ripples.
!  In this limit, tauCWmax=tauCW=tauC
!            
            tauCWmax(i)=tauC(i)
            tauW(i)=0.0
!!          tauW(i,j)=tau_cs
!!          tauCW(i,j)=tau_cs

            if(MB_Z0RIP)then    
!!             IF (tauC(i,j).gt.tau_up) THEN
               IF (tauC(i).gt.tau_up) THEN
                 rhgt=0.0
                 rlen=0.0
!!             ELSE IF (tauC(i,j).lt.tau_cb) THEN
               ELSE IF (tauC(i).lt.tau_cb) THEN
                 rlen=bottom(i,irlen)
                 rhgt=bottom(i,irhgt)
               ELSE
                 rlen=1000.0*d50                       ! Yalin (1964)
                 rhgt=0.0308*(rlen**1.19)
!!               rhgt=100.0_r8*0.074_r8*                                   &
!!   &                (0.01_r8*rlen)**1.19_r8             ! Allen (1970)
               END IF
               if(MB_CALC_ZNOT)then
                 ZnotC=ZnotC+0.92*rhgt*rhgt/(MAX(rlen,RLmin))
               endif
            endif 
            cff1=vonKar/LOG(Zr(i)/ZnotC)
            cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
!!          tauc(i,j)=cff2*Umag(i,j)*Umag(i,j)
            tauCW(i)=cff2*Umag(i)*Umag(i)
          END IF 
!
!-----------------------------------------------------------------------
!  Load variables for output purposes.
!-----------------------------------------------------------------------
!
          bottom(i,ibwav)=Ab
          bottom(i,irlen)=rlen
          bottom(i,irhgt)=rhgt
          bottom(i,izdef)=Znot                       ! Zob(ng)
          bottom(i,izapp)=ZnotC
   end do
!
!
!-----------------------------------------------------------------------
!  Compute kinematic bottom stress (m2/s2) components for flow due to
!  combined current and wind-induced waves.
!-----------------------------------------------------------------------
!
   do i=1,m
          angleC=Ur_mb(i)/Umag(i)
          bustr(i)=tauCW(i)*angleC
          angleC=Vr_mb(i)/Umag(i)
          bvstr(i)=tauCW(i)*angleC
   end do

   do i=1,m
          angleC=Ucur(i)/Umag(i)
          angleW=COS(Dwave(i))
          bustrc(i)=tauCW(i)*angleC
          bustrw(i)=tauW(i)*angleW
          bustrcwmax(i)=tauCWmax(i)*angleW
          Ubot(i)=Ub(i)*angleW
          Ur(i)=Ucur(i)
!
          angleC=Vcur(i)/Umag(i)
          angleW=SIN(Dwave(i))
          bvstrc(i)=tauCW(i)*angleC
          bvstrw(i)=tauW(i)*angleW
          bvstrcwmax(i)=tauCWmax(i)*angleW
          Vbot(i)=Ub(i)*angleW
          Vr(i)=Vcur(i)
	  
# if defined (WET_DRY)	  
          if(iswetn(i)/=1)then
             bustr(i)       = 0.0_sp
             bvstr(i)       = 0.0_sp
             bustrc(i)      = 0.0_sp
             bvstrc(i)      = 0.0_sp
             bustrw(i)      = 0.0_sp
             bvstrw(i)      = 0.0_sp
             bustrcwmax(i)  = 0.0_sp
             bvstrcwmax(i)  = 0.0_sp
             Ubot(i)        = 0.0_sp
             Vbot(i)        = 0.0_sp
             Ur(i)          = 0.0_sp
             Vr(i)          = 0.0_sp
             taucwmax(i)    = 0.0_sp
          end if
# endif	  
   end do


# if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustr     , bvstr      )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustrc    , bvstrc     )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustrw    , bvstrw     )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustrcwmax, bvstrcwmax )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Ubot      , Vbot       )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Ur        , Vr         )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,tauCWmax               ) 
# if defined (ORIG_SED)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bottom                 )
# elif defined (CSTMS_SED)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,sedbed%bottom                 )
# endif
# endif

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: mb_bbl"

  return

end subroutine mb_bbl




subroutine ssw_bbl
!=======================================================================
!                                                                      !
!  This routine compute bottom stresses for the case when the wave     !
!  solution in the wave boundary layer is based on a  2-layer eddy     !
!  viscosity that is linear increasing above Zo and constant above     !
!  Z1.                                                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Styles, R. and S.M. glenn,  2000: Modeling stratified wave and      !
!    current bottom boundary layers in the continental shelf, JGR,     !
!    105, 24119-24139.                                                 !
!                                                                      !
!=======================================================================
!
      logical :: ITERATE

      integer :: Iter, i, j, k

      real(sp), parameter :: eps = 1.0E-10

      real(sp) :: Kbh, Kbh2, Kdh
      real(sp) :: taucr, wsedr, tstar, coef_st
      real(sp) :: coef_b1, coef_b2, coef_b3, d0
      real(sp) :: dolam, dolam1, doeta1, doeta2, fdo_etaano
      real(sp) :: lamorb, lamanorb
      real(sp) :: m_ubr, m_wr, m_ucr, m_zr, m_phicw, m_kb
      real(sp) :: m_ustrc, m_ustrwm, m_ustrr, m_fwc, m_zoa
      real(sp) :: zo
      real(sp) :: Kb2, Kdelta, Ustr
      real(sp) :: anglec, anglew
      real(sp) :: cff, cff1, cff2, cff3, og, fac, fac1, fac2
      real(sp) :: sg_ab, sg_abokb, sg_a1, sg_b1, sg_chi, sg_c1, d50
      real(sp) :: sg_epsilon, ssw_eta, sg_fofa, sg_fofb, sg_fofc, sg_fwm
      real(sp) :: sg_kbs, ssw_lambda, sg_mu, sg_phicw, sg_ro, sg_row
      real(sp) :: sg_shdnrm, sg_shld, sg_shldcr, sg_scf, rhos, sg_star
      real(sp) :: sg_ub, sg_ubokur, sg_ubouc, sg_ubouwm, sg_ur
      real(sp) :: sg_ustarc, sg_ustarcw, sg_ustarwm, sg_znot, sg_znotp
      real(sp) :: sg_zr, sg_zrozn, sg_z1, sg_z1ozn, sg_z2, twopi, zd1, zd2
      real(sp) :: zoMIN, zoMAX
      real(sp) :: coef_fd,temp

      real(sp), parameter :: absolute_zoMIN = 5.0d-5  ! in Harris-Wiberg
!!    real(r8), parameter :: absolute_zoMIN = 5.0d-8  ! in Harris-Wiberg
      real(sp), parameter ::  Cd_fd = 0.5

      real(sp), parameter :: K1 = 0.6666666666     ! Coefficients for
      real(sp), parameter :: K2 = 0.3555555555     ! explicit
      real(sp), parameter :: K3 = 0.1608465608     ! wavenumber
      real(sp), parameter :: K4 = 0.0632098765     ! calculation
      real(sp), parameter :: K5 = 0.0217540484     ! (Dean and
      real(sp), parameter :: K6 = 0.0065407983     !  Dalrymple, 1991)

      real(sp), parameter :: coef_a1=0.095         ! Coefficients for
      real(sp), parameter :: coef_a2=0.442         ! ripple predictor
      real(sp), parameter :: coef_a3=2.280         ! (Wiberg-Harris)


      real(sp), parameter :: ar_gm82 = 27.7/30.0     ! Grant-Madsen (1982)
      real(sp), parameter :: ar_n92  = 0.267         ! Nielsen (1992)
      real(sp), parameter :: ar_r88  = 0.533         ! Raudkivi (1988)

      integer :: jj,cnt,cc
!#if defined GM82_RIPRUF
!      real(r8), parameter :: ar = 27.7/30.0  ! Grant-Madsen (1982)
!#elif defined N92_RIPRUF
!      real(r8), parameter :: ar = 0.267         ! Nielsen (1992)
!#elif defined R88_RIPRUF
!      real(r8), parameter :: ar = 0.533         ! Raudkivi (1988)
!#else
!      No ripple roughness coeff. chosen
!#endif

      real(sp), dimension(1:kbm1) :: Un,Vn

      real(sp), dimension(0:mt) :: Ab
      real(sp), dimension(0:mt) :: Fwave_bot
      real(sp), dimension(0:mt) :: Tauc
      real(sp), dimension(0:mt) :: Tauw
      real(sp), dimension(0:mt) :: Tauo
!      real(sp), dimension(0:mt) :: Taucwmax
      real(sp), dimension(0:mt) :: Ur_sg
      real(sp), dimension(0:mt) :: Vr_sg
      real(sp), dimension(0:mt) :: Ub
      real(sp), dimension(0:mt) :: Ucur
      real(sp), dimension(0:mt) :: Umag
      real(sp), dimension(0:mt) :: Vcur
      real(sp), dimension(0:mt) :: Zr
      real(sp), dimension(0:mt) :: phic
      real(sp), dimension(0:mt) :: phicw
      real(sp), dimension(0:mt) :: rheight
      real(sp), dimension(0:mt) :: rlength
      real(sp), dimension(0:mt) :: u100
      real(sp), dimension(0:mt) :: znot
      real(sp), dimension(0:mt) :: znotc
      real(sp), dimension(0:mt) :: zoN
      real(sp), dimension(0:mt) :: zoST
      real(sp), dimension(0:mt) :: zoBF
      real(sp), dimension(0:mt) :: zoDEF
      real(sp), dimension(0:mt) :: zoBIO

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: ssw_bbl"

!
!-----------------------------------------------------------------------
!  Set currents above bed.
!-----------------------------------------------------------------------
!
  sg_mp=CMPLX(SQRT(1.0/(2.0*sg_z1p)),SQRT(1.0/(2.0*sg_z1p)))
  twopi=2.0*pi

  do i=1,m
     Zr(i)=dt(i)*(zz(i,kbm1)-z(i,kb))
     Ur_sg(i) = 0.0
     Vr_sg(i) = 0.0
!     Un       = 0.0
!     Vn       = 0.0
     cnt = 0 
     do jj=1,ntve(i)
        cc = nbve(i,jj)
         !if(xc(cc) < vx(i))then !cell is downtream
           cnt =cnt + 1
           Ur_sg(i) = Ur_sg(i) + u(cc,kbm1)
           Vr_sg(i) = Vr_sg(i) + v(cc,kbm1)
!           Un(:)    = Un(:)    + u(cc,:)
!           Vn(:)    = Vn(:)    + v(cc,:)
         !endif
     end do
     Ur_sg(i) = Ur_sg(i) / cnt
     Vr_sg(i) = Vr_sg(i) / cnt
!     Un(:)    = Un(:)    / cnt
!     Vn(:)    = Vn(:)    / cnt
!     Un(:)    = sum(u(nbve(i,1:ntve(i)),1:kbm1),dim=1)/float(ntve(i))
!     Vn(:)    = sum(v(nbve(i,1:ntve(i)),1:kbm1),dim=1)/float(ntve(i))
!     Ur_sg(i) = Un(kbm1)
!     Vr_sg(i) = Vn(kbm1)

     if(SSW_LOGINT)then
         
!
!  If current height is less than z1ur, interpolate logarithmically
!  to z1ur. (This has not been updated wrt ssw...uses outdated zo defs)
!
          Un(:)    = sum(u(nbve(i,1:ntve(i)),1:kbm1),dim=1)/float(ntve(i))
          Vn(:)    = sum(v(nbve(i,1:ntve(i)),1:kbm1),dim=1)/float(ntve(i))
          IF (Zr(i).lt.sg_z1min) THEN
            DO k=kbm2,1,-1
              zd1=dt(i)*(zz(i,k+1)-z(i,kb))
              zd2=dt(i)*(zz(i,k  )-z(i,kb))
              IF ((zd1.lt.sg_z1min).and.(sg_z1min.lt.zd2)) THEN
                fac=1.0/LOG(zd2/zd1)
                fac1=fac*LOG(zd2/sg_z1min)
                fac2=fac*LOG(sg_z1min/zd1)
                Ur_sg(i)=fac1*Un(k+1)+fac2*Un(k)
                Vr_sg(i)=fac1*Vn(k+1)+fac2*Vn(k)
                Zr(i)=sg_z1min
              END IF
            END DO
          END IF
     endif
  end do


!
!-----------------------------------------------------------------------
!  Compute bottom stresses.
!-----------------------------------------------------------------------
!
  do i=1,m 
!
!  Set bed wave orbital velocity and excursion amplitude.  Use data
!  from wave models (SWAN) or use Dean and Dalrymple (1991) 6th-degree
!  polynomial to approximate wave number on shoaling water.

          Fwave_bot(i)=twopi/MAX(Pwave_bot(i),0.05)
          if(SSW_CALC_UB)then
            Kdh=h(i)*Fwave_bot(i)**2/grav
            Kbh2=Kdh*Kdh+                                     &
     &         Kdh/(1.0+Kdh*(K1+Kdh*(K2+Kdh*(K3+Kdh*(K4+      &
     &              Kdh*(K5+K6*Kdh))))))
            Kbh=SQRT(Kbh2)
            Ab(i)=0.5*HSC1(i)/SINH(Kbh)+eps
            Ub(i)=Fwave_bot(i)*Ab(i)+eps
          else
            Ub(i)=MAX(Ub_swan(i),0.0)+eps
            Ab(i)=Ub(i)/Fwave_bot(i)+eps
          endif
!
!  Compute bottom current magnitude at RHO-points.
!
          Ucur(i)=Ur_sg(i)
          Vcur(i)=Vr_sg(i)
          Umag(i)=SQRT(Ucur(i)*Ucur(i)+Vcur(i)*Vcur(i)+eps)
!
!  Compute angle between currents and waves (radians)
!
          IF (Ucur(i).eq.0.0) THEN
            phic(i)=0.5*pi*SIGN(1.0,Vcur(i))
          ELSE
            phic(i)=ATAN2(Vcur(i),Ucur(i))
          ENDIF
          phicw(i)=Dwave(i)-phic(i)
  end do


  do i=1,m
!
!  Load sediment properties, stresses, and roughness from previous time
!  step (stresses are in m2 s-2).
!
          d50=bottom(i,isd50)
          rhos=bottom(i,idens)/(1000.0*rho1(i,kbm1)+1000.0)
          wsedr=bottom(i,iwset)
          Taucr=bottom(i,itauc)/1025.
          Tauc(i)=SQRT(bustrc(i)**2+bvstrc(i)**2)
          Tauw(i)=SQRT(bustrw(i)**2+bvstrw(i)**2)
          Taucwmax(i)=SQRT( bustrcwmax(i)**2+bvstrcwmax(i)**2)
!
          rheight(i)=bottom(i,irhgt)
          rlength(i)=bottom(i,irlen)
          zoMAX=0.9*Zr(i)
          zoMIN=MAX(absolute_zoMIN,2.5*d50/30.0) 
!
!  Initialize arrays.
!
          zoN(i)=MIN(MAX(2.5*d50/30.0, zoMIN ),zoMAX)
          zoST(i)=0.0
          zoBF(i)=0.0
          zoBIO(i)=0.0
!!        zoDEF(i,j)=bottom(i,j,izdef)
          zoDEF(i)=BOTTOM_ROUGHNESS_LENGTHSCALE
        
         if(SSW_CALC_ZNOT)then
!
!  Calculate components of roughness and sum: zo = zoN + zoST + zoBF
!  Determine whether sediment is in motion. Use Shields criterion to
!  determine if sediment is mobile.
!
! print*,'Taucwmax=',Taucwmax(i),'Taucr=',Taucr
            tstar=Taucwmax(i)/(Taucr+eps)
            IF (tstar.lt.1.0) THEN                         ! no motion
              zoST(i)=0.0
              if(GM82_RIPRUF)then
                 zoBF(i)=ar_gm82*rheight(i)**2/rlength(i)
              elseif(N92_RIPRUF)then
                 zoBF(i)=ar_n92*rheight(i)**2/rlength(i)
              elseif(R88_RIPRUF)then
                 zoBF(i)=ar_r88*rheight(i)**2/rlength(i)
              else
                if(msr)write(ipt,*)'one of GM82_RIPRUF, N92_RIPRUF and&
                     & R88_RIPRUF needs be activated, stopping......  &
                     &'
                call pstop
              end if
            ELSE
!
!  Threshold of motion exceeded - calculate new zoST and zoBF
!  Calculate saltation roughness according to Wiberg & Rubin (1989)
!  (Eqn. 11 in Harris & Wiberg, 2001)
!  (d50 is in m, but this formula needs cm)  
!
               coef_st=0.0204*LOG(100.0*d50)**2+                    &
     &                 0.0220*LOG(100.0*d50)+0.0709
               zoST(i)=0.056*d50*0.68*tstar/                      &
     &                 (1.0+coef_st*tstar)
               IF (zoST(i).lt.0.0) THEN
                 IF (MSR) THEN
                   PRINT *, ' Warning: zoST<0  tstar, d50, coef_st:'
                   PRINT *, tstar,d50,coef_st
                 END IF
               END IF
!
!  Calculate ripple height and wavelength.
!  Use Malarkey & Davies (2003) explict version of Wiberg & Harris.
!
               coef_b1=1.0/coef_a1
               coef_b2=0.5*(1.0 + coef_a2)*coef_b1
               coef_b3=coef_b2**2-coef_a3*coef_b1
               d0=2.0*Ab(i)
!print*,d50,d0/d50,13000
               IF ((d0/d50).gt.13000.0) THEN              ! sheet flow
                 rheight(i)=0.0
                 rlength(i)=535.0*d50        ! does not matter since
               ELSE                               ! rheight=0 
                 dolam1=d0/(535.0*d50)
                 doeta1=EXP(coef_b2-SQRT(coef_b3-coef_b1*LOG(dolam1)))
                 lamorb=0.62*d0
                 lamanorb=535.0*d50
                 IF (doeta1.lt.20.0) THEN
                   dolam=1.0/0.62
                 ELSE IF (doeta1.gt.100.0) THEN
                   dolam=dolam1
                 ELSE
                   fdo_etaano=-LOG(lamorb/lamanorb)*                      &
     &                         LOG(0.01*doeta1)/LOG(5.0)
                   dolam=dolam1*EXP(-fdo_etaano)
                 END IF
                 doeta2=EXP(coef_b2-SQRT(coef_b3-coef_b1*LOG(dolam)))
                 rheight(i)=d0/doeta2
                 rlength(i)=d0/dolam
               END IF
!print*,rlength(i)
!
!  Value of ar can range from 0.3 to 3 (Soulsby, 1997, p. 124)
!
              if(GM82_RIPRUF)then
                 zoBF(i)=ar_gm82*rheight(i)**2/rlength(i)
              elseif(N92_RIPRUF)then
                 zoBF(i)=ar_n92*rheight(i)**2/rlength(i)
              elseif(R88_RIPRUF)then
                 zoBF(i)=ar_r88*rheight(i)**2/rlength(i)
              end if
!               zoBF(i)=ar*rheight(i)**2/rlength(i)

            END IF
            zo=zoN(i)
            if(SSW_ZOBL)then
              zo=zo+zoST(i)
            endif
            if(SSW_ZORIP)then
              zo=zo+zoBF(i)
            endif
            if(SSW_ZOBIO)then
              zo=zo+zoBIO(i)
            endif
    else
            IF (zoDEF(i).lt.absolute_zoMIN) THEN
               temp = zoDEF(i)
               zoDEF(i)=absolute_zoMIN
               IF (MSR) THEN
                 PRINT *, ' Warning: default zo (',temp,') <',absolute_zoMIN,' mm, replaced with:',&
     &                 zoDEF(i)
                END IF
            END IF
            zo=zoDEF(i)
     endif
!
!  Compute stresses.
!
!  Default stress calcs for pure currents
!
          zo=MIN(MAX(zo,zoMIN),zoMAX)
!
          cff1=vonKar/LOG(Zr(i)/zo)
          cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
          Tauc(i)=cff2*Umag(i)*Umag(i)
          Tauo(i)=Tauc(i)
          Tauw(i)=0.0
          Taucwmax(i)=Tauc(i)
          znot(i)=zo
          znotc(i)=zo
!
          IF ((Umag(i).le.eps).and.(Ub(i).gt.eps)) THEN
!
!  Pure waves - use wave friction factor approach from Madsen
!  (1994, eqns 32-33).
!
            sg_abokb=Ab(i)/(30.0*zo)
            sg_fwm=0.3
            IF ((sg_abokb.gt.0.2).and.(sg_abokb.le.100.0)) THEN
              sg_fwm=EXP(-8.82+7.02*sg_abokb**(-0.078))
            ELSE IF (sg_abokb.gt.100.0)THEN
              sg_fwm=EXP(-7.30+5.61*sg_abokb**(-0.109))
            END IF
            Tauc(i)= 0.0
            Tauw(i)= 0.5*sg_fwm*Ub(i)*Ub(i)
            Taucwmax(i)=Tauw(i)
            znot(i)=zo
            znotc(i)=zo
          ELSE IF ((Umag(i).gt.0.0).and.(Ub(i).gt.eps).and.      &
     &             ((Zr(i)/zo).le.1.0)) THEN
!
!  Waves and currents, but zr <= zo.
!
            IF (MSR) THEN
              PRINT *,' Warning: w-c calcs ignored because zr <= zo'
            END IF
          ELSE IF ((Umag(i).gt.0.0).and.(Ub(i).gt.eps).and.      &
     &             ((Zr(i)/zo).gt.1.0)) THEN
!
!  Waves and currents, zr > zo.
!
            if(SGWC)then
              sg_zrozn=Zr(i)/zo
              sg_ubokur=Ub(i)/(sg_kappa*Umag(i))
              sg_row=Ab(i)/zo
              sg_a1=1.0d-6
              sg_phicw=phicw(i)
              CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur,     &
     &                       sg_a1, sg_mu, sg_epsilon, sg_ro, sg_fofa)
              sg_abokb=Ab(i)/(30.0*zo)
              IF (sg_abokb.le.100.0) THEN
                sg_fwm=EXP(-8.82+7.02*sg_abokb**(-0.078))
              ELSE
                sg_fwm=EXP(-7.30+5.61*sg_abokb**(-0.109))
              END IF
              sg_ubouwm=SQRT(2.0/sg_fwm)
!
!  Determine the maximum ratio of wave over combined shear stresses,
!  sg_ubouwm (ub/ustarwm).
!
              CALL sg_purewave (sg_row, sg_ubouwm, sg_znotp, sg_ro)
!
!  Set initial guess of the ratio of wave over shear stress, sg_c1
!  (ub/ustarc).
!
              sg_b1=sg_ubouwm
              sg_fofb=-sg_fofa
              sg_c1=0.5*(sg_a1+sg_b1)
              CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur,     &
     &                       sg_c1, sg_mu, sg_epsilon, sg_ro, sg_fofc)
!
!  Solve PDE via bi-section method.
!
              ITERATE=.true.
              DO Iter=1,sg_n
                IF (ITERATE) THEN
                  IF ((sg_fofb*sg_fofc).lt.0.0) THEN
                    sg_a1=sg_c1
                  ELSE
                    sg_b1=sg_c1
                  END IF
                  sg_c1=0.5*(sg_a1+sg_b1)
                  CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur, &
     &                           sg_c1, sg_mu, sg_epsilon, sg_ro,       &
     &                           sg_fofc)
                  ITERATE=(sg_b1-sg_c1) .ge. sg_tol
                  IF (ITERATE) Iconv(i)=Iter
                END IF
              END DO
              sg_ubouc=sg_c1
!
!  Compute bottom shear stress magnitude (m/s).
!
              sg_ustarcw=Ub(i)/sg_ubouc
              sg_ustarwm=sg_mu*sg_ustarcw
!!            sg_ustarc=MIN(sg_ustarcdef,sg_epsilon*sg_ustarcw)  !original
              sg_ustarc=MAX(SQRT(Tauc(i)),sg_epsilon*sg_ustarcw)
              Tauc(i)=sg_ustarc*sg_ustarc
              Tauw(i)=sg_ustarwm*sg_ustarwm
              Taucwmax(i)=SQRT((Tauc(i)+                              &
     &                            Tauw(i)*COS(phicw(i)))**2+          &
     &                           (Tauw(i)*SIN(phicw(i)))**2)
!
!  Compute apparent hydraulic roughness (m).
!
              IF (sg_epsilon.gt.0.0) THEN
                sg_z1=sg_alpha*sg_kappa*Ab(i)/sg_ubouc
                sg_z2=sg_z1/sg_epsilon
                sg_z1ozn=sg_z1/zo
                znotc(i)=sg_z2*                                         &
     &                     EXP(-(1.0-sg_epsilon+                       &
     &                         sg_epsilon*LOG(sg_z1ozn)))
!
!  Compute mean (m/s) current at 100 cm above the bottom.
!
                IF (sg_z100.gt.sg_z2) THEN
                  u100(i)=sg_ustarc*                                    &
     &                    (LOG(sg_z100/sg_z2)+1.0-sg_epsilon+        &
     &                    sg_epsilon*LOG(sg_z1ozn))/sg_kappa
                ELSE IF ((sg_z100.le.sg_z2).and.(Zr(i).gt.sg_z1)) THEN
                  u100(i)=sg_ustarc*sg_epsilon*                         &
     &                    (sg_z100/sg_z1-1.0+LOG(sg_z1ozn))/sg_kappa
                ELSE
                  u100(i)=sg_ustarc*sg_epsilon*                         &
     &                    LOG(sg_z100/zo)/sg_kappa
                END IF
              END IF
            elseif(M94WC)then
              m_ubr=Ub(i)
              m_wr=Fwave_bot(i)
              m_ucr=Umag(i)
              m_zr=Zr(i)
              m_phicw=phicw(i)
              m_kb=30.0*zo
              CALL madsen94 (m_ubr, m_wr, m_ucr,                          &
     &                      m_zr, m_phicw, m_kb,                         &
     &                      m_ustrc, m_ustrwm, m_ustrr, m_fwc, m_zoa)
              Tauc(i)=m_ustrc*m_ustrc
              Tauw(i)=m_ustrwm*m_ustrwm
              Taucwmax(i)=m_ustrr*m_ustrr
              znotc(i)=min( m_zoa, zoMAX )
              u100(i)=(m_ustrc/vonKar)*LOG(1.0/m_zoa)
            endif
            if(SSW_FORM_DRAG_COR)then
              IF (rheight(i).gt.(zoN(i)+zoST(i))) THEN
                coef_fd=0.5*Cd_fd*(rheight(i)/rlength(i))*         &
     &                (1.0/(vonKar*vonKar))*                         &
     &                (LOG(rheight(i)/                                &
     &                 (zoN(i)+zoST(i))-1.0))**2
                 Taucwmax(i)=Taucwmax(i)/(1.0+coef_fd)
                 Taucwmax(i)=Taucwmax(i)*(1.0+8.0*              &
     &                       rheight(i)/rlength(i))
              END IF
            endif
          END IF
   end do

!-----------------------------------------------------------------------
!  Compute kinematic bottom stress components due current and wind-
!  induced waves.
!-----------------------------------------------------------------------
!
  do i=1,m
          anglec=Ur_sg(i)/Umag(i)
          bustr(i)=Tauc(i)*anglec
          !bustro(i)=Tauo(i)*anglec
          anglec=Vr_sg(i)/Umag(i)
          bvstr(i)=Tauc(i)*anglec
          !bvstro(i)=Tauo(i)*anglec
  end do

  do i=1,m
          anglec=Ucur(i)/Umag(i)
          anglew=COS(Dwave(i))
          bustrc(i)=Tauc(i)*anglec
          bustrw(i)=Tauw(i)*anglew
          bustrcwmax(i)=Taucwmax(i)*anglew
          Ubot(i)=Ub(i)*anglew
          Ur(i)=Ucur(i)
!
          anglec=Vcur(i)/Umag(i)
          anglew=SIN(Dwave(i))
          bvstrc(i)=Tauc(i)*anglec
          bvstrw(i)=Tauw(i)*anglew
          bvstrcwmax(i)=Taucwmax(i)*anglew
          Vbot(i)=Ub(i)*anglew
          Vr(i)=Vcur(i)
!
          bottom(i,ibwav)=Ab(i)
          bottom(i,irhgt)=rheight(i)
          bottom(i,irlen)=rlength(i)
          bottom(i,izdef)=zoDEF(i)
          bottom(i,izapp)=znotc(i)
          bottom(i,izNik)=zoN(i)
          bottom(i,izbio)=zoBIO(i)
          bottom(i,izbfm)=zoBF(i)
          bottom(i,izbld)=zoST(i)
          bottom(i,izwbl)=znot(i)
	  
# if defined (WET_DRY)	  
          if(iswetn(i)/=1)then
             bustr(i)       = 0.0_sp
             bvstr(i)       = 0.0_sp
             bustrc(i)      = 0.0_sp
             bvstrc(i)      = 0.0_sp
             bustrw(i)      = 0.0_sp
             bvstrw(i)      = 0.0_sp
             bustrcwmax(i)  = 0.0_sp
             bvstrcwmax(i)  = 0.0_sp
             Ubot(i)        = 0.0_sp
             Vbot(i)        = 0.0_sp
             Ur(i)          = 0.0_sp
             Vr(i)          = 0.0_sp
             taucwmax(i)    = 0.0_sp
          end if
# endif	  
  end do

# if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustr     , bvstr      )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustrc    , bvstrc     )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustrw    , bvstrw     )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bustrcwmax, bvstrcwmax )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Ubot      , Vbot       )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Ur        , Vr         )
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,taucwmax               ) 
# if defined (ORIG_SED)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bottom                 )
# elif defined (CSTMS_SED)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,sedbed%bottom                 )
# endif
# endif

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: ssw_bbl"

  return

end subroutine ssw_bbl




 SUBROUTINE sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur,     &
         &              sg_ubouc, sg_mu, sg_epsilon, sg_ro,        &
         &              sg_fofx)
!
!=======================================================================
!                                                                      !
!  This routine computes bottom stresses via bottom boundary layer     !
!  formulation of Styles and Glenn (1999).                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     sg_row      Ratio of wave excursion amplitude over roughness.    !
!     sg_zrozn    Ratio of height of current over roughness.           !
!     sg_phiwc    Angle between wave and currents (radians).           !
!     sg_ubokur   Ratio of wave over current velocity:                 !
!                   ub/(vonKar*ur)                                     !
!     sg_ubouc    Ratio of bed wave orbital over bottom shear stress   !
!                   (ub/ustarc), first guess.                          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     sg_ubouc    Ratio of bed wave orbital over bottom shear stress   !
!                   (ub/ustarc), iterated value.                       !
!     sg_mu       Ratio between wave and current bottom shear          !
!                   stresses (ustarwm/ustarc).                         !
!     sg_epsilon  Ratio between combined (wave and current) and        !
!                   current bottom shear stresses (ustarc/ustarcw).    !
!     sg_ro       Internal friction Rossby number:                     !
!                   ustarc/(omega*znot)                                !
!     sg_fofx     Root of PDE used for convergence.                    !
!                                                                      !
!=======================================================================
!

!
!  Imported variable declarations.
!
      real(sp), intent(in) :: sg_row, sg_zrozn, sg_phicw, sg_ubokur

      real(sp), intent(inout) :: sg_ubouc

      real(sp), intent(out) :: sg_mu, sg_epsilon, sg_ro, sg_fofx
!
!  Local variable declarations.
!
      logical :: ITERATE

      integer :: Iter

      real(sp) :: cff, sg_bei, sg_beip, sg_ber, sg_berp, sg_cosphi
      real(sp) :: sg_eps2, sg_kei, sg_keip, sg_ker, sg_kerp, sg_mu2
      real(sp) :: sg_phi, sg_ror, sg_x, sg_z2p, sg_znotp, sg_zroz1
      real(sp) :: sg_zroz2, sg_z1ozn, sg_z2ozn

      complex  :: sg_argi, sg_bnot, sg_bnotp, sg_b1, sg_b1p
      complex  :: sg_gammai, sg_knot, sg_knotp, sg_k1, sg_k1p
      complex  :: sg_ll, sg_nn

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: sg_bstress"
!
!-----------------------------------------------------------------------
!  Compute bottom stresses.
!-----------------------------------------------------------------------
!
!  Compute nondimensional bottom wave shear, phi.  Iterate to make
!  sure that there is an upper limit in "ubouc".  It usually requires
!  only one pass.
!
      ITERATE=.TRUE.
      DO Iter=1,sg_n
        IF (ITERATE) THEN
          sg_ro=sg_row/sg_ubouc
          sg_znotp=1.0/(sg_kappa*sg_ro)
          IF ((sg_z1p/sg_znotp).gt.1.0) THEN
            sg_x=2.0*SQRT(sg_znotp)
            IF (sg_x.le.8.0) THEN
              CALL sg_kelvin8m (sg_x, sg_ber, sg_bei, sg_ker, sg_kei,   &
     &                          sg_berp, sg_beip, sg_kerp, sg_keip)
            ELSE
              CALL sg_kelvin8p (sg_x, sg_ker, sg_kei, sg_ber, sg_bei,   &
     &                          sg_kerp, sg_keip, sg_berp, sg_beip)
            END IF
            cff=1.0/SQRT(sg_znotp)
            sg_bnot =CMPLX(sg_ber,sg_bei)
            sg_knot =CMPLX(sg_ker,sg_kei)
            sg_bnotp=CMPLX(sg_berp,sg_beip)*cff
            sg_knotp=CMPLX(sg_kerp,sg_keip)*cff
!
            sg_x=2.0*SQRT(sg_z1p)
            IF (sg_x.le.8.0) THEN
              CALL sg_kelvin8m (sg_x, sg_ber, sg_bei, sg_ker, sg_kei,   &
     &                          sg_berp, sg_beip, sg_kerp, sg_keip)
            ELSE
              CALL sg_kelvin8p (sg_x, sg_ker, sg_kei, sg_ber, sg_bei,   &
     &                          sg_kerp, sg_keip, sg_berp, sg_beip)
            END IF
            cff=1.0/SQRT(sg_z1p)
            sg_b1 =CMPLX(sg_ber,sg_bei)
            sg_k1 =CMPLX(sg_ker,sg_kei)
            sg_b1p=CMPLX(sg_berp,sg_beip)*cff
            sg_k1p=CMPLX(sg_kerp,sg_keip)*cff
!
            sg_ll=sg_mp*sg_b1+sg_b1p
            sg_nn=sg_mp*sg_k1+sg_k1p
            sg_argi=sg_bnotp*sg_nn/(sg_bnot*sg_nn-sg_knot*sg_ll)+       &
     &              sg_knotp*sg_ll/(sg_knot*sg_ll-sg_bnot*sg_nn)
            sg_gammai=-sg_kappa*sg_znotp*sg_argi
            sg_phi=CABS(sg_gammai)
          ELSE
            sg_gammai=-sg_kappa*sg_z1p*sg_mp
            sg_phi=CABS(sg_gammai)
          END IF
!
          IF (sg_ubouc.gt.(1.0/sg_phi)) THEN
            sg_ubouc=1.0/sg_phi
          ELSE
            ITERATE=.FALSE.
          END IF
        END IF
      END DO
!
!  Compute ratio of wave over current bottom shear stresses.
!
      sg_mu=SQRT(sg_ubouc*sg_phi)
!
!  Compute ratio of current over combined bottom shear stresses.
!
      IF (sg_mu.eq.1.0) THEN
        sg_epsilon=0.0
      ELSE
        sg_mu2=sg_mu*sg_mu
        sg_cosphi=ABS(COS(sg_phicw))
        sg_eps2=-sg_mu2*sg_cosphi+                                      &
     &          SQRT(1.0+sg_mu2*sg_mu2*(sg_cosphi*sg_cosphi-1.0))
        sg_epsilon=SQRT(sg_eps2)
      END IF
!
!  Determine root of PDE used for convergence.
!
      IF (sg_epsilon.ne.0.0) THEN
        sg_z2p=sg_z1p/sg_epsilon
        sg_ror=sg_ro/sg_zrozn
        sg_zroz1=1.0/(sg_alpha*sg_kappa*sg_ror)
        sg_zroz2=sg_epsilon*sg_zroz1
        sg_z1ozn=sg_alpha*sg_kappa*sg_ro
        sg_z2ozn=sg_z1ozn/sg_epsilon
!
        IF ((sg_zroz2.gt.1.0).and.(sg_z1ozn.gt.1.0)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*                       &
     &                      (LOG(sg_zroz2)+1.0-sg_epsilon+           &
     &                       sg_epsilon*LOG(sg_z1ozn))
        ELSE IF ((sg_zroz2.le.1.0).and.(sg_zroz1.gt.1.0).and.     &
     &          (sg_z1ozn.gt.1.0)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*sg_epsilon*            &
     &                      (sg_zroz1-1.0+LOG(sg_z1ozn))
        ELSE IF ((sg_zroz1.le.1.0).and.(sg_z1ozn.gt.1.0)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*sg_epsilon*            &
     &                      LOG(sg_zrozn)
        ELSE IF ((sg_zroz2.gt.1.0).and.(sg_z1ozn.le.1.0).and.     &
     &          (sg_z2ozn.gt.1.0)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*                       &
     &                      (LOG(sg_zroz2)+1.0-1.0/sg_z2ozn)
        ELSE IF ((sg_zroz2.le.1.0).and.(sg_zroz1.gt.1.0).and.     &
     &          (sg_z1ozn.le.1.0).and.(sg_z2ozn.gt.1.0)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*sg_epsilon*            &
     &                      (sg_zroz1-1.0/sg_z1ozn)
        ELSE IF ((sg_zroz2.gt.1.0).and.(sg_z2ozn.le.1.0)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*LOG(sg_zrozn)
        END IF
      END IF

     if(dbg_set(dbg_sbr)) write(ipt,*) "End: sg_bstress"

      RETURN
      END SUBROUTINE sg_bstress

      SUBROUTINE  sg_purewave (sg_row, sg_ubouwm, sg_znotp, sg_ro)
!
!=======================================================================
!                                                                      !
!  This routine determines the maximum ratio of waves over combined    !
!  bottom shear stress.                                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     sg_row      Ratio of wave excursion amplitude over roughness.    !
!     sg_ubouwm   Maximum ratio of waves over combined bottom shear    !
!                   stress.                                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     sg_ubouwm   Maximum ratio of waves over combined bottom shear    !
!                   stress.                                            !
!     sg_znotp    Ratio of hydraulic roughness over scaled height      !
!                   of bottom boundary layer.                          !
!     sg_ro       Internal friction Rossby number:                     !
!                   ustarc/(omega*znot)                                !
!                                                                      !
!=======================================================================
!

!
!  Imported variable declarations.
!
      real(sp), intent(in) :: sg_row

      real(sp), intent(inout) :: sg_ubouwm

      real(sp), intent(out) :: sg_znotp, sg_ro
!
!  Local variable declarations.
!
      integer :: Iter

      real(sp) :: cff, sg_bei, sg_beip, sg_ber, sg_berp, sg_kei
      real(sp) :: sg_keip, sg_ker, sg_kerp, sg_phi, sg_ubouwmn, sg_x

      complex  :: sg_argi, sg_bnot, sg_bnotp, sg_b1, sg_b1p
      complex  :: sg_gammai, sg_knot, sg_knotp, sg_k1, sg_k1p
      complex  :: sg_ll, sg_nn

     if(dbg_set(dbg_sbr)) write(ipt,*) "Start: sg_purewave"
!
!-----------------------------------------------------------------------
!  Compute wind-induced wave stress.
!-----------------------------------------------------------------------
!
      DO Iter=1,sg_n
        sg_ro=sg_row/sg_ubouwm
        sg_znotp=1.0/(sg_kappa*sg_ro)
        IF (sg_z1p/sg_znotp.gt.1.0) THEN
          sg_x=2.0*SQRT(sg_znotp)
          IF (sg_x.le.8.0) THEN
            CALL sg_kelvin8m (sg_x, sg_ber, sg_bei, sg_ker, sg_kei,     &
     &                        sg_berp, sg_beip, sg_kerp, sg_keip)
          ELSE
            CALL sg_kelvin8p (sg_x, sg_ker, sg_kei, sg_ber, sg_bei,     &
     &                        sg_kerp, sg_keip, sg_berp, sg_beip)
          END IF
          cff=1.0/SQRT(sg_znotp)
          sg_bnot =CMPLX(sg_ber,sg_bei)
          sg_knot =CMPLX(sg_ker,sg_kei)
          sg_bnotp=CMPLX(sg_berp,sg_beip)*cff
          sg_knotp=CMPLX(sg_kerp,sg_keip)*cff
!
          sg_x=2.0*SQRT(sg_z1p)
          IF (sg_x.le.8.0) THEN
            CALL sg_kelvin8m (sg_x, sg_ber, sg_bei, sg_ker, sg_kei,     &
     &                        sg_berp, sg_beip, sg_kerp, sg_keip)
          ELSE
            CALL sg_kelvin8p (sg_x, sg_ker, sg_kei, sg_ber, sg_bei,     &
     &                        sg_kerp, sg_keip, sg_berp, sg_beip)
          END IF
          cff=1.0/SQRT(sg_z1p)
          sg_b1 =CMPLX(sg_ber,sg_bei)
          sg_k1 =CMPLX(sg_ker,sg_kei)
          sg_b1p=CMPLX(sg_berp,sg_beip)*cff
          sg_k1p=CMPLX(sg_kerp,sg_keip)*cff
!
          sg_ll=sg_mp*sg_b1+sg_b1p
          sg_nn=sg_mp*sg_k1+sg_k1p
          sg_argi=sg_bnotp*sg_nn/(sg_bnot*sg_nn-sg_knot*sg_ll)+         &
     &            sg_knotp*sg_ll/(sg_knot*sg_ll-sg_bnot*sg_nn)
          sg_gammai=-sg_kappa*sg_znotp*sg_argi
          sg_phi=CABS(sg_gammai)
        ELSE
          sg_gammai=-sg_kappa*sg_z1p*sg_mp
          sg_phi=CABS(sg_gammai)
        END IF
!
        sg_ubouwmn=1.0/sg_phi
        IF (abs((sg_ubouwmn-sg_ubouwm)/sg_ubouwmn).le.sg_tol) THEN
          sg_ubouwm=sg_ubouwmn
          RETURN
        ELSE
          sg_ubouwm=sg_ubouwmn
        END IF
      END DO

      if(dbg_set(dbg_sbr)) write(ipt,*) "End: sg_purewave"

      RETURN
      END SUBROUTINE  sg_purewave

      SUBROUTINE sg_kelvin8m (x, ber, bei, ker, kei, berp, beip,        &
     &                        kerp, keip)
!
!=======================================================================
!                                                                      !
! This rotuine computes the Kelvin functions for arguments less        !
! than eight.                                                          !
!                                                                      !
!=======================================================================
!
!
!  Imported variable declarations.
!
      real(sp), intent(in) :: x
      real(sp), intent(out) :: ber, bei, ker, kei
      real(sp), intent(out) :: berp, beip, kerp, keip
!
!  Local variable declarations.
!
      integer :: i

      real(sp) :: cff, xhalf

      real(sp), dimension(28) :: xp
!
!-----------------------------------------------------------------------
!  Compute Kelvin functions.
!-----------------------------------------------------------------------
!
      cff=0.125*x
      xp(1)=cff
      DO i=2,28
        xp(i)=xp(i-1)*cff
      END DO
      xhalf=0.5*x
!
      ber=1.0-                                                       &
     &    64.0*xp(4)+113.77777774*xp(8)-                          &
     &    32.36345652*xp(12)+2.64191397*xp(16)-                   &
     &    0.08349609*xp(20)+0.00122552*xp(24)-                    &
     &    0.00000901*xp(28)
      bei=16.0*xp(2)-113.77777774*xp(6)+                          &
     &    72.81777742*xp(10)-10.56765779*xp(14)+                     &
     &    0.52185615*xp(18)-0.01103667*xp(22)+                    &
     &    0.00011346*xp(26)
!
      ker=-ber*LOG(xhalf)+0.25*pi*bei-                               &
     &    0.57721566-59.05819744*xp(4)+                              &
     &    171.36272133*xp(8)-60.60977451*xp(12)+                  &
     &    5.65539121*xp(16)-0.19636347*xp(20)+                    &
     &    0.00309699*xp(24)-0.00002458*xp(28)
      kei=-bei*LOG(xhalf)-0.25*pi*ber+                               &
     &    6.76454936*xp(2)-142.91827687*xp(6)+                    &
     &    124.23569650*xp(10)-21.30060904*xp(14)+                 &
     &    1.17509064*xp(18)-0.02695875*xp(22)+                    &
     &    0.00029532*xp(26)
!
      berp=x*(-4.0*xp(2)+14.22222222*xp(6)-                       &
     &        6.06814810*xp(10)+0.66047849*xp(14)-                &
     &        0.02609253*xp(18)+0.00045957*xp(22)-                &
     &        0.00000394*xp(26))
      beip=x*(0.5-10.66666666*xp(4)+11.37777772*xp(8)-         &
     &        2.31167514*xp(12)+0.14677204*xp(16)-                &
     &        0.00379386*xp(20)+0.00004609*xp(24))
!
      kerp=-berp*LOG(xhalf)-ber/x+0.25*pi*beip+                         &
     &     x*(-3.69113734*xp(2)+21.42034017*xp(6)-                &
     &        11.36433272*xp(10)+1.41384780*xp(14)-               &
     &        0.06136358*xp(18)+0.00116137*xp(22)-                &
     &        0.00001075*xp(26))
      keip=-beip*LOG(xhalf)-bei/x-0.25*pi*berp+                      &
     &     x*(0.21139217-13.39858846*xp(4)+                       &
     &        19.41182758*xp(8)-4.65950823*xp(12)+                &
     &        0.33049424*xp(16)-0.00926707*xp(20)+                &
     &        0.00011997*xp(24))
      RETURN
      END SUBROUTINE sg_kelvin8m

      SUBROUTINE sg_kelvin8p (x, ker, kei, ber, bei, kerp, keip,        &
     &                        berp, beip)
!
!=======================================================================
!                                                                      !
! This rotuine computes the Kelvin functions for arguments greater     !
! than eight.                                                          !
!                                                                      !
!=======================================================================
!

!
!  Imported variable declarations.
!
      real(sp), intent(in) :: x
      real(sp), intent(out) :: ker, kei, ber, bei
      real(sp), intent(out) :: kerp, keip, berp, beip
!
!  Local variable declarations.
!
      integer :: i

      real(sp) :: cff, xhalf

      real(sp), dimension(6) :: xm, xp

      complex :: argm, argp, fofx, gofx, phim, phip, thetam, thetap
!
!-----------------------------------------------------------------------
!  Compute Kelvin functions.
!-----------------------------------------------------------------------
!
      cff=8.0/x
      xp(1)=cff
      xm(1)=-cff
      DO i=2,6
        xp(i)=xp(i-1)*cff
        xm(i)=-xm(i-1)*cff
      END DO
!
      thetap=CMPLX(0.0,-0.3926991)+                               &
     &       CMPLX(0.0110486,-0.0110485)*xp(1)+                   &
     &       CMPLX(0.0,-0.0009765)*xp(2)+                         &
     &       CMPLX(-0.0000906,-0.0000901)*xp(3)+                  &
     &       CMPLX(-0.0000252,0.0)*xp(4)+                         &
     &       CMPLX(-0.0000034,0.0000051)*xp(5)+                   &
     &       CMPLX(0.0000006,0.0000019)*xp(6)
      thetam=CMPLX(0.0,-0.3926991)+                               &
     &       CMPLX(0.0110486,-0.0110485)*xm(1)+                   &
     &       CMPLX(0.0,-0.0009765)*xm(2)+                         &
     &       CMPLX(-0.0000906,-0.0000901)*xm(3)+                  &
     &       CMPLX(-0.0000252,0.0)*xm(4)+                         &
     &       CMPLX(-0.0000034,0.0000051)*xm(5)+                   &
     &       CMPLX(0.0000006,0.0000019)*xm(6)
!
      phip=CMPLX(0.7071068,0.7071068)+                            &
     &     CMPLX(-0.0625001,-0.0000001)*xp(1)+                    &
     &     CMPLX(-0.0013813,0.0013811)*xp(2)+                     &
     &     CMPLX(0.0000005,0.0002452)*xp(3)+                      &
     &     CMPLX(0.0000346,0.0000338)*xp(4)+                      &
     &     CMPLX(0.0000117,-0.0000024)*xp(5)+                     &
     &     CMPLX(0.0000016,-0.0000032)*xp(6)
      phim=CMPLX(0.7071068,0.7071068)+                            &
     &     CMPLX(-0.0625001,-0.0000001)*xm(1)+                    &
     &     CMPLX(-0.0013813,0.0013811)*xm(2)+                     &
     &     CMPLX(0.0000005,0.0002452)*xm(3)+                      &
     &     CMPLX(0.0000346,0.0000338)*xm(4)+                      &
     &     CMPLX(0.0000117,-0.0000024)*xm(5)+                     &
     &     CMPLX(0.0000016,-0.0000032)*xm(6)
!
      cff=x/SQRT(2.0)
      argm=-cff*CMPLX(1.0,1.0)+thetam
      fofx=SQRT(pi/(2.0*x))*CEXP(argm)
      ker=REAL(fofx)
      kei=AIMAG(fofx)
!
      argp=cff*CMPLX(1.0,1.0)+thetap
      gofx=1.0/SQRT(2.0*pi*x)*CEXP(argp)
      ber=REAL(gofx)-kei/pi
      bei=AIMAG(gofx)+ker/pi
!
      kerp=REAL(-fofx*phim)
      keip=AIMAG(-fofx*phim)
!
      berp=REAL(gofx*phip)-keip/pi
      beip=AIMAG(gofx*phip)+kerp/pi
      RETURN
      END SUBROUTINE sg_kelvin8p


      SUBROUTINE madsen94 (ubr, wr, ucr, zr, phiwc, kN,                 &
     &                     ustrc, ustrwm, ustrr, fwc, zoa)
!
!=======================================================================
!                                                                      !
!  Grant-Madsen model from Madsen (1994).                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ubr     Rep. wave-orbital velocity amplitude outside WBL (m/s).  !
!     wr      Rep. angular wave frequency,  2* pi/T (rad/s).           !
!     ucr     Current velocity at height zr (m/s).                     !
!     zr      Reference height for current velocity (m).               !
!     phiwc   Angle between currents and waves at zr (radians).        !
!     kN      Bottom roughness height, like Nikuradse k, (m).          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ustrc   Current friction velocity, u*c (m/s).                    !
!     ustrwm  Wave maximum friction velocity, u*wm (m/s).              !
!     ustrr   Wave-current combined friction velocity, u*r (m/s).      !
!     fwc     Wave friction factor (nondimensional).                   !
!     zoa     Apparent bottom roughness (m).                           !
!                                                                      !
!=======================================================================
!

!
!  Imported variable declarations.
!
      real(sp), intent(in) ::  ubr, wr, ucr, zr, phiwc, kN
      real(sp), intent(out) ::  ustrc, ustrwm, ustrr, fwc, zoa
!
!  Local variable declarations.
!
      integer, parameter :: MAXIT = 20
      integer :: i, nit
      integer :: iverbose = 1

      real(sp) :: bigsqr, cosphiwc, cukw, diff, lndw, lnln, lnzr
      real(sp) :: phicwc, zo
      real(sp) :: dval = 99.99

      real(sp), dimension(MAXIT) :: Cmu
      real(sp), dimension(MAXIT) :: dwc
      real(sp), dimension(MAXIT) :: fwci
      real(sp), dimension(MAXIT) :: rmu
      real(sp), dimension(MAXIT) :: ustrci
      real(sp), dimension(MAXIT) :: ustrr2
      real(sp), dimension(MAXIT) :: ustrwm2
!
!-----------------------------------------------------------------------
!  Compute bottom friction velocities and roughness.
!-----------------------------------------------------------------------
!
!  Set special default values.
!
      ustrc=dval
      ustrwm=dval
      ustrr=dval
      fwc=0.4
      zoa=kN/30.0
      phicwc=phiwc

      zo = kN/30.0
  
      IF (ubr.le.0.01) THEN
        IF (ucr.le. 0.01) THEN          ! no waves or currents
          ustrc=0.0
          ustrwm=0.0
          ustrr=0.0
          RETURN
        END IF
        ustrc=ucr*vonKar/LOG(zr/zo)        ! no waves
        ustrwm=0.0
        ustrr=ustrc
        RETURN
      END IF
!
!  Iterate to compute friction velocities, roughness, and wave friction
!  factor.  Notice that the computation of the wave friction factor
!  has been inlined for efficiency.
!
      cosphiwc=ABS(COS(phiwc))
      rmu(1)=0.0
      Cmu(1)=1.0

      cukw=Cmu(1)*ubr/(kN*wr)
      IF ((cukw.gt.0.2).and.(cukw.le.100.0)) THEN       ! Eq 32/33
        fwci(1)=Cmu(1)*EXP(7.02*cukw**(-0.078)-8.82)
      ELSE IF ((cukw.gt.100.).and.(cukw.le.10000.0)) THEN
        fwci(1)=Cmu(1)*EXP(5.61*cukw**(-0.109)-7.30)
      ELSE IF (cukw.gt.10000.0 ) THEN
        fwci(1)=Cmu(1)*EXP(5.61*10000.0**(-0.109)-7.30)
      ELSE
        fwci(1)=Cmu(1)*0.43
      END IF
      ustrwm2(1)=0.5*fwci(1)*ubr*ubr                       ! Eq 29
      ustrr2(1)=Cmu(1)*ustrwm2(1)                             ! Eq 26
      ustrr=SQRT(ustrr2(1))
      IF (cukw.ge.8.0) THEN
        dwc(1)=2.0*vonKar*ustrr/wr                         ! Eq 36
      ELSE
        dwc(1)=kN
      END IF
      lnzr=LOG(zr/dwc(1))
      lndw=LOG(dwc(1)/zo)
      lnln=lnzr/lndw
      bigsqr=-1.0+SQRT(1.0+((4.0*vonKar*lndw)/                 &
     &                            (lnzr*lnzr))*ucr/ustrr)
      ustrci(1)=0.5*ustrr*lnln*bigsqr
!
      i=1
      diff=1.0
      DO WHILE ((i.lt.MAXIT).and.(diff.gt.0.000005))
        i=i+1
        rmu(i)=ustrci(i-1)*ustrci(i-1)/ustrwm2(i-1)
        Cmu(i)=SQRT(1.0+                                             &
     &              2.0*rmu(i)*cosphiwc+rmu(i)*rmu(i))     ! Eq 27
        cukw=Cmu(i)*ubr/(kN*wr)
        IF ((cukw.gt.0.2).and.(cukw.le.100.0)) THEN     ! Eq 32/33
          fwci(i)=Cmu(i)*EXP(7.02*cukw**(-0.078)-8.82)
        ELSE IF ((cukw.gt.100.).and.(cukw.le.10000.0)) THEN
          fwci(i)=Cmu(i)*EXP(5.61*cukw**(-0.109)-7.30)
        ELSE IF (cukw.gt.10000.0 ) THEN
          fwci(i)=Cmu(i)*EXP(5.61*10000.0**(-0.109)-7.30)
        ELSE
          fwci(i)=Cmu(i)*0.43
        END IF
        ustrwm2(i)=0.5*fwci(i)*ubr*ubr                     ! Eq 29
        ustrr2(i)=Cmu(i)*ustrwm2(i)                           ! Eq 26
        ustrr=SQRT(ustrr2(i))
!!      IF ((Cmu(1)*ubr/(kN*wr)).ge.8.0_r8) THEN  ! HGA Why 1?
        IF (cukw.ge.8.0) THEN
          dwc(i)=2.0*vonKar*ustrr/wr                       ! Eq 36
        ELSE
          dwc(i)=kN
        END IF
        lnzr=LOG(zr/dwc(i))
        lndw=LOG(dwc(i)/zo)
        lnln=lnzr/lndw
        bigsqr=-1.0+SQRT(1.0+((4.0*vonKar*lndw)/               &
     &                              (lnzr*lnzr))*ucr/ustrr)
        ustrci(i)=0.5*ustrr*lnln*bigsqr                    ! Eq 38
        diff=ABS((fwci(i)-fwci(i-1))/fwci(i))
      END DO
      ustrwm=SQRT(ustrwm2(i))
      ustrc=ustrci(i)
      ustrr=SQRT(ustrr2(i))
      phicwc=phiwc
      zoa=EXP(LOG(dwc(i))-(ustrc/ustrr)*LOG(dwc(i)/zo))       ! Eq 11
      fwc=fwci(i)

      RETURN
      END SUBROUTINE madsen94
# endif

End Module Mod_BBL
