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
! FVCOM Sediment Module 
!
! Copyright:    2005(c)
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
! Based on the Community Sediment Transport Model (CSTM) as implemented
!     in ROMS by J. Warner (USGS)
!
! Comments:     Sediment Dynamics Module 
!
! Current FVCOM dependency
!
!   init_sed.F:  - user defined sediment model initial conditions
!   mod_ncdio.F: - netcdf output includes concentration/bottom/bed fields
!   fvcom.F:     - main calls sediment setup 
!   internal.F:  - calls sediment advance
!
! History
!   Feb 7, 2008: added initialization of bottom(:,:) to 0 (w/ T. Hamada)
!              : fixed loop bounds in hot start and archive for conc (w/ T. Hamada)
!              : added comments describing theoretical bases of dynamics
!   Feb 14,2008: added non-constant settling velocity for cohesive sediments (w/ T. Hamada) 
!              : updated vertical flux routine to handle non-constant vertical velocity (w/ T. Hamada)
!              : added a user-defined routine to calculate settling velocity based on concentration (w/ T. Hamada)
!              : added a user-defined routine to calculate erosion for a general case (w/ T. Hamada)
!
!  PLEASE NOTE!!!!!!!!!!! 
!  Do NOT USE INTEL FORTRAN COMPILER VERSION 11.0 IT HAS KNOWN BUGS WHEN DEALING WITH TYPES WITH ALLOCATABLE
!    COMPONENTS.  YOU WILL SEE WEIRD BEHAVIOR.  VERSION 11.1 IS OK.
!   
!
!  Later
!   1.) Modify vertical flux routines to work with general vertical coordinate
!   2.) Add divergence term for bedload transport calc 
!   3.) Add ripple roughness calculation
!   4.) Add morphological change (bathymetry + vertical velocity condition) 
!   5.) Eliminate excess divisions and recalcs
!
!=======================================================================
Module Mod_Sed  
#if defined (SEDIMENT) && (ORIG_SED)
Use Mod_Par
Use Mod_Prec 
Use Mod_Types
Use Mod_wd
Use Control, only : seddis
Use all_vars,only : CNSTNT,UNIFORM,SEDIMENT_PARAMETER_TYPE
# if defined (WAVE_CURRENT_INTERACTION)
USE MOD_WAVE_CURRENT_INTERACTION
# endif

# if defined (FLUID_MUD)
USE MOD_FLUID_MUD
# endif

implicit none 

!--------------------------------------------------
!Sediment Type                         
!
! sname        => sediment name (silt,clay,etc)     
! stype        => sediment type: 'cohesive'/'non-cohesive'
! Sd50         => sediment mean diameter (mm)
! Wset         => mean sediment settling velocity (mm/s)
! tau_ce       => critical shear stress for erosion (N/m^2)
! tau_cd       => critical shear stress for deposition 
! Srho         => sediment density (kg/m^3)
! Spor         => sediment porosity (dimensionless)
! erate        => surface erosion mass flux         [kg m^-2 s^-1] 
! conc         => sed concentration in water column [kg m^-3]
! cnew         => sed concentration during update   [kg m^-3]
! mass         => sediment mass in bed layers       [kg m^-2]
! frac         => sediment fraction in bed layers   [-]
! bflx         => bedload sediment flux (kg/m^2)  (+out of bed) 
! eflx         => suspended sediment erosive flux (+out of bed)
! dflx         => suspended sed depositional flux (+into bed)
! cdis         => concentration at river source 
! cflx         => store advective flux at open bndry 
! cobc         => user specd open bndry concentration 
! depm         => store deposited mass
! arraysize    => spatial size of sediment arrays
!--------------------------------------------------
type sed_type
  character(len=20) :: sname    
  character(len=20) :: sname2    
  character(len=20) :: stype 
  real(sp)          :: Sd50
  real(sp)          :: Wset  
  real(sp)          :: w0   ! settling velocity without floccuation
  real(sp)          :: tau_ce
  real(sp)          :: tau_cd
  real(sp)          :: Srho 
  real(sp)          :: Spor 
! Jianzhong Ge 03/05/2013
  real(sp)          :: Chin ! concentration for hindered settling 
  real(sp)          :: Wrdc ! reduction scale for hindered settling velocity
! Jianzhong Ge 03/05/2013
  real(sp)          :: erate 
  real(dp)          :: cmax 
  real(dp)          :: cmin 
  real(dp)          :: crms 
  real(sp), allocatable :: conc(:,:)
  real(sp), allocatable :: cnew(:,:)
  real(sp), allocatable :: mass(:,:)
  real(sp), allocatable :: frac(:,:)
  real(sp), allocatable :: bflx(:)        
  real(sp), allocatable :: eflx(:)        
  real(sp), allocatable :: dflx(:)        
  real(sp), allocatable :: cdis(:)
  real(sp), allocatable :: cflx(:,:)
  real(sp), allocatable :: cobc(:)
  real(sp), allocatable :: depm(:)   
! Jianzhong Ge 03/05/2013
  real(sp), allocatable :: wset0(:,:) ! vertical settling velocity
! Jianzhong Ge 03/05/2013
!-------------------Jianzhong Ge-------------------
  real(sp), allocatable :: t_cd(:)   ! =tau_cd
  real(sp), allocatable :: t_ce(:)   ! =tau_ce
  real(sp), allocatable :: rate(:)   ! =erate
!--------------------------------------------------
  integer           :: arraysize
end type
 
public

!--------------------------------------------------
!Global Model Parameters
!
! sedfile   : sediment input parameter control file
! nsed      : number of sediment classes
! nbed      : number of layers in sediment bed
! min_Srho  : minimum Sediment density
! inf_bed   : true if bed has infinite sediment supply
! bedload   : true if bedload is to be considered
! susload   : true if suspended load is to be considered
! DTsed     : sediment model time step (seconds)
! T_model   : model time (seconds)
! taub      : array to hold bottom shear stress
! rho0      : mean density parameter used to convert tau/rho to tau
! thck_cr   : critical thickness for initiating new surface layer (m)
! n_report  : iteration interval for statistics printing
! sed_start : start interval for sed model
! sed_its   : sediment model iteration counter
! sed_nudge : flag for activiating sediment nudging on obc
! sed_alpha : sediment nudging relaxation factor 
! sed_ramp  : number of iterations over which to ramp sed nudging
! sed_source: flag for activiating sediment point sources
!--------------------------------------------------
character(len=120), parameter :: sed_model_version = "cstms 1.0: non-cohesive only"
character(len=120) :: sedfile
integer  :: nsed
integer  :: nbed
integer  :: Nobc
integer  :: Nobc_gl
real(sp) :: min_Srho
logical  :: inf_bed
logical  :: sed_nudge
real(sp) :: sed_alpha
integer  :: sed_ramp  
logical  :: sed_dumpbed
logical  :: sed_dumpbot 
logical  :: sed_source
!Jianzhong Ge 03/05/2013
logical  :: vert_hindered
!Jianzhong Ge 03/05/2013
logical  :: morpho_model = .false.
real(sp) :: morpho_factor = 1.0
integer  :: morpho_incr = 1
integer  :: morpho_strt = 0
real(sp) :: rho_water = 1025.
logical  :: bedload
logical  :: oned_model
logical  :: susload
integer  :: n_report
integer  :: sed_start
logical  :: sed_hot_start
integer  :: sed_its
real(sp) :: DTsed
real(sp) :: T_model
real(sp) :: settle_cfl = 1.0
integer  :: settle_limiter = 2   
real(sp) :: qlim  = 1.0
real(sp), allocatable :: taub(:)
real(sp), parameter   :: rho0 = 1025.
real(sp), parameter   :: thck_cr = .005 
!--------------Jianzhong Ge------------------------
real(sp), parameter   :: mf=0.012,Vcons=1.e-6
real(sp), allocatable :: csed(:,:)
real(sp) :: deposition_d50
real(sp) :: sed_rhos
real(sp) :: sed_wset
!--------------------------------------------------
logical, parameter    :: debug_sed = .false.

!--------------------------------------------------
!Initial condition parameters
!  Read in sediment input file
!  Can be overwritten in init_sed.F or through 
!  restart
!
!  init_bed_porosity:  initial bed porosity [0,1]
!  init_bed_thickness: initial bed thickness [m] [> 0]
!  init_bed_fraction:  initial bed fraction [-] [> 0]
!--------------------------------------------------
 
real(sp) :: init_bed_porosity  = 1.0
real(sp), allocatable :: init_bed_thickness(:) 
real(sp), allocatable :: init_bed_fraction(:) 

!--------------------------------------------------
!Bedload 
!
!Meyer-Peter Muller Bedload Formulation Parameters
!
! Shield_Cr_MPM => Effective Critical Shields number
! Gamma_MPM     => MPM Power Law Coefficient
! k_MPM         => MPM Multiplier
!
!Generic Params for all bedload
!
! bedload_rate  => bedload rate coefficient
!
!--------------------------------------------------

real(sp) :: Shield_Cr_MPM  = 0.047  !  default: 0.047 
real(sp) :: Gamma_MPM      = 1.5    !  default: 1.5
real(sp) :: k_MPM          = 8.0    !  default: 8.0
real(sp) :: bedload_rate   = 0.1    !  default: 0.1
logical  :: bedload_smooth = .false.!  default: 0.1


!--------------------------------------------------
!Bed 
! data --> bed(horizontal,vertical,proptype)
! size --> bed(0:mt , nbed , n_bed_chars))
!
! with proptype=:
!   ithck  => layer thickness
!   iaged  => layer age 
!   iporo  => layer porosity  [volume of voids/total volume]
!
! global bed characteristics 
!
! n_bed_chars => number of bed characteristics 
!-------------------------------------------------

integer, parameter    :: n_bed_chars = 3
!character(len=80), dimension(n_bed_chars) :: bed_snames = [character(len=80) :: & 
!     "bed_thick","bed_age","bed_por"]
!character(len=80), dimension(n_bed_chars) :: bed_lnames = [character(len=80) :: & 
!     "bed thickness","bed age","bed porosity"]
!character(len=80), dimension(n_bed_chars) :: bed_units =  [character(len=80) :: & 
!     "m","days","-"]
character(len=80), dimension(n_bed_chars) :: bed_snames = (/"bed_thick","bed_age","bed_por"/)
character(len=80), dimension(n_bed_chars) :: bed_lnames = (/"bed thickness","bed age","bed porosity"/)
character(len=80), dimension(n_bed_chars) :: bed_units = (/"m","days","-"/)

integer, parameter    :: ithck = 1
integer, parameter    :: iaged = 2
integer, parameter    :: iporo = 3
real(sp), allocatable, target :: bed(:,:,:)

!--------------------------------------------------
!Bottom (Exposed Sediment Layer)
! data --> bed(0:mt,n_bot_vars)
!
! with proptype=:
!   isd50  => mean grain diameter 
!   isphi  => mean grain diameter in phi values
!   idens  => mean grain density
!   iwset  => mean settling velocity
!   itauc  => mean critical erosion stress [Pa, N/m^2]
!   iactv  => active layer thickness            
!   nthck  => new total thickness of sediment layer
!   lthck  => last thickness of sediment layer 
!   dthck  => accumulated delta of layer thickness [m]
!   tmass  => total mass in sediment layer [kg/m^2]
!   morph  => differential in thickness of sediment layer
!
! global bottom characteristics 
!   n_bot_vars => number of bottom characteristics
!
!-------------------------------------------------


!  BOTTOM properties indices:                                          !
!                                                                      !
!   MBOTP    Number of bottom properties (array dimension).            !
!   idBott   IO indices for bottom properties variables.               !
!   isd50    Median sediment grain diameter (m).                       !
!   isphi    mean grain diameter in phi values (-)                     !
!   idens    Median sediment grain density (kg/m3).                    !
!   iwsed    Mean settling velocity (m/s).                             !
!   itauc    Mean critical erosion stress (m2/s2).                     !
!   irlen    Sediment ripple length (m).                               !
!   irhgt    Sediment ripple height (m).                               !
!   ibwav    Bed wave excursion amplitude (m).                         !
!   izdef    Default bottom roughness (m).                             !
!   izapp    Apparent bottom roughness (m).                            !
!   izNik    Nikuradse bottom roughness (m).                           !
!   izbio    Biological bottom roughness (m).                          !
!   izbfm    Bed form bottom roughness (m).                            !
!   izbld    Bed load bottom roughness (m).                            !
!   izwbl    Bottom roughness used wave BBL (m).                       ! 
!   iactv    Active layer thickness for erosive potential (m).         !
!   ishgt    Sediment saltation height (m).                            !
!   idefx    Erosion flux                                              !
!   idnet    Erosion or deposition                                     !
!   idoff    Offset for calculation of dmix erodibility profile (m)    !
!   idslp    Slope for calculation of dmix or erodibility profile      !
!   idtim    Time scale for restoring erodibility profile (s)          !
!   idbmx    Bed biodifusivity maximum                                 !
!   idbmm    Bed biodifusivity minimum                                 !
!   idbzs    Bed biodifusivity zs                                      !
!   idbzm    Bed biodifusivity zm                                      !
!   idbzp    Bed biodifusivity phi                                     !
!   idprp    Cohesive behavior                                         !
!   bflag    Bed Flag (=0, erosion/deposition disabled)                !
!                                                                      !
!=======================================================================


integer, parameter    :: n_bot_vars  = 22
!character(len=80), dimension(n_bot_vars) :: bot_snames = [character(len=80) :: &
!      "bot_sd50","bot_dens","bot_wset","bot_tauc","bot_actv","bot_nthck",&
!      "bot_lthck","bot_dthck","bot_tmass","morph","rip_length","&
!       rip_height","wav_excur","def_rough","app_rough","Niku_rough","b&
!       io_rough","bfm_rough","bld_rough","wave_bbl","bot_phi","bedflag"]
!
!character(len=80), dimension(n_bot_vars) :: bot_lnames = [character(len=80) :: &
!      "bot_sd50","bot_dens","bot_wset","bot_tauc","bot_actv","sedlayer_depth_current_iteration",& 
!      "sedlayer_dapth_last_iteration","sedlayer_depth_net_change","bot_tmass","morph","rip_length","&
!       rip_height","wav_excur","def_rough","app_rough","Niku_rough","b&
!       io_rough","bfm_rough","bld_rough","wave_bbl","bot_phi","bedflag"]
!
!character(len=80), dimension(n_bot_vars) :: bot_units = [character(len=80) :: &  
!       "m","kg/m^3","m/s","N/m^2","days","m","m","m","kg","m","m","m","m","m"&
!      ,"m","m","m","m","m","m","-","-"] 
character(len=80), dimension(n_bot_vars) :: bot_snames = &
      (/"bot_sd50","bot_dens","bot_wset","bot_tauc","bot_actv","bot_nthck",&
      "bot_lthck","bot_dthck","bot_tmass","morph","rip_length","&
      &rip_height","wav_excur","def_rough","app_rough","Niku_rough","b&
      &io_rough","bfm_rough","bld_rough","wave_bbl","bot_phi","bedflag"/)
character(len=80), dimension(n_bot_vars) :: bot_lnames = &
      (/"bot_sd50","bot_dens","bot_wset","bot_tauc","bot_actv","sedlayer_depth_current_iteration",&
      "sedlayer_dapth_last_iteration","sedlayer_depth_net_change","bot_tmass","morph","rip_length","&
      &rip_height","wav_excur","def_rough","app_rough","Niku_rough","b&
      &io_rough","bfm_rough","bld_rough","wave_bbl","bot_phi","bedflag"/)
character(len=80), dimension(n_bot_vars) :: bot_units = (/"m","kg&
     &/m^3","m/s","N/m^2","days","m","m","m","kg","m","m","m","m","m"&
     &,"m","m","m","m","m","m","-","-"/)

integer, parameter    :: isd50  = 1  
integer, parameter    :: idens  = 2
integer, parameter    :: iwset  = 3
integer, parameter    :: itauc  = 4
integer, parameter    :: iactv  = 5 
integer, parameter    :: nthck  = 6 
integer, parameter    :: lthck  = 7 
integer, parameter    :: dthck  = 8 
integer, parameter    :: tmass  = 9 
integer, parameter    :: morph  = 10
integer, parameter    :: irlen  = 11
integer, parameter    :: irhgt  = 12
integer, parameter    :: ibwav  = 13
integer, parameter    :: izdef  = 14
integer, parameter    :: izapp  = 15
integer, parameter    :: izNik  = 16
integer, parameter    :: izbio  = 17
integer, parameter    :: izbfm  = 18
integer, parameter    :: izbld  = 19
integer, parameter    :: izwbl  = 20
integer, parameter    :: isphi  = 21
integer, parameter    :: bflag  = 22

real(sp), allocatable, target :: bottom(:,:)

!--------------------------------------------------
!Sediment Point Source Data 
!   sbc_tm:  time map for source data
!   seddis:  source data
!--------------------------------------------------

!--------------------------------------------------
!Sediment Array                  
!--------------------------------------------------

type(sed_type), allocatable, target :: sed(:)

!J. Ge for tracer advection
type(sed_type), allocatable, target :: sed0(:),sed2(:)
!J. Ge for tracer advection

contains


!==========================================================================
! Allocate Data, Initialize Concentration Fields and Bed Parameters       
!==========================================================================
  Subroutine Setup_Sed(fname,Nin,Nin_gl,N_sed,N_sed_Max,Sed_Names) 
  use control, only : msr,casename,input_dir
  use lims,    only : m
  Use Mod_Utils, only: pstop
  implicit none
  integer, intent(in)  :: Nin,Nin_gl   
  character(len=80)    :: fname
  integer, intent(out) :: N_sed
  integer, intent(in ) :: N_sed_max   
  character(len=20)    :: Sed_Names(N_sed_max)
  !-Local--------------------------
  logical             :: fexist
  integer             :: ised
  integer             :: i,k

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: setup_sed" 

  Nobc = Nin
  Nobc_gl = Nin_gl
!----------------------------------------------------------------
! check arguments and ensure necessary files exist
!----------------------------------------------------------------

  !initialize sediment model iteration counter
  sed_its = 0

  !ensure sediment parameter file exists
  if(fname == "none" .or. fname == "NONE" .or. fname == "None")then
    sedfile = "./"//trim(input_dir)//"/"//trim(casename)//'_sediment.inp'
  else
    sedfile = "./"//trim(input_dir)//"/"//trim(fname)
  endif
    
  inquire(file=trim(sedfile),exist=fexist)
  if(.not.fexist)then
    write(*,*)'sediment parameter file: ',trim(sedfile),' does not exist'
    write(*,*)'stopping'
    call pstop
  end if

!----------------------------------------------------------------
! read in sediment parameter file                    
!----------------------------------------------------------------
  call read_sed_params

! J. Ge for fluid mud layer
# if defined (FLUID_MUD)
!----------------------------------------------------------------
! read in fluid mud layer parameter file
!----------------------------------------------------------------
  call initialize_fluid_mud
# endif
! J. Ge for fluid mud layer


!----------------------------------------------------------------
! allocate data                                      
!----------------------------------------------------------------
  call alloc_sed_vars

!----------------------------------------------------------------
! initialize fields                                  
!   these will be overwritten if this is a restart case
!----------------------------------------------------------------
  call init_sed 
!----------------------------------------------------------------
! setup open boundary condition nudging values       
!----------------------------------------------------------------
  if(sed_nudge) call setup_sed_obc 

!----------------------------------------------------------------
! setup point source river forcing data 
!----------------------------------------------------------------
!  if(numqbc_gl .and. sed_source) call setup_sed_PTsource
!  if(sed_source) call setup_sed_PTsource

  !--------------------------------------------------
  !Initialize Bed_Mass properties
  !--------------------------------------------------
   Do k=1,Nbed
     Do i=1,m
       Do ised=1,Nsed
         sed(ised)%mass(i,k) = bed(i,k,ithck)*       &
                              (1.0-bed(i,k,iporo))*  &
                               sed(ised)%Srho*       &
                               sed(ised)%frac(i,k)
       end do
     end do
  end do

  !--------------------------------------------------
  !Calculate Bottom Thicknesses
  !--------------------------------------------------

  do i=1,m
    do k=1,nbed
      bottom(i,nthck) = bottom(i,nthck) + bed(i,k,ithck)
    end do
  end do


!----------------------------------------------------------------
! report sediment parameters and statistics to screen
!----------------------------------------------------------------
  call report_sed_setup

!----------------------------------------------------------------
! convert parameters to MKS 
!----------------------------------------------------------------
  call convert_sed_params 

!----------------------------------------------------------------
! update bottom properties                           
!----------------------------------------------------------------
  call update_bottom_properties

!----------------------------------------------------------------
! set arguments for sed names and N_sed for FVCOM to use 
! in setting up river forcing
!----------------------------------------------------------------
  N_sed = Nsed
  do i=1,Nsed
    sed_names(i) = sed(i)%sname
  end do
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: setup_sed" 

End Subroutine Setup_Sed

!=======================================================================
!Convert Sed Params from nonstandard (mm, etc) to mks units
!=======================================================================
  Subroutine Convert_Sed_Params 
  Implicit None
  integer :: i

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Convert_Sed_Params" 

  !convert settling velocity to m/s from input mm/s
  do i=1,nsed
    sed(i)%Wset = sed(i)%Wset*.001
    sed(i)%w0   = sed(i)%w0*0.001
  end do

  !convert mean grain diameter from mm to m 
  do i=1,nsed
    sed(i)%Sd50 = sed(i)%Sd50*.001
  end do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Convert_Sed_Params" 

  End Subroutine Convert_Sed_Params

!=======================================================================
!Read Sediment Parameters From Sediment Input File
!=======================================================================
  Subroutine Read_Sed_Params
  Use Input_Util
  Implicit None
  integer linenum,i,k1,iscan
  real(sp)           :: ftemp
  character(len=120) :: stemp

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Read_Sed_Params " 

  !check if we are running 1D model
  Call Get_Val(oned_model,sedfile,'SED_ONED',line=linenum,echo=.false.)

  !read in number of sediment classes 
  Call Get_Val(nsed,sedfile,'NSED',line=linenum,echo=.false.)  

  !read in start interval for sed model 
  Call Get_Val(sed_start,sedfile,'SED_START',line=linenum,echo=.false.)  

  !hot start option
  Call Get_Val(sed_hot_start,sedfile,'SED_HOT_START',line=linenum,echo=.false.)  

  !read in morphology parameters
  Call Get_Val(morpho_model,sedfile,'MORPHO_MODEL',line=linenum,echo=.false.)
  Call Get_Val(morpho_factor,sedfile,'MORPHO_FACTOR',line=linenum,echo=.false.)  
  if(morpho_model)then
    Call Get_Val(morpho_incr,sedfile,'MORPHO_INCR',line=linenum,echo=.false.)  
    Call Get_Val(morpho_strt,sedfile,'MORPHO_STRT',line=linenum,echo=.false.)  
  else
    morpho_incr   = 1
  endif
  
  !deactivate morpho model - not ready
!  morpho_model = .false.
!  morpho_factor = 1.0

  !read in interation interval for reporting 
  Call Get_Val(n_report,sedfile,'N_REPORT',line=linenum,echo=.false.)  

  !read in number of bed layers 
  Call Get_Val(nbed,sedfile,'NBED',line=linenum,echo=.false.)

  !read in logical for infinite sediment supply 
  Call Get_Val(inf_bed,sedfile,'INF_BED',line=linenum,echo=.false.)

  !read in initial bed porosity  
  Call Get_Val(init_bed_porosity,sedfile,'INIT_BED_POROSITY',line=linenum,echo=.false.)
  if(init_bed_porosity < 0. .or. init_bed_porosity > 1)then
    write(*,*)'error in init_bed_porosity in sed param file'
    write(*,*)'you entered: ',init_bed_porosity
    write(*,*)'must be in range [0,1]'
  endif

  !read in initial bed thickness
  allocate(init_bed_thickness(nbed))
  if(nbed==1)then
    Call Get_Val(ftemp,sedfile,'INIT_BED_THICKNESS',line=linenum,echo=.false.)
    init_bed_thickness(1) = ftemp
  else
    open(unit=999,file=trim(sedfile),form='formatted')
    iscan = scan_file(999,"INIT_BED_THICKNESS",fvec = init_bed_thickness,nsze = nbed)
    if(iscan /= 0)then
      write(*,*)'problem reading INIT_BED_THICKNESS from sediment param file'
      stop
    endif
    close(999)
  endif
  if(minval(init_bed_thickness) <= 0.)then
    write(*,*)'error in init_bed_thickness in sed param file'
    write(*,*)'you entered: ',init_bed_thickness
    write(*,*)'must be > 0'
  endif

  !read in initial bed fraction   
  allocate(init_bed_fraction(nsed))
  if(nsed==1)then
    Call Get_Val(ftemp,sedfile,'INIT_BED_FRACTION',line=linenum,echo=.false.)
    init_bed_fraction(1) = ftemp
  else
    open(unit=999,file=trim(sedfile),form='formatted')
    iscan = scan_file(999,"INIT_BED_FRACTION",fvec = init_bed_fraction,nsze = nsed)
    if(iscan /= 0)then
      write(*,*)'problem reading INIT_BED_FRACTION from sediment param file'
      stop
    endif
    close(999)
  endif
  if(minval(init_bed_fraction) < 0. .or. maxval(init_bed_fraction) > 1.)then
    write(*,*)'error in init_bed_fraction in sed param file'
    write(*,*)'also error in your chosen career'
    write(*,*)'you entered: ',init_bed_fraction  
    write(*,*)'must be >= 0 and <= 1'
  endif

  !read in logical for bedload calculation 
  Call Get_Val(bedload,sedfile,'BEDLOAD',line=linenum,echo=.false.)

  !shut down bedload if model is oned
  if(oned_model) bedload = .false.

  !read in logical for suspended load calculation 
  Call Get_Val(susload,sedfile,'SUSLOAD',line=linenum,echo=.false.)

  !read in minumum sediment density
  Call Get_Val(min_Srho,sedfile,'MIN_SRHO',line=linenum,echo=.false.)

  !read in nudging switch 
  Call Get_Val(sed_nudge,sedfile,'SED_NUDGE',line=linenum,echo=.false.)

  !read in params controlling output 
  Call Get_Val(sed_dumpbed,sedfile,'SED_DUMPBED',line=linenum,echo=.false.)
  Call Get_Val(sed_dumpbot,sedfile,'SED_DUMPBOT',line=linenum,echo=.false.)

  !read in nudging relaxation factor 
  if(sed_nudge) then
  Call Get_Val(sed_alpha,sedfile,'SED_ALPHA',line=linenum,echo=.false.)
  Call Get_Val(sed_ramp ,sedfile,'SED_RAMP ',line=linenum,echo=.false.)
  endif

  !settling CFL
  Call Get_Val(settle_cfl ,sedfile,'SETTLE_CFL',line=linenum,echo=.false.)
  if(settle_cfl <= 0 .or. settle_cfl > 1.)then
    write(*,*)'error in mod_sed, settle_cfl must be in range [>0,1]'
    write(*,*)'you entered: ',settle_cfl
    stop
  endif

  !limiter type
  Call Get_Val(settle_limiter,sedfile,'SETTLE_LIMITER',line=linenum,echo=.false.)
  if(settle_limiter < 1 .or. settle_limiter > 3)then
    write(*,*)'error in mod_sed, settle_limiter must be [1,2,3]'
    write(*,*)'you entered: ',settle_limiter
    stop
  endif
  if(settle_limiter==1) qlim = 0.
  if(settle_limiter==2) qlim = 1.
  if(settle_limiter==3) qlim = 2.
 
  !hindered effect on vertical diffusion
  Call Get_Val(VERT_HINDERED,sedfile,'VERT_HINDERED',line=linenum,echo=.false.)

  !read in point source switch 
  Call Get_Val(sed_source,sedfile,'SED_PTSOURCE',line=linenum,echo=.false.)

  !check values
  if(nsed < 1)then
    write(*,*)'sediment input file must have at least one sediment class'
    write(*,*)'currently nsed = ',nsed
    stop
  endif
  if(nbed < 1)then
    write(*,*)'sediment input file must have at least one bed layer'
    write(*,*)'currently nbed = ',nbed
    stop
  endif
  if(morpho_factor < 0. )then
    write(*,*)'sediment morpho_factor must be positive'
    write(*,*)'currently is = ',morpho_factor
    stop
  endif
  if(morpho_incr < 0 )then
    write(*,*)'sediment morpho_incr must be positive'
    write(*,*)'currently is = ',morpho_incr  
    stop
  endif

  if(morpho_strt < 0 )then
    write(*,*)'sediment morpho_strt must be non-negative'
    write(*,*)'currently is = ',morpho_strt  
    stop
  endif

  !allocate sediment data space
  Allocate(sed(nsed))

!J. Ge for tracer advection
  Allocate(sed0(nsed))
  Allocate(sed2(nsed))
!J. Ge for tracer advection


  !read in sediment parameters
  k1 = 1

  do i=1,nsed

    !read SED_NAME and mark position
    Call Get_Val(stemp,sedfile,'SED_NAME',line=linenum,echo=.false.,start=k1)
    sed(i)%sname = stemp
    sed(i)%sname2 = trim(sed(i)%sname)//'_bload'
    k1        = linenum+1

    !read type
    Call Get_Val(stemp,sedfile,'SED_TYPE',line=linenum,echo=.false.,start=k1)
    sed(i)%stype = stemp

    !read mean diameter
    Call Get_Val(ftemp,sedfile,'SED_SD50',line=linenum,echo=.false.,start=k1)
!Jianzhong GE
    sed(i)%Sd50 = ftemp
    deposition_d50 = ftemp*0.001
!------------------Jianzhong---------------
!    sed(i)%w0 = -4.*4.27*Vcons/(ftemp*1.0e-3)/1.22+ &
!                sqrt((4.*4.27/1.22*Vcons/(ftemp*1.0e-3))**2+ &
!                4./(3.*1.22)*1.65*9.8*(ftemp*1.0e-3)  )
!------------------------------------------
    !read sediment density 
    Call Get_Val(ftemp,sedfile,'SED_SRHO',line=linenum,echo=.false.,start=k1)
    sed(i)%Srho = ftemp 
!Jianzhong GE
    sed_rhos = ftemp

    !read sediment settling rate 
    Call Get_Val(ftemp,sedfile,'SED_WSET',line=linenum,echo=.false.,start=k1)
    sed(i)%Wset = ftemp
    sed(i)%w0   = ftemp
!Jianzhong GE
    sed_wset = ftemp*0.001

    !read sediment surface erosion rate 
    Call Get_Val(ftemp,sedfile,'SED_ERAT',line=linenum,echo=.false.,start=k1)
    sed(i)%erate = ftemp

    !read sediment critical erosive shear stress 
    Call Get_Val(ftemp,sedfile,'SED_TAUE',line=linenum,echo=.false.,start=k1)
    sed(i)%tau_ce = ftemp

    !read sediment critical depositional shear stress 
    Call Get_Val(ftemp,sedfile,'SED_TAUD',line=linenum,echo=.false.,start=k1)
    sed(i)%tau_cd = ftemp

    !read sediment porosity
    Call Get_Val(ftemp,sedfile,'SED_PORS',line=linenum,echo=.false.,start=k1)
    sed(i)%Spor = ftemp

!Jianzhong Ge 03/05/2013
    if(VERT_HINDERED)then
      !read concentration for hindered settling  
      if(sed(i)%stype=='cohesive')then
        Call Get_Val(ftemp,sedfile,'SED_CHIN',line=linenum,echo=.false.,start=k1)
        sed(i)%Chin = ftemp
      end if

      !read reduction scale for settling velocity
      if(sed(i)%stype=='cohesive')then
        Call Get_Val(ftemp,sedfile,'SED_WRDC',line=linenum,echo=.false.,start=k1)
        sed(i)%Wrdc = ftemp
      end if
    end if
!Jianzhong Ge 03/05/2013

  end do

  ! read in bedload function parameters
  Call Get_Val(Shield_Cr_MPM,sedfile,'MPM_CS',line=linenum,echo=.false.)
  Call Get_Val(    Gamma_MPM,sedfile,'MPM_GM',line=linenum,echo=.false.)
  Call Get_Val(        k_MPM,sedfile,'MPM_K ',line=linenum,echo=.false.)
  Call Get_Val(bedload_rate ,sedfile,'BEDLOAD_RATE ',line=linenum,echo=.false.)
  Call Get_Val(bedload_smooth ,sedfile,'BEDLOAD_SMOOTH',line=linenum,echo=.false.)
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Read_Sed_Params " 
 
  End Subroutine Read_Sed_Params 
     
!=======================================================================
! Advance Sediment Model by time step DTSED 
!=======================================================================
  Subroutine Advance_Sed(DTin,Tin,taub_in)  
  Use Input_Util
  Use Scalar
  Use Control, only : ireport,iint,msr,par
  Use Lims,    only : m,mt,nt,kbm1,numqbc,kb,nprocs,myid
  Use All_Vars,only : BACKWARD_ADVECTION,BACKWARD_STEP
!  Use Mod_OBCS,only : iobcn
# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif

  Implicit None
  real(sp), intent(in ) :: DTin,Tin
  real(sp), intent(in ) :: taub_in(m)

  integer :: i,k,ised,l1,l2,ierr,d_cdis,d_cflx
  character(len=4) :: fnum
  real(sp) :: fact,ufact,temp

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Advance_Sed "

  !------------------------------------------------------
  !Check for Initialization
  !------------------------------------------------------
  if(iint < sed_start) return

  !------------------------------------------------------
  !Increment Sed Model Counter and Report Sed Model Init 
  !------------------------------------------------------
  sed_its = sed_its + 1
  if(msr .and.  sed_its == 1)then
    write(*,*)'========Sed Model Initiated======='
  endif

  !------------------------------------------------------
  !set bed age to current model time if first sed call 
  !------------------------------------------------------
  if(sed_its ==1)then
    bed(:,:,iaged) = T_model
  endif

  !------------------------------------------------------
  !Set up Model Forcing  -> time step and bottom shear 
  !Convert taux/tauy to node-based values
  !Convert Model tau/rho to tau [N/m^2]
  !------------------------------------------------------
  DTsed = DTin
  T_model = Tin
  taub(1:m) = taub_in(1:m)*rho_water

  !------------------------------------------------------
  !Calculate Active Layer Thickness --> bottom(:,iactv)
  !------------------------------------------------------
  call calc_active_layer

  !------------------------------------------------------
  !Vertical Sediment Dynamics --> Sediment 'Model'
  !------------------------------------------------------
  if(bedload)then
    !calculate bedload flux => sed(1:nsed)%bflx 
    !call calc_bedload_flux
    call calc_bedload_flux2

    !adjust surface bed properties
    call update_surface_bed('bedload')

    !add new surface layer if necessary
    call add_new_layer

    !update bed layer stratigraphy
    call update_stratigraphy

    !adjust bottom properties
    call update_bottom_properties
  endif

  if(susload)then
    !calculated depositional flux from settling
    call calc_deposition

    !calculate erosive flux
    call calc_erosion

! J. Ge for fluid mud layer

# if defined (FLUID_MUD)

!# if defined (WAVE_CURRENT_INTERACTION)
!    call fluid_mud_source_sink(DTsed,sed(1)%Srho,bottom(:,bflag),bed(:,1,ithck),bed(:,1,iporo),sed(1)%frac(:,1),bustrcwmax,bvstrcwmax,,sed(1)%t_ce,sed(1)%rate)
!# else
!    call fluid_mud_source_sink(DTsed,sed(1)%Srho,bottom(:,bflag),bed(:,1,ithck),bed(:,1,iporo),sed(1)%frac(:,1),wubot,wvbot,sed(1)%t_ce,sed(1)%rate)
!# endif
    call fluid_mud_source_sink(DTsed,sed(1)%Srho,bottom(:,bflag),bed(:,1,ithck),bed(:,1,iporo),sed(1)%frac(:,1),taub,sed(1)%t_ce,sed(1)%rate)

    call internal_step_fluid_mud

    do ised=1,nsed
      do i=1,m
# if defined (WET_DRY)
        if(iswetn(i)==1)then
# endif
          sed(ised)%conc(i,kbm1) = sed(ised)%conc(i,kbm1) + Entrainment(i)/(d(i)*dz(i,kbm1))
          sed(ised)%dflx(i) = dewatering(i)
          sed(ised)%eflx(i) = BedErosion(i) 
# if defined (WET_DRY)
        else
          sed(ised)%conc(i,kbm1) = 0.0_sp
          sed(ised)%dflx(i) = 0.0_sp
          sed(ised)%eflx(i) = 0.0_sp
        end if
# endif
      end do
    end do

# endif

    !adjust surface bed properties
    call update_surface_bed('susload')

    !add new surface layer if necessary
    call add_new_layer

    !update bed layer stratigraphy
    call update_stratigraphy

    !adjust bottom properties
    call update_bottom_properties
  endif
   
  !------------------------------------------------------
  !Horizontal Advection of Sediment Concentration 
  !------------------------------------------------------
  if(.not. oned_model)then
    d_cdis = 1
     if(numqbc > 0) d_cdis = numqbc
     d_cflx = 1
     if(Nobc > 0) d_cflx = Nobc
      
     do i=1,nsed
      
      !set discharge sediment concentrations
      if(numqbc > 0 .and. sed_source)then
        sed(i)%cdis(1:numqbc) = seddis(1:numqbc,i)
      endif
  
      !3D advection

!      call adv_scal(sed(i)%conc, &
!                      sed(i)%cnew, &
!                      d_cdis,      &
!                      sed(i)%cdis(1:numqbc), &
!                      d_cflx,      &
!                      sed(i)%cflx, &
!                      DTsed ,sed_source)
  

!J. Ge for tracer advection
      call adv_scal(sed(i)%conc,sed0(i)%conc,sed2(i)%conc, &
                      sed(i)%cnew, &
                      d_cdis,      &
                      sed(i)%cdis(1:numqbc), &
                      d_cflx,      &
                      sed(i)%cflx, &
                      DTsed ,sed_source)

!J. Ge for tracer advection


      call fct_sed(sed(i)%conc,sed(i)%cnew)
  
    end do
  else
   do i=1,nsed
     sed(i)%cnew = sed(i)%conc
   end do
  endif 


  !------------------------------------------------------
  !Vertical Diffusion of Sediment Concentrations
  !------------------------------------------------------

!Jianzhong GE 03/05/2013
  if(vert_hindered)call hindered_diffusion
!Jianzhong GE 03/05/2013

  do i=1,nsed
    call vdif_scal(sed(i)%cnew,DTsed) 
  end do

  !------------------------------------------------------
  !Exchange Concentration on Processor Boundaries
  !------------------------------------------------------
# if defined (MULTIPROCESSOR)
  if(par)then
    do i=1,nsed
      call aexchange(nc,myid,nprocs,sed(i)%cnew)
    end do
  endif
# endif

  !------------------------------------------------------
  ! Set Open Boundary Conditions (nudging to background)
  !------------------------------------------------------
  if(.not. oned_model)then
    if(sed_ramp == 0) sed_ramp = 1
    temp = min(float(sed_its/sed_ramp), 1.0)*sed_alpha
    if(.not.sed_nudge) temp = 0.0
    do i=1,nsed 
       call bcond_scal_OBC(sed(i)%conc, &
                           sed(i)%cnew, &
                           sed(i)%cflx, &
                           sed(i)%cobc, &
                           DTsed,       &
                           temp) 
    end do
  endif

  !------------------------------------------------------
  ! Point Source Boundary Conditions 
  !------------------------------------------------------
  if(sed_source .and. .not. oned_model)then
    do i=1,nsed
      call bcond_scal_PTsource(sed(i)%conc, &
                               sed(i)%cnew, &
                               sed(i)%cdis) 
    end do
  endif

  !------------------------------------------------------
  ! Update Concentration Variables to next time level
  !------------------------------------------------------
  do i=1,nsed
    sed(i)%conc = sed(i)%cnew
  end do

  !------------------------------------------------------
  ! Limit to Non-Negative (modify to Venkat Limiter) 
  !------------------------------------------------------
  do i=1,nsed
    sed(i)%conc = max(sed(i)%cnew,0.0)
  end do
  
  !-----------------------Jianzhong----------------------  
  csed = 0.0
  do i=1,nsed
    csed = csed + sed(i)%conc
  end do
  !------------------------------------------------------ 
  
  !------------------------------------------------------
  !Exchange Concentration on Processor Boundaries
  !------------------------------------------------------
# if defined (MULTIPROCESSOR)
  if(par)then
    do i=1,nsed
      call aexchange(nc,myid,nprocs,sed(i)%cnew)
    end do
    call aexchange(nc,myid,nprocs,taub)
  endif
# endif

!J. Ge for tracer advection
  IF(BACKWARD_ADVECTION==.TRUE.)THEN
    IF(BACKWARD_STEP==1)THEN
      do i=1,nsed
        SED0(i)%conc = sed(i)%cnew
      end do
    ELSEIF(BACKWARD_STEP==2)THEN
      do i=1,nsed
        SED2(i)%conc  = SED0(i)%conc 
        SED0(i)%conc  = sed(i)%cnew 
      end do
    ENDIF
  ENDIF
!J. Ge for tracer advection

  !------------------------------------------------------
  ! Update Depth Delta and Morphology Array
  ! this array is positive if material has accumulated,
  ! negative if eroded
  !------------------------------------------------------
  call update_thickness_delta

  !------------------------------------------------------
  ! Report 
  !------------------------------------------------------
  if(mod(sed_its,n_report)==0)call sed_report    

!  call layer_report(1)
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Advance_Sed "

  End Subroutine Advance_Sed

!=======================================================================
! Allocate Sediment Variables                      
!=======================================================================
  Subroutine Alloc_Sed_Vars
  use lims,     only: nt,mt,kb,numqbc,kbm1
!  use mod_obcs, only: iobcn
  implicit none
  integer i,tmp1,tmp2,tmp3

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Alloc_Sed_Vars "

  !ensure arrays have nonzero dimension
  tmp1 = max(numqbc,1)
  tmp2 = max(Nobc+1,1)
  tmp3 = max(Nobc  ,1)

  !allocate suspended sediment arrays
  allocate(csed(0:mt,kb)) ; csed = 0.0

  do i=1,nsed
    sed(i)%arraysize = mt+1
    allocate(sed(i)%conc(0:mt,kb  ))      ; sed(i)%conc = 0.0
    allocate(sed(i)%cnew(0:mt,kb  ))      ; sed(i)%cnew = 0.0
    allocate(sed(i)%mass(0:mt,nbed))      ; sed(i)%mass = 0.0
    allocate(sed(i)%frac(0:mt,nbed))      ; sed(i)%frac = 0.0
    allocate(sed(i)%bflx(0:mt     ))      ; sed(i)%bflx = 0.0
    allocate(sed(i)%eflx(0:mt     ))      ; sed(i)%eflx = 0.0
    allocate(sed(i)%dflx(0:mt     ))      ; sed(i)%dflx = 0.0
    allocate(sed(i)%cdis(500      ))      ; sed(i)%cdis = 0.0
    allocate(sed(i)%cflx(tmp2,kbm1))      ; sed(i)%cflx = 0.0
    allocate(sed(i)%cobc(tmp3     ))      ; sed(i)%cobc = 0.0
    allocate(sed(i)%depm(0:mt     ))      ; sed(i)%depm = 0.0
! Jianzhong Ge 03/05/2013
    allocate(sed(i)%wset0(0:mt,kb  ))     ; sed(i)%wset0 = 0.0
    sed(i)%wset0 = sed_wset
! Jianzhong Ge 03/05/2013

!    if(SEDIMENT_PARAMETER_TYPE/=UNIFORM)then
    allocate(sed(i)%t_cd(0:mt     ))      ; sed(i)%t_cd = 0.0
    allocate(sed(i)%t_ce(0:mt     ))      ; sed(i)%t_ce = 0.0
    allocate(sed(i)%rate(0:mt     ))      ; sed(i)%rate = 0.0
!    end if

!J. Ge for tracer advection
!    if(backward_advection==.true.)then
      allocate(sed0(i)%conc(0:mt,kb  ))     ; sed0(i)%conc = 0.0
      allocate(sed2(i)%conc(0:mt,kb  ))     ; sed2(i)%conc = 0.0
!    endif
!J. Ge for tracer advection
  end do

  !allocate bottom shear stress array 
  allocate(taub(0:mt)) ; taub = 0.0

  !allocate bed data
  allocate(bed(0:mt , nbed , n_bed_chars)) ; bed = 0.0

  !allocate bottom data
  allocate(bottom(0:mt , n_bot_vars)) ; bottom = 0.0

  !allocate discharge data
!  write(*,*)'allocating discharage arrays: ',numqbc
!  pause
!  allocate(seddis(numqbc,10)) ; seddis = 0.0 

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Alloc_Sed_Vars "

  End Subroutine Alloc_Sed_Vars


!==========================================================================
! Report Sediment Setup To Screen 
!==========================================================================

  Subroutine Report_Sed_Setup
  use control, only: msr 
  implicit none
  integer :: i

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Report_Sed_Setup "

  ! echo sediment model parameters to screen
  if(msr)then

  write(*,*)
  write(*,*)'------------------ Sediment Model Setup --------------------'
  if(oned_model)then
    write(*,*)'!  model                 : oneD   '
    write(*,*)'!  1D model: bedload     : deactivated '
  endif
  write(*,*)'!  nbed                  :',nbed    
  write(*,*)'!  nsed                  :',nsed    
  write(*,*)'!  n_report              :',n_report
  if(sed_start > 0)then
    write(*,*)'!  sed_start             :',sed_start
  else
    write(*,*)'!  sed_start             : no delay'  
  endif
  if(sed_hot_start)then
    write(*,*)'!  sed_hot_start         :',sed_hot_start
  endif
    
  write(*,*)'!  min_Srho              :',min_Srho
  if(sed_nudge)then
    write(*,*)'!  open boundary nudging :  active'
    write(*,*)'!  nudging relax factor  :',sed_alpha 
    write(*,*)'!  # nudging ramp its    :',sed_ramp  
  else
    write(*,*)'!  open boundary nudging :  not active'
  endif
  if(sed_source)then
    write(*,*)'!  point source forcing  :  active'
  else
    write(*,*)'!  point source forcing  :  not active'
  endif
  if(morpho_model)then
    write(*,*)'!  morphodynamics        :  active'
    write(*,*)'!  morpho_incr           :',morpho_incr   
    write(*,*)'!  morpho_strt           :',morpho_strt   
  else
    write(*,*)'!  morphodynamics        :  not active'
  endif
  write(*,*)'!  morpho_factor         :',morpho_factor
  if(bedload)then
    write(*,*)'!  bedload dynamics      :  active'
    write(*,*)'!  Critical Shields      :',Shield_Cr_MPM
    write(*,*)'!  Bedload Exponent      :',Gamma_MPM     
    write(*,*)'!  Bedload Constant      :',k_MPM     
    write(*,*)'!  Bedload Rate Coeff.   :',bedload_rate
    if(bedload_smooth)then
      write(*,*)'!  Bedload Smoother      :  active'
    endif 
  else
    write(*,*)'!  bedload dynamics      :  not active'
  endif
  if(susload)then
    write(*,*)'!  susload dynamics      :  active'
  else
    write(*,*)'!  susload dynamics      :  not active'
  endif
  if(inf_bed)then
    write(*,*)'!  sediment supply       :  infinite'
  else
    write(*,*)'!  sediment supply       :  finite'
  endif
  write(*,*)'!  Settling CFL          :',settle_cfl
  write(*,*)'!'
  write(*,*)'!   class        type   Sd50    Wset   tau_ce   '//&
            'tau_cd  Erate   Spor    Srho'
  do i=1,nsed
    write(*,'(1X,A1,1X,A12,1X,A6,6(F7.4,1X),F8.2)')'!', &
              trim(sed(i)%sname),trim(sed(i)%stype(1:6)), &
              sed(i)%Sd50,sed(i)%Wset, &
              sed(i)%tau_ce,sed(i)%tau_cd,sed(i)%erate,sed(i)%Spor, &
              sed(i)%Srho  
  end do

  end if

  ! echo sediment model initial conditions to screen
  call sed_report  

  write(*,*)'------------------------------------------------------------'
  write(*,*)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Report_Sed_Setup "

  End Subroutine Report_Sed_Setup

!==========================================================================
! Monitor Sediment Model: Compute and Write Statistics to Screen
!==========================================================================

  Subroutine Sed_Report  
  use control, only: msr ,par, mpi_fvcom_group
  use lims   , only: m,mt,kbm1, msrid
  implicit none
  real(dp) :: tau_min,tau_max,tmp1(3),tmp2(3),dthck_max
  integer :: i,dim1,dim2,ierr
#  if defined (FLUID_MUD)
  real(dp) :: temp1(5),temp2(5),mud_max,ua_max,va_max,ua_min,va_min
  REAL(DP), DIMENSION(3) :: SBUF,RBUF1,RBUF2,RBUF3
  real(dp) :: ESTOT,NSTOT
#  endif

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Sed_Report "

  !set up limits
  dim1 = mt 
  dim2 = kbm1
  
  !calculate max/min/rms concentrations
  do i=1,nsed
    sed(i)%cmax = maxval(sed(i)%conc(1:dim1,1:dim2))
    sed(i)%cmin = minval(sed(i)%conc(1:dim1,1:dim2))
    sed(i)%crms = sum(dble(abs(sed(i)%conc(1:dim1,1:dim2))))
    sed(i)%crms = sed(i)%crms/(dim1*dim2)
  end do

  dthck_max = maxval(abs(bottom(1:m,dthck)))
  tau_max   = maxval(taub(1:m))
  tau_min   = minval(taub(1:m))

#  if defined (FLUID_MUD)
  ESTOT = DBLE(NGL)
  NSTOT = DBLE(MGL)
  !calculate/max/min/rms mud layer
  mud_max = maxval(abs(D_FMLcm))
  ua_max = maxval(ua_fml)
  ua_min = minval(ua_fml)
  va_max = maxval(va_fml)
  va_min = minval(va_fml)

  SBUF(1)  = SUM(DBLE(UA_fml(1:n)))
  SBUF(2)  = SUM(DBLE(VA_fml(1:n)))
  SBUF(3)  = SUM(DBLE(D_FML(1:m)))

  RBUF1 = SBUF

#  endif
 
# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
  do i=1,nsed
    tmp1(1:3) = (/sed(i)%cmax,sed(i)%cmin,sed(i)%crms/)
    CALL MPI_REDUCE(tmp1,tmp2,3,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
    if(msr)then
      sed(i)%cmax = tmp2(1)
      sed(i)%cmin = tmp2(2)
      sed(i)%crms = tmp2(3)
    endif
  end do

  tmp1(1:3) = (/dthck_max,tau_min,tau_max/)
  CALL MPI_REDUCE(tmp1(1),tmp2(1),1,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
  CALL MPI_REDUCE(tmp1(2),tmp2(2),1,MPI_DP,MPI_MIN,MSRID-1,MPI_FVCOM_GROUP,IERR)
  CALL MPI_REDUCE(tmp1(3),tmp2(3),1,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
  dthck_max = tmp2(1) 
  tau_min   = tmp2(2) 
  tau_max   = tmp2(3) 

#  if defined (FLUID_MUD)
  temp1(1:5) = (/mud_max,ua_max,ua_min,va_max,va_min/)
  CALL MPI_REDUCE(temp1(1),temp2(1),1,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
  CALL MPI_REDUCE(temp1(2),temp2(2),1,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
  CALL MPI_REDUCE(temp1(3),temp2(3),1,MPI_DP,MPI_MIN,MSRID-1,MPI_FVCOM_GROUP,IERR)
  CALL MPI_REDUCE(temp1(4),temp2(4),1,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
  CALL MPI_REDUCE(temp1(5),temp2(5),1,MPI_DP,MPI_MIN,MSRID-1,MPI_FVCOM_GROUP,IERR)
  mud_max = temp2(1)
  ua_max  = temp2(2)
  ua_min  = temp2(3)
  va_max  = temp2(4)
  va_min  = temp2(5)
  CALL MPI_REDUCE(SBUF,RBUF1,3,MPI_DP,MPI_SUM,MSRID-1,MPI_FVCOM_GROUP,IERR)
  IF(ISNAN(RBUF1(1)).or.ISNAN(RBUF1(2)).or.ISNAN(RBUF1(3)))THEN
     CALL FATAL_ERROR("fluid mud caclulation fails with unstable iteration","Please check the model")
   END IF
#  endif

  ENDIF
# endif

  ! write statistics to screen
  if(.not.msr)return
  write(*,*  )'====================Sediment Model Stats======================'
  write(*,*  )'!  quantity              :     avg(g/L)      max(g/L)      '
  do i=1,nsed
    write(*,100)'!',trim(sed(i)%sname),'    :',sed(i)%crms,sed(i)%cmax
  end do
  write(*,*  )'!  max sed thick change  :     ',dthck_max
  write(*,*  )'!  max bottom stress     :     ',tau_max
  write(*,*  )'!  min bottom stress     :     ',tau_min
#  if defined (FLUID_MUD)
  write(*,*  )'!  max mud thickness     :     ',mud_max
  write(*,*  )'!  max u-speed of mud    :     ',ua_max
  write(*,*  )'!  min u-speed of mud    :     ',ua_min
  write(*,*  )'!  max v-speed of mud    :     ',va_max
  write(*,*  )'!  min v-speed of mud    :     ',va_min
  write(*,*  )'!  avg u-speed of mud    :     ',RBUF1(1)/ESTOT
  write(*,*  )'!  avg v-speed of mud    :     ',RBUF1(2)/ESTOT
  write(*,*  )'!  avg thickness of mud  :     ',RBUF1(3)/NSTOT
#  endif  

  100 format(1x,a1,a20,a5,3f12.6)
  101 format(1x,a26,f12.6)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Sed_Report "

  End Subroutine Sed_Report  

!==========================================================================
! Calculate Bed Load Transport using Meyer-Peter Muller Formulation   
!   Modify for influence of bathy slope
!   See eqns 24-27+35 + Warner et al, Comp. Geos.  2008 
!==========================================================================

  Real(sp) Function Bedload_MPM(rho_sed,D_50,tau_b) 
  use control, only: grav
  implicit none
  real(sp), intent(in) :: rho_sed
  real(sp), intent(in) :: D_50 
  real(sp), intent(in) :: tau_b  
  !--------------------------------
  real(sp) :: tau_nd,bedld_nd,R

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Bedload_MPM "

  R = (rho_sed/1025.-1)

  !nondimensional bottom stress 
  tau_nd = tau_b/(R*grav*D_50) 

  !compute bedload flux (PHI, eq. 25)
  bedld_nd = k_MPM*(MAX((tau_nd - Shield_Cr_MPM),0.0)**Gamma_MPM)

  !compute bedload (kg/(m-s)) (q_bl eq. 24)
  bedload_MPM = sqrt(R*grav*D_50)*D_50*rho_sed*bedld_nd

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Bedload_MPM "

  End Function Bedload_MPM

!==========================================================================
! Given Bed Load Transports, Calculate the flux using divergence 
!   - calculate bedload transport in elements
!   - modify transport using bed slope
!   - exchange at interprocessor boundaries
!   - calculate divergence to give flux out of / into bed
!==========================================================================
  Subroutine  Calc_Bedload_Flux() 

  !use all_vars, only : niec,ntrg,ncv,nt,dltxe,dltye,nprocs,n,msr,par,wubot,wvbot,myid,art1,m
  use all_vars
  use bcs, only : iobcn,i_obc_n
  use mod_obcs, only: next_obc
# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif
  integer  :: i1,ia,ib,ns,i,n1,n2,n3
  real(sp) :: qbl,ootaub_e,taub_e,flux
  real(sp) :: dbdx,dbdy,beta,slopex,slopey,sed_slope,max_slope,bld_dr
  real(sp), parameter :: sed_angle = 33.0  !sed "friction angle" in degrees
  real(sp), allocatable, target :: btx(:) 
  real(sp), allocatable, target :: bty(:) 
  real(sp), parameter :: taub_min = tiny(1.0_sp)

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Bedload_Flux "

  !initialize working bedload transport arrays
  allocate(btx(0:nt))
  allocate(bty(0:nt))

  !main loop over sed types (cohesive?)
  do ns = 1,nsed
    sed(ns)%bflx = 0.0
    btx = 0.0
    bty = 0.0

    !loop over elements, calculate bedload transport in x (btx) and  y (bty) at each
    ! => btx (qblx, eq 28) [after loop is kg/m]
    ! => bty (qbly, eq 29) [after loop is kg/m]
    do i=1,n
  
      !magnitude (units stress/density) of bottom shear stress at elements
      taub_e = sqrt(wubot(i)**2 + wvbot(i)**2)
    
      !calculate transport using Meyer-Peter Mueller formulation
      !qbl has units of kg-m^2/s
      qbl = bedload_mpm(sed(ns)%Srho,sed(ns)%Sd50,taub_e)
      if(taub_e > taub_min)then
        ootaub_e = 1./taub_e
      else
        ootaub_e = 0.0
      endif


      !calculate the bedload transport (kg/m)
      !note the negative sign: bottom stress in FVCOM (wubot/wvbot) is stress
      !on the fluid.  So a + current will produce a - bed stress
      !we want this transport to be in the direction of the current, thus the 
      !change of sign
      btx(i) = -DTsed*qbl*wubot(i)*ootaub_e  
      bty(i) = -DTsed*qbl*wvbot(i)*ootaub_e 
  
    end do
  
    !add bed slope contribution to our horizontal fluxes
    !Eq 35, Warner etal, Comp & Geo, 2008 => qbl_slope = slopex,slopey 
    !minimize the maximum downslope, but do not truncate upslope
    !upslope is limited by Lesser formula
    !range of slopefactor (slopex)  = [.55,5.64] for ..
    !..  upslope(dhdx<0)=-infty -> downslope(dhdx>0) = infty
    !dbdx > 0 is a downslope (see bld_dr to enforce this)
    sed_slope = tan(sed_angle*3.14159/180.)
    max_slope = 0.8*sed_slope
    do i=1,n
      n1 = nv(i,1) ; n2 = nv(i,2) ; n3 = nv(i,3)

      bld_dr = sign(1.0_sp,btx(i))
      dbdx   = min( (awx(i,1)*h(n1)+awx(i,2)*h(n2)+awx(i,3)*h(n3))*bld_dr ,max_slope)
      beta   = atan(dbdx)
      slopex = sed_slope/( (sed_slope-dbdx)*cos(beta)) 
      btx(i) = btx(i)*slopex

      bld_dr = sign(1.0_sp,bty(i))
      dbdy   = min( (awy(i,1)*h(n1)+awy(i,2)*h(n2)+awy(i,3)*h(n3))*bld_dr ,max_slope)
      beta   = atan(dbdy)
      slopey = sed_slope/( (sed_slope-dbdy)*cos(beta)) 
      bty(i) = bty(i)*slopey
      
    end do
  
    ! multiply by morphological time scale (to speed up time)
!    if(morpho_model)then
       btx = btx*morpho_factor
       bty = bty*morpho_factor
!    endif
  
    ! multiply by bedload transport coefficient (fudge factor) 
     btx = btx*bedload_rate 
     bty = bty*bedload_rate 
  
    ! limit transport to available bed mass
    ! ROMS does this, but doesn't make sense, mass is linked to divergence of  
    ! transport, not horizontal fluxes , check with J. Warner
     
  
    !exchange btx/bty across interprocessor boundaries 
# if defined (MULTIPROCESSOR)
    if(par)then
      call aexchange(ec,myid,nprocs,btx,bty)
    endif
# endif
  
    !calculate the divergence to give flux (kg) on the nodes
    do i=1,ncv
       i1   = ntrg(i)
       ia   = niec(i,1)
       ib   = niec(i,2)
       flux = (-btx(i1)*dltye(i) + bty(i1)*dltxe(i))
       sed(ns)%bflx(ia) = sed(ns)%bflx(ia)-flux
       sed(ns)%bflx(ib) = sed(ns)%bflx(ib)+flux
     end do

     !update the final bedload flux (kg/m^2) on the nodes
     !this flux is + out of the bed
     !if flow dir is positive and bed_stress is increasing in x
     !this should produce a net positive (out of bed) flux
     do i=1,m
       sed(ns)%bflx(i) = sed(ns)%bflx(i)/art1(i)
     end do

     !set the flux on the open boundary
     do i=1,iobcn
       sed(ns)%bflx( i_obc_n(i) ) = sed(ns)%bflx(next_obc(i))
     end do

     !shutdown bedload if node is locally an inactive node
     do i=1,m
       sed(ns)%bflx(i) = sed(ns)%bflx(i)*bottom(i,bflag)
     end do

     !if bathy is rough - choose smooth_bedload, averaging filter
     if(bedload_smooth)then
       call n2e2d(sed(ns)%bflx,btx)
       call e2n2d(btx,sed(ns)%bflx)
     endif
       
  end do !loop over sediment types

  deallocate(btx)
  deallocate(bty)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Bedload_Flux "  

  End Subroutine Calc_Bedload_Flux

!==========================================================================
! Given Bed Load Transports, Calculate the flux using divergence 
!   - calculate bedload transport in elements
!   - modify transport using bed slope
!   - exchange at interprocessor boundaries
!   - calculate flux on elements, upwinding using sign of tau dot n 
!   - convert bedload flux to nodes for morphology
!==========================================================================
  Subroutine  Calc_Bedload_Flux2() 

  !use all_vars, only : niec,ntrg,ncv,nt,dltxe,dltye,nprocs,n,msr,par,wubot,wvbot,myid,art1,m
  use all_vars
  use bcs, only : iobcn,i_obc_n
  use mod_obcs, only: next_obc
# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif
  integer  :: i1,ia,ib,ns,i,n1,n2,n3,j1,j2,cnt,jj,cc
  real(sp) :: qbl,ootaub_e,taub_e,flux
  real(sp) :: dbdx,dbdy,beta,slopex,slopey,sed_slope,max_slope,bld_dr
  real(sp) :: tau_x,tau_y,flowdir,bldx,bldy
  real(sp), parameter :: one = 1.0_sp
  real(sp), parameter :: eps = .5    
  real(sp), parameter :: sed_angle = 33.0  !sed "friction angle" in degrees
  real(sp), allocatable, target :: btx(:) 
  real(sp), allocatable, target :: bty(:) 
  real(sp), allocatable, target :: bflxe(:) 
  real(sp), parameter :: taub_min = tiny(1.0_sp)

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Bedload_Flux2 "

  !initialize working bedload transport arrays
  allocate(btx(0:nt))
  allocate(bty(0:nt))
  allocate(bflxe(0:nt))

  !main loop over sed types (cohesive?)
  do ns = 1,nsed
    sed(ns)%bflx = 0.0
    bflxe = 0.0
    btx   = 0.0
    bty   = 0.0

    !loop over elements, calculate bedload transport in x (btx) and  y (bty) at each
    ! => btx (qblx, eq 28) [after loop is kg/m]
    ! => bty (qbly, eq 29) [after loop is kg/m]
    do i=1,n
  
      !magnitude (units stress/density) of bottom shear stress at elements
      taub_e = sqrt(wubot(i)**2 + wvbot(i)**2)
    
      !calculate transport using Meyer-Peter Mueller formulation
      qbl = bedload_mpm(sed(ns)%Srho,sed(ns)%Sd50,taub_e)
      if(taub_e > taub_min)then
        ootaub_e = 1./taub_e
      else
        ootaub_e = 0.0
      endif


      !calculate the bedload transport (kg/m)
      !note the negative sign: bottom stress in FVCOM (wubot/wvbot) is stress
      !on the fluid.  So a + current will produce a - bed stress
      !we want this transport to be in the direction of the current, thus the 
      !change of sign
      btx(i) = -DTsed*qbl*wubot(i)*ootaub_e  
      bty(i) = -DTsed*qbl*wvbot(i)*ootaub_e 
  
    end do
  
    !add bed slope contribution to our horizontal fluxes
    !Eq 35, Warner etal, Comp & Geo, 2008 => qbl_slope = slopex,slopey 
    !minimize the maximum downslope, but do not truncate upslope
    !upslope is limited by Lesser formula
    !range of slopefactor (slopex)  = [.55,5.64] for ..
    !..  upslope(dhdx<0)=-infty -> downslope(dhdx>0) = infty
    !dbdx > 0 is a downslope (see bld_dr to enforce this)
    sed_slope = tan(sed_angle*3.14159/180.)
    max_slope = 0.8*sed_slope
    do i=1,n
      n1 = nv(i,1) ; n2 = nv(i,2) ; n3 = nv(i,3)

      bld_dr = sign(1.0_sp,btx(i))
      dbdx   = min( (awx(i,1)*h(n1)+awx(i,2)*h(n2)+awx(i,3)*h(n3))*bld_dr ,max_slope)
      beta   = atan(dbdx)
      slopex = sed_slope/( (sed_slope-dbdx)*cos(beta)) 
      btx(i) = btx(i)*slopex

      bld_dr = sign(1.0_sp,bty(i))
      dbdy   = min( (awy(i,1)*h(n1)+awy(i,2)*h(n2)+awy(i,3)*h(n3))*bld_dr ,max_slope)
      beta   = atan(dbdy)
      slopey = sed_slope/( (sed_slope-dbdy)*cos(beta)) 
      bty(i) = bty(i)*slopey
      
    end do
  
    ! multiply by morphological time scale (to speed up time)
!    if(morpho_model)then
       btx = btx*morpho_factor
       bty = bty*morpho_factor
!    endif
  
    ! multiply by bedload transport coefficient (fudge factor) 
     btx = btx*bedload_rate 
     bty = bty*bedload_rate 
  
    ! limit transport to available bed mass
    ! ROMS does this, but doesn't make sense, mass is linked to divergence of  
    ! transport, not horizontal fluxes , check with J. Warner
     
  
    !exchange btx/bty across interprocessor boundaries 
# if defined (MULTIPROCESSOR)
     if(par)then
       call aexchange(ec,myid,nprocs,btx,bty)
     endif
# endif

     !set dummy values for btx/bty to check divergence calc
     !do i=1,n
     !  btx(i) = 0.0 !-xc(i)
     !  bty(i) = -yc(i) 
     !end do
  
     !calculate the flux on the elements using first order upwind formulation 
     do i=1,ne
       ia=iec(i,1)
       ib=iec(i,2)
       j1=ienode(i,1)
       j2=ienode(i,2)
       !bed stress in x and y averaged to an edge , used for upwind selection
       tau_x = 0.5*(wubot(ia)+wubot(ib))  
       tau_y = 0.5*(wvbot(ia)+wvbot(ib))  
       !dot the bed stress vector with length-weighted edge normal 
       !note bed stress opposes flow dir so we need neg sign
       !the vector (-dltyc,dltxc) points from ia to ib
       flowdir  = -(-tau_x*dltyc(i) + tau_y*dltxc(i)) ! -(tau dot dL)
       !upwind select the x-dir and y-dir bedload transport 
       bldx  = ((one-sign(one,flowdir))*btx(ib)+(one+sign(one,flowdir))*btx(ia))*0.5_sp 
       bldy  = ((one-sign(one,flowdir))*bty(ib)+(one+sign(one,flowdir))*bty(ia))*0.5_sp 
       !accum contributions to bedload transport divergence calculations at elements ia/ib
       !on either side of edge i
       bflxe(ia) = bflxe(ia) - bldx*dltyc(i) + bldy*dltxc(i)
       bflxe(ib) = bflxe(ib) + bldx*dltyc(i) - bldy*dltxc(i)
     end do

     !divide by element area to complete divergence => bedload in kg/m^2
     !if there was a divergence of bedload transport, this should be a positive number
     do i=1,n
       bflxe(i) = bflxe(i)/art(i)
     end do
      
     !exchange bflxe across interprocessor boundaries
#    if defined (MULTIPROCESSOR)
     if(par)then
       call aexchange(ec,myid,nprocs,bflxe)
     endif
#    endif

     !interpolate element-based blfx to the nodes
     call e2n2d(bflxe,sed(ns)%bflx)

     !smooth using explicit smoother
     do i=1,m
       do jj=1,ntsn(i)
         cc = nbsn(i,jj)
         sed(ns)%bflx(i) = sed(ns)%bflx(i) + eps*sed(ns)%bflx(cc)
       end do
       sed(ns)%bflx(i) = sed(ns)%bflx(i)/(1. + eps*ntsn(i))
     end do

     !set the flux on the open boundary
     do i=1,iobcn
       sed(ns)%bflx( i_obc_n(i) ) = sed(ns)%bflx(next_obc(i))
     end do

     !shutdown bedload if node is locally an inactive node
     do i=1,m
       sed(ns)%bflx(i) = sed(ns)%bflx(i)*bottom(i,bflag)
     end do

  end do !loop over sediment types

  !cleanup 
  deallocate(btx)
  deallocate(bty)
  deallocate(bflxe)
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Bedload_Flux2 "

  End Subroutine Calc_Bedload_Flux2



!==========================================================================
! Calculate Erosion Rate (kg/m^2) 
! Calculate erosion rate using the method of Ariathurai and Arulanandan (1978)
!    -Erosion rates of cohesive soils. J. Hydraulics Division, ASCE, 104(2),279-282.-
!    surface_mass_flux[i] = erate[i]*(1-porosity)*bfrac[i]*(tau_w/tau_crit[i] - 1)
!      where
!      i               = index of sediment type
!      bfrac[i]        = fraction of bed composed by sediment i
!      erate[i]        = bed erodibility constant
!      porosity        = volume of voids/total volume in the top layer
!      tau_w           = shear stress on the bed
!      tau_crit[i]     = critical shear stress of 
!==========================================================================

  Subroutine Calc_Erosion
  use all_vars
  implicit none
  integer  :: ised ,i
  real(sp) :: bed_por,bed_frac,tau_w,tau_ce,erate,dep,active,erosion_orig
  real(sp) :: dx,dy 
  integer  :: cnt,jj,cc
  real(sp),dimension(mt) :: unode,vnode

  ! ------ new: Karsten Lettmann ---
  real(sp) :: dz_bed, eflux_avail
  ! --------------------------------

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Erosion "
  ! Compute bottom velocity at node
   if(.not.oned_model)then
   CALL E2N2D(u(:,kbm1),unode)
   CALL E2N2D(v(:,kbm1),vnode)
   endif

  !Compute Erosive Flux (kg/m^2), limited by available material in active layer 
  do ised=1,nsed
    do i=1,m
# if defined (WET_DRY)
      if(iswetn(i)/=1)cycle  ! no calculation when dry node
# endif      
      if(SEDIMENT_PARAMETER_TYPE/=UNIFORM)then
        erate  = sed(ised)%rate(i)
        tau_ce = sed(ised)%t_ce(i) 
      else
        erate  = sed(ised)%erate
        tau_ce = sed(ised)%tau_ce 
      end if 
      !calc taub from upwind or downwind elements - beta
# if defined (WAVE_CURRENT_INTERACTION)
      if(MB_BBL_USE.or.SG_BBL_USE.or.SSW_BBL_USE)then
         taub(i) = 1025.0 * taucwmax(i)
      else
# endif
        taub(i) = 0.
        cnt = 0.
        do jj=1,ntve(i)
          cc = nbve(i,jj)
          !if(xc(cc) < vx(i))then !cell is upstream, only for eastward flow
          !if(xc(cc) > vx(i))then !cell is downtream, only for eastward flow
          !dx=xc(cc)-vx(i);dy=yc(cc)-vy(i);
          !if((dx*unode(i)-dy*vnode(i))>0.) then ! arbitrary direction
            cnt = cnt + 1
            taub(i) = taub(i) + sqrt(wubot(cc)**2 + wvbot(cc)**2) 
          !endif
        end do
        taub(i) = 1025.*taub(i)/cnt
# if defined (WAVE_CURRENT_INTERACTION)
      end if
# endif
     
! J. Ge for fluid mud layer
# if !defined (FLUID_MUD)
      bed_por  = bed(i,1,iporo)
      bed_frac = sed(ised)%frac(i,1)
      tau_w    = taub(i)
      dep      = sed(ised)%dflx(i)
      active   = (1.0-bed_por)*Bed_Frac*(tau_w/tau_ce-1.0)*bottom(i,iactv) + dep
      !erosion_orig  = MAX(DTsed*Erate*(1.0-bed_por)*Bed_Frac*(tau_w/tau_ce-1.0),0.0)
      !sed(ised)%eflx(i) = erosion_orig 
      sed(ised)%eflx(i) = DTsed*erosion(erate,bed_por,bed_frac,tau_w,tau_ce) !min(erosion,active) !geoff
      !sed(ised)%eflx(i) = DTsed*min(erosion(erate,bed_por,bed_frac,tau_w,tau_ce),active) !min(erosion,active) !geoff

      !shutdown erosion if node is inactive 
      sed(ised)%eflx(i) = sed(ised)%eflx(i)*bottom(i,bflag)

      ! ========= new: Karsten Lettmann, Okt. 2012 ==========
      ! limit the eflux by the thickness of the top layer_report
      ! how much is available for erosion
      dz_bed = bed(i,1,ithck)
      eflux_avail = bed_frac * sed(ised)%srho*(1.0-bed_por) * dz_bed + dep

      ! take only 99% of possible material
      sed(ised)%eflx(i) = min(sed(ised)%eflx(i) , eflux_avail)
      ! ===============  end new ================================

!      if(i==1)print*,'calc_erosion:',dep,bed(i,1,ithck),eflux_avail,sed(ised)%eflx(i)

! J. Ge for fluid mud layer
# endif

    end do
  end do

! J. Ge for fluid mud layer

# if !defined (FLUID_MUD)
  !Update Concentration in Bottom of Water Column
  do ised=1,nsed
    do i=1,m
# if defined (WET_DRY)
      if(iswetn(i)==1)then
# endif
        sed(ised)%conc(i,kbm1) = sed(ised)%conc(i,kbm1) + sed(ised)%eflx(i)/(d(i)*dz(i,kbm1))
# if defined (WET_DRY)
      else
        sed(ised)%conc(i,kbm1) = 0.0_sp
      end if
# endif
    end do
  end do
# endif
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Erosion "

  End Subroutine Calc_Erosion

!==========================================================================
! Calculate Erosion (kgm^-2s^-1) using user-defined formula 
!==========================================================================
  Real(sp) Function Erosion(erate,bed_por,bed_frac,tau_w,tau_ce)
  implicit none
  real(sp), intent(in) ::  erate 
  real(sp), intent(in) ::  bed_por 
  real(sp), intent(in) ::  bed_frac 
  real(sp), intent(in) ::  tau_w 
  real(sp), intent(in) ::  tau_ce 
!Jianzhong
  real(sp),parameter :: a1=-0.144,a2=0.904,a3=-0.823,a4=0.204  
  real(sp)::ratio
!JIanzhong

  !standard CSTM formulation 
!  Erosion = MAX(Erate*(1.0-bed_por)*Bed_Frac*(tau_w/tau_ce-1.0),0.0)

!Jianzhong
  !Prooijen-Winterwerp formulation
   ratio=tau_w/tau_ce
   if(ratio<0.52)then
     Erosion=0.0
   elseif(ratio>1.7)then
     Erosion = Erate*(1.0-bed_por)*Bed_Frac*(tau_w/tau_ce-1.0)    
   else
     Erosion = Erate*(1.0-bed_por)*Bed_Frac*(a1*ratio**3+a2*ratio**2&
          &+a3*ratio+a4)
   endif
!JIanzhong

  !Winterwerp
  !Erosion =  MAX(Erate*Bed_Frac*(tau_w/tau_ce-1.0),0.0)
  
  End Function Erosion
!==========================================================================
! Calculate Depositional Flux and Update Concentration from Settling
!==========================================================================
  Subroutine Calc_Deposition
  use all_vars
  implicit none
  integer  :: ised ,i,mcyc,ncyc,k
  real(sp) :: wset(kbm1) ,c(kbm1),dx(kbm1),flux,DTmax,eps,DTdep
  real(sp) :: sal(kbm1),g1(kbm1),ds50,je,tau_cd,cb
  integer  :: cnt,jj,cc
  real(sp),dimension(0:mt,1:kb) :: unode,vnode
  real(sp),dimension(0:mt)      :: cbcn
  real(sp) :: ud


  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Deposition "

  ! Compute bottom velocity at node
  if(.not.oned_model)then
  CALL E2N3D(u,unode)
  CALL E2N3D(v,vnode)
  endif

  eps = epsilon(eps)

  !Loop over Sediment Class       
  do ised=1,nsed
    do i=1,m
      if(bottom(i,bflag) > 0.0)then
  
      ! no calculation when dry nodes
# if defined (WET_DRY)
      if(iswetn(i)/=1)cycle  
# endif

      !initialize flux
      flux = 0.0

      !setup 1D concentration and grid arrays
      c(1:kbm1) = sed(ised)%conc(i,1:kbm1)
      dx(1:kbm1) = dz(i,1:kbm1)*d(i)
      
      !calculate the settling velocity if sediment is cohesive
      if(sed(ised)%stype == 'cohesive')then
!        call calc_wset(kbm1,c,wset,sed(ised)%Wset)   !Mahta's method
        ds50 = sed(ised)%sd50
        sal(1:kbm1)=s1(i,1:kbm1)
      !calc taub from upwind or downwind elements - beta
# if defined (WAVE_CURRENT_INTERACTION)
      if(MB_BBL_USE.or.SG_BBL_USE.or.SSW_BBL_USE)then
         taub(i) = 1025.0 * taucwmax(i)
      else
# endif
        cnt = 0
        je  = 0.
        taub(i) = 0.
        do jj=1,ntve(i)
           cc=nbve(i,jj)
          !if(xc(cc) < vx(i))then !cell is upstream
          !if(xc(cc) > vx(i))then !cell is downtream
            cnt = cnt + 1
            taub(i) = taub(i) + sqrt(wubot(cc)**2 + wvbot(cc)**2) 
          ! endif
        end do
        taub(i) = 1025.*taub(i)/cnt
# if defined (WAVE_CURRENT_INTERACTION)
        end if
# endif
        g1 = 0.0
        do k=1,kbm1
          if(d(i)<=5.0)then
            je = mf**2*(unode(i,k)**2 + vnode(i,k)**2)/5**(4./3.) 
          else
            je = mf**2*(unode(i,k)**2 + vnode(i,k)**2)/d(i)**(4./3.) 
          end if          
          g1(k)=sqrt(grav*sqrt(unode(i,k)**2+vnode(i,k)**2)*je/Vcons)
        end do
        call calc_wset_floc(kbm1,c,sal,g1,ds50,sed(ised)%w0,wset)
!        if(d(i)<1.0)then
!          ud = 0.812*(ds50*1000.0)**(0.4)*(wset(kbm1)*1000.0*1.0)**(0.2)
!        else
!          ud = 0.812*(ds50*1000.0)**(0.4)*(wset(kbm1)*1000.0*d(i))**(0.2)
!        end if
!        tau_cd = 1025*cbcn(i)*ud**2
        tau_cd = sed(ised)%t_cd(i)
      else
        wset = sed(ised)%Wset
        if(SEDIMENT_PARAMETER_TYPE/=UNIFORM)then
           tau_cd = sed(ised)%t_cd(i)
        else
           tau_cd = sed(ised)%tau_cd
        end if
      endif

      !set up cycles (use max(CFL) == 1)
      DTmax = Settle_CFL*minval(dx(1:kbm1)/wset(1:kbm1))
      mcyc = int(DTsed/DTmax + 1. - eps)
      DTdep = DTsed/float(mcyc)

      !call flux-limited settling equation
      do ncyc = 1,mcyc
        if(sed(ised)%stype == 'cohesive' )then
          call settle_flux_floc(kbm1,c,dx,wset,flux,DTdep,tau_cd,taub(i))
        else
          call settle_flux(kbm1,c,dx,wset,flux,DTdep)
        end if
      end do

! Jianzhong Ge 03/05/2013
      sed(ised)%wset0(i,1:kbm1)=wset(1:kbm1)
! Jianzhong Ge 03/05/2013

      !store solution in 3D array
      sed(ised)%conc(i,1:kbm1) = c(1:kbm1)

      !store depositional flux
      sed(ised)%dflx(i) = flux 

      !shutdown erosion if node is inactive 
      sed(ised)%dflx(i) = sed(ised)%dflx(i)
    else
      sed(ised)%dflx(i) = 0.0
    endif
    end do
  end do

! J. Ge for fluid mud layer
# if defined (FLUID_MUD)
  Settling = 0.0_SP
  do ised=1,nsed
    do i=1,m
# if defined (WET_DRY)
      if(iswetn(i)==1)then
# endif
        Settling(i) = Settling(i) + sed(ised)%dflx(i)
# if defined (WET_DRY)
      else
        Settling(i) = 0.0_sp
      end if
# endif
    end do
  end do  
# endif

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Deposition "

  End Subroutine Calc_Deposition

!==========================================================================
! Calculate settling velocity of cohesive sediment using the method of
!           Qingyu Sha (including floccuation process)
! Refer : 1.Pingxing Ding et al, 2003, Numerical simulation of total
!         sediment under wave and current in Changjiang Estuary,
!         Acta Oceanologica Sinica, pp114-124,Vol.25, No.5 
!         2.Kelin Hu,2003,2-D Suspended Sediment numerical Simulation
!           under waves and currents in Yangtze Estuary,PHD doctoral
!           dissertation 2003 East China Normail University
!    or
!    User Defined
!==========================================================================
  Subroutine Calc_Wset_Floc(n,c,s,g1,ds50,w0,wset) 
  implicit none
  integer , intent(in   ) :: n
  real(sp), intent(in   ) :: ds50,w0, c(n),s(n),g1(n)
  real(sp), intent(inout) :: wset(n)
 
  integer :: i

  !setup for Changjiang Estuary
  real(sp) :: a = 0.274 ! proportional constant
  real(sp) :: b = 0.48  ! exponential constant of suspended sediment concertration
  real(sp) :: d = 0.03  ! exponential constant of salinity
  real(sp) :: e = 0.22  ! exponential constant of turblence intensity
  real(sp) :: f = 0.58  ! exponential constant of Ds50
  real(sp) :: factor

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Wset_Floc"

  do i=1,n
     factor = a*c(i)**b*s(i)**d*g1(i)**e/((ds50*1000.0)**f)
     if(factor<1.0.or.s(i)<5.0.or.s(i)>15.0)factor = 1.0
     wset(i)=w0*factor
     if(wset(i)>0.5e-3)wset(i)=0.5e-3
  end do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Wset_Floc"

  End Subroutine Calc_Wset_Floc

  Subroutine Calc_Wset_hindered(n,c,w0,wset,cmax,scale) 
  implicit none
  integer , intent(in   ) :: n
  real(sp), intent(in   ) :: w0,cmax,scale,c(n)
  real(sp), intent(inout) :: wset(n)
 
  integer :: i

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Wset_Floc2"

  do i=1,n
    if(c(i) < cmax)then
       wset(i) = w0
    else
      wset(i)=w0*(1.0-c(i)/1650.0)**5
      !wset(i) = w0*scale
      !wset(i) = w0*(cmax/c(i))**3
    endif
  end do


  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Wset_Floc2"

  End Subroutine Calc_Wset_hindered

!Jianzhong Ge 03/05/2013
  Subroutine Hindered_Diffusion
  use lims, only: m
  use all_vars,only: kh,kbm1
  implicit none
  real(sp) :: sed_conc(1:m,1:kbm1)  !mass concentration
  real(sp) :: vol_conc(1:m,1:kbm1)  !volume concentration
  integer :: i,ii,k
  real,parameter :: cs0=0.65
  real :: phi

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Hindered_Diffusion"

  sed_conc = 0.0  
  do ii=1,nsed
     sed_conc=sed_conc+sed(ii)%cnew
  end do

  do i=1,m
    do k=1,kbm1
      if(sed_conc(i,k)>2.5)then
!        vol_conc(i,k)=sed_conc(i,k)/2650.0
!        phi=1+(vol_conc(i,k)/Cs0)**0.8-2.0*(vol_conc(i,k)/Cs0)**0.4
!        kh(i,k)=kh(i,k)*phi**4.0
         kh(i,k)=kh(i,k)*0.02
!        print*,kh(i,k),phi**4.0,vol_conc(i,k),sed_conc(i,k)
      end if
    end do
  end do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Hindered_Diffusion"

  End Subroutine Hindered_Diffusion
!Jianzhong Ge 03/05/2013

!==========================================================================
! Calculate settling velocity of cohesive sediment using the method of
!    Mehta, 1986 or
!    Hindered settlement strategy
!    or
!    User Defined
!==========================================================================
  Subroutine Calc_Wset(n,c,wset,Wset_Mean) 
  implicit none
  integer , intent(in   ) :: n
  real(sp), intent(in   ) ::  c(n)
  real(sp), intent(inout) :: wset(n)
  real(sp), intent(in   ) :: Wset_Mean
  integer :: i

  !setup for Ariake Bay
  real(sp) :: a = 0.5 !Mehta proportional constant
  real(sp) :: b = 2.0 !Mehta exponential constant
  real(sp) :: cmax = 1.0 !Concentration for hindered settling

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Wset"

  do i=1,n
    if(c(i) < cmax)then
       wset(i) = a*(c(i)**b)  
    else
       !hindered settling here
    endif
  end do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Wset"

  End Subroutine Calc_Wset




!==========================================================================
! Calculate settling velocity of cohesive sediment using the method of
!    Richardson & Zaki (1954)
!    Hindered settlement strategy
!    or
!    User Defined
!==========================================================================
  Subroutine Calc_Wset_delft(n,c,Wset_max,wset) 
  implicit none
  integer , intent(in   ) :: n
  real(sp), intent(in   ) :: c(n)
  real(sp), intent(inout) :: wset(n)
  real(sp), intent(in   ) :: Wset_Max
  integer :: i

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Wset_delft"

  do i=1,n
     wset(i)=(1-c(i)/1600.0)**5.0*Wset_max
  end do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Wset_delft"

  End Subroutine Calc_Wset_delft

!==========================================================================
! Calculate Settling Flux and Update Sediment Concentration in Water Column 
!   use the SLIP limiter
!   Jameson, A., Analysis and Design of Numerical Schemes for Gas Dynamics I
!      Artificial Diffusion, Upwind Biasing, Limiters, and their effects
!      on Accuracy and Multigrid Convergence, International Journal of 
!      Computational Fluid Dynamics, Vol 4, 171-218, 1995.
!
!      -second order accurate in smooth regions
!      -guaranteed local extremum diminishing
!      
!   conv:  convective flux
!   diss:  dissipative flux
!   
!==========================================================================
  Subroutine Settle_Flux(n,c,dx,wset,flux,deltat) 
  implicit none
  integer , intent(in   ) :: n
	  real(sp), intent(inout) ::  c(n)
	  real(sp), intent(in   ) :: dx(n)
	  real(sp), intent(in   ) :: wset(n)
	  real(sp), intent(inout) :: flux 
	  real(sp), intent(in   ) :: deltat 
	  real(sp) :: conv(n+1),diss(n+1)
	  real(sp) :: cin(-1:n+2)
	  real(sp) :: win(-1:n+2)
	  real(sp) :: dis4,wvel
  integer  :: i

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Settle_Flux "

  !transfer to working array
  cin(1:n) = c(1:n)
  win(1:n) = wset(1:n) 

  !surface bcs (no flux)
  cin(0)  =  -cin(1) 
  cin(-1) =  -cin(1)
  win(0)  =  -wset(1)
  win(-1) =  -wset(1)

  !bottom bcs (extrapolate)
  cin(n+1) = cin(n) 
  cin(n+2) = cin(n)
  win(n+1) = win(n)
  win(n+2) = win(n)

  !flux computation
  do i=1,n+1
    wvel    = .5*(win(i)+win(i-1))  !settle velocity at interface
    dis4    = wvel/2.
    conv(i) = wvel*(cin(i)+cin(i-1))/2. 
    diss(i) = dis4*(cin(i)-cin(i-1)-lim(cin(i+1)-cin(i),cin(i-1)-cin(i-2))) 
  end do

  !zero out surface flux
  conv(1) = 0. ; diss(1) = 0.

  !update
  do i=1,n
    c(i) = cin(i) + (deltat/dx(i))*(-conv(i+1)+conv(i) + diss(i+1)-diss(i)) 
  end do

  !set bottom flux (> 0 = deposition)
  flux = flux + deltat*(conv(n+1)-diss(n+1))

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Settle_Flux "

 End Subroutine Settle_Flux


!==========================================================================
! Calculate Settling Flux and Update Sediment Concentration in Water Column 
!           use cohesive sediment method with floccuation process
! Refer : 1.Pingxing Ding et al, 2003, Numerical simulation of total
!         sediment under wave and current in Changjiang Estuary,
!         Acta Oceanologica Sinica, pp114-124,Vol.25, No.5 
!         2.Kelin Hu,2003,2-D Suspended Sediment numerical Simulation
!           under waves and currents in Yangtze Estuary,PHD doctoral
!           dissertation 2003 East China Normail University      
!
!   conv:  convective flux
!   diss:  dissipative flux
!   
!==========================================================================
  Subroutine Settle_Flux_Floc(n,c,dx,wset,flux,deltat,tau_cd,taub) 
  implicit none
  integer , intent(in   ) :: n
	  real(sp), intent(inout) ::  c(n)
	  real(sp), intent(in   ) :: dx(n)
	  real(sp), intent(in   ) :: wset(n)
	  real(sp), intent(inout) :: flux 
	  real(sp), intent(in   ) :: deltat 
          real(sp), intent(in   ) :: tau_cd,taub
	  real(sp) :: conv(n+1),diss(n+1)
	  real(sp) :: cin(-1:n+2)
	  real(sp) :: win(-1:n+2)
          real(sp) :: wset2(n)
	  real(sp) :: dis4,wvel
 
  integer  :: i

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Settle_Flux_Floc "
  
  !apply the critical shear stress for deposition
  wset2 = wset
!  if(taub>=tau_cd)then
!    wset2(n)=0.0
!  end if

  !transfer to working array
  cin(1:n) = c(1:n)
  win(1:n) = wset2(1:n) 

  !surface bcs (no flux)
  cin(0)  =  -cin(1) 
  cin(-1) =  -cin(1)
  win(0)  =  -wset2(1)
  win(-1) =  -wset2(1)

  !bottom bcs (extrapolate)
  cin(n+1) = cin(n) 
  cin(n+2) = cin(n)
  win(n+1) = win(n)
  win(n+2) = win(n)
  !flux computation
  do i=1,n+1
    wvel    = .5*(win(i)+win(i-1))  !settle velocity at interface
    dis4    = wvel/2.
    conv(i) = wvel*(cin(i)+cin(i-1))/2. 
    diss(i) = dis4*(cin(i)-cin(i-1)-lim(cin(i+1)-cin(i),cin(i-1)-cin(i-2))) 
  end do
  
  !zero out surface flux
  conv(1) = 0. ; diss(1) = 0.
  
  !update
  do i=1,n
    c(i) = cin(i) + (deltat/dx(i))*(-conv(i+1)+conv(i) + diss(i+1)-diss(i)) 
  end do

  !set bottom flux (> 0 = deposition)
  flux = flux + deltat*(conv(n+1)-diss(n+1))

  if(flux<0)call fatal_error("Error: Deposition Flux is negative!!!")

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Settle_Flux_Floc "

 End Subroutine Settle_Flux_Floc





!==========================================================================
! Calculate LED Limiter L(u,v)  
!==========================================================================
  Function Lim(a,b)
  real(sp) lim,a,b
  real(sp) q,R
  real(sp) eps
  eps = epsilon(eps)
  
 ! exponent
 ! qlim = 0. !1st order
 ! qlim = 1. !minmod
 ! qlim = 2. !van leer

  R = abs(   (a-b)/(abs(a)+abs(b)+eps) )**qlim
  lim = .5*(1-R)*(a+b)

  End Function Lim

!==========================================================================
! Update Surface Bed Properties Due to Bedload/SusLoad Flux 
! Change:
!   sed(:)%mass(:,1)   to reflect increased/decrease mass (kg)
!   bed(:,1,itchk)     to reflect changed bed thickness properties
!==========================================================================

  Subroutine Update_Surface_Bed(fluxtype)
  use lims, only: m
  implicit none
  character(len=*) :: fluxtype
  real(sp) :: dz,accum,flux,t_mass,eps,junk
  integer  :: i,ised

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Update_Surface_Bed "

  eps = epsilon(eps)

  !======= adjust for bedload flux ===================
  if(fluxtype == 'bedload')then

  !top layer mass and thickness
  !bflx > 0 means material is lost from the bed.
  do i=1,m
    do ised=1,nsed
      flux = sed(ised)%bflx(i)

      !change top layer bed mass of each sed type 
      accum  = sed(ised)%mass(i,1) - flux 
      sed(ised)%mass(i,1) = max(accum,0.0) 

      !change top layer bed thickness 
      dz = flux/(sed(ised)%srho*(1.0-bed(i,1,iporo)))
      junk = bed(i,1,ithck)
      bed(i,1,ithck) = max(bed(i,1,ithck)-dz,0.0)
      bottom(i,morph) = bottom(i,morph) - dz 

    end do
  end do

  !top layer fractions 
  do i=1,m
    t_mass = 0.0_sp
    do ised=1,nsed
      t_mass = t_mass + sed(ised)%mass(i,1)  
    end do
    do ised=1,nsed
      sed(ised)%frac(i,1) = sed(ised)%mass(i,1)/MAX(t_mass,eps)
    end do
  end do

  !======= adjust for susload flux ===================
  elseif(fluxtype =='susload')then
  
  do i=1,m
    do ised=1,nsed
      flux = (sed(ised)%eflx(i)-sed(ised)%dflx(i))*morpho_factor
      !if depositional + first deposit time step, save mass and loop 
      !http://woodshole.er.usgs.gov/project-pages/sediment-transport/ 
      sed(ised)%depm(i) = 0.0
      if(flux < 0)then                                  !!depositional
        if(T_model > bed(i,1,iaged)+1.1*DTsed .and. &   !!first deposit time step
           bed(i,1,ithck) > thck_cr)then                !!thickness surpasses critical
          sed(ised)%depm(i) = -flux
          cycle
        else                                            !!update age of surface layer
          bed(i,1,iaged) = T_model
        endif
      endif

      !change top layer bed mass of each sed type 
      accum  = sed(ised)%mass(i,1) - flux 
      sed(ised)%mass(i,1) = max(accum,0.0) 

      !change top layer bed thickness 
      dz = flux/(sed(ised)%srho*(1.0-bed(i,1,iporo)))
      bed(i,1,ithck) = max(bed(i,1,ithck)-dz,0.0)
      bottom(i,morph) = bottom(i,morph) - dz 

    end do
  end do

  !top layer fractions 
  do i=1,m
    t_mass = 0.0_sp
    do ised=1,nsed
      t_mass = t_mass + sed(ised)%mass(i,1)  
    end do
    do ised=1,nsed
      sed(ised)%frac(i,1) = sed(ised)%mass(i,1)/MAX(t_mass,eps)
    end do
  end do

  !======= error ======================
  else
    write(*,*)'argument to Update_Surface_Bed must be: "susload" or "bedload"'
    stop
  endif !fluxtype

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Update_Surface_Bed "

  End Subroutine Update_Surface_Bed


!==========================================================================
! Update Bottom Properties due to changes in surface layer sed fractions 
! Use geometric mean, good for log normal distribution
!==========================================================================
  Subroutine Update_Bottom_Properties
  use lims, only: m
  implicit none
  integer  :: i,ised
  real(sp) :: temp,sum1,sum2,sum3,sum4,eps

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Update_Bottom_Properties "

  ! initialize
  eps = epsilon(temp) 
  
  ! update bottom properties using new bed fractions
  do i=1,m
    sum1 = 1.0 ; sum2 = 1.0 ; sum3 = 1.0 ; sum4 = 1.0
    do ised=1,nsed
      if(SEDIMENT_PARAMETER_TYPE/=UNIFORM)then
        sum1 = sum1*(sed(ised)%t_ce(i))**sed(ised)%frac(i,1) 
      else
        sum1 = sum1*(sed(ised)%tau_ce)**sed(ised)%frac(i,1)
      end if  
      sum2 = sum2*(sed(ised)%Sd50   )**sed(ised)%frac(i,1)
      sum3 = sum3*(sed(ised)%wset   )**sed(ised)%frac(i,1)
      sum4 = sum4*(sed(ised)%Srho   )**sed(ised)%frac(i,1)
    end do
    bottom(i,itauc) = sum1
    bottom(i,isd50) = sum2
    bottom(i,isphi) = -log(bottom(i,isd50)*1000.)/log(2.0)
    bottom(i,iwset) = sum3
    bottom(i,idens) = max(sum4,min_Srho) 
  end do 
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Update_Bottom_Properties "

  End Subroutine Update_Bottom_Properties 

!==========================================================================
! Update Active Layer Thickness 
!   Use method of Harris and Wiberg (1997):
!   -Approaches to quantifying long-term continental shelf sediment transport
!   -with an example from the northern California STRESS mid-shelf site.
!   [Continental Shelf REsearch, 17, 1389-1418]
!
!   z_a = max[k_1*(tau_bot - tau_ce)*rho_0 , 0] + k_2 * D_50
!       where
!       k_1 = .007 (empirical constant)
!       k_2 = 6.0  (empirical constant)
!       tau_bot = shear stress on the bed
!       tau_ce  = critical shear stress for erosion of the bed (averaged over
!                 sediment classes)
!       D_50    = median grain diameter of the surface sediment
!       note: our shear stresses are nondimensionalized by rho, so
!             rho does not appear in actual calculation (below)
!==========================================================================
  Subroutine Calc_Active_Layer 
  use lims,     only: m
  implicit none
  integer  :: i
  real(sp) :: mfac

  mfac = dim(morpho_factor,1.0) 
  !calculate active layer thickness 
  do i=1,m
    bottom(i,iactv) = mfac*max(0.,.007*(taub(i)-bottom(i,itauc)))+6.0*bottom(i,isd50)
  end do

  End Subroutine Calc_Active_Layer
  
!==========================================================================
! Update Stratigraphy       
!   -Calculate empirical active layer thickness from bottom stress
!   -Expand top layer to this thickness by shifting necessary mass up
!    through the bed.
!==========================================================================

  Subroutine Update_Stratigraphy
  use lims,     only: m
  implicit none
  integer :: i,k,ksed,ks,ised
  real(sp) :: thck_avail,thck_to_add,eps,tmp1,tmp2,tmp3,top_layer_mass
  real(sp) :: delta(m) 
 
  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Update_Stratigraphy "

  !initialize
  eps = epsilon(tmp1) 
  delta = 0.0

  !calculate deficit between active layer thickness and surface layer thickness
  do i=1,m
    delta(i) = bottom(i,iactv)-bed(i,1,ithck)
  end do
  
  ! (single layer system) ensure top layer > active layer thickness 
  if(nbed == 1)then 
    do i=1,m
       if(delta(i) > 0) bottom(i,iactv) = bed(i,1,ithck)
    end do
    return
  endif
   
  ! (multi layer system) ensure top layer > active layer thickness 
  do i=1,m
    thck_to_add = delta(i)
    if(thck_to_add > 0.0)then !must redistribute layers 
      thck_avail=0.0

      !find fractional layer (below which there will be no mass after forming active layer) 
      Ksed=1 
      do k=2,nbed
        if (thck_avail < thck_to_add) then
          thck_avail=thck_avail+bed(i,k,ithck)
          Ksed=k
        end if
      end do

      !if not enough material in bed to satisfy active layer req, use all available 
      if (thck_avail < thck_to_add) then
        bottom(i,iactv)=bed(i,1,ithck)+thck_avail
        thck_to_add=thck_avail
      end if
  
      !update the bed mass of surface and fractional layers
      tmp2=MAX(thck_avail-thck_to_add,0.0)/MAX(bed(i,Ksed,ithck),eps)
      do ised=1,nsed
        tmp1=0.0 
        do k=1,Ksed
          tmp1=tmp1+sed(ised)%mass(i,k)
        end do
        sed(ised)%mass(i,1   ) = tmp1-sed(ised)%mass(i,Ksed)*tmp2
        sed(ised)%mass(i,Ksed) = sed(ised)%mass(i,Ksed)*tmp2
      end do

      !update thickness of fractional layer Ksed
      bed(i,Ksed,ithck)=MAX(thck_avail-thck_to_add,0.0)

      !update top layer bed fraction 
      top_layer_mass = 0.0
      do ised=1,nsed
        top_layer_mass = top_layer_mass + sed(ised)%mass(i,1)
      end do
      do ised=1,nsed
        sed(ised)%frac(i,1)=sed(ised)%mass(i,1)/MAX(top_layer_mass,eps) 
      end do

      !update bed thickness of top layer
      bed(i,1,ithck)=bottom(i,iactv)

      !pull layers Ksed to Bottom up to fill layer 2 down 
      do k=Ksed,Nbed
        ks=Ksed-2
        bed(i,k-ks,ithck)=bed(i,k,ithck)
        bed(i,k-ks,iporo)=bed(i,k,iporo)
        bed(i,k-ks,iaged)=bed(i,k,iaged)
        do ised=1,nsed
          sed(ised)%frac(i,k-ks) = sed(ised)%frac(i,k)
          sed(ised)%mass(i,k-ks) = sed(ised)%mass(i,k)
        end do
      end do
 
      !split what was in the bottom layer to fill empty bottom cells
      !note:  porosity of nbed (bed(i,nbed,iporo)) does not change.
      ks=Ksed-2
      tmp3=1.0/real(ks+1)
      do k=Nbed,Nbed-ks,-1
        bed(i,k,ithck)=bed(i,Nbed-ks,ithck)*tmp3
        bed(i,k,iaged)=bed(i,Nbed-ks,iaged)
        do ised=1,nsed
          sed(ised)%frac(i,k)=sed(ised)%frac(i,Nbed-ks)
          sed(ised)%mass(i,k)=sed(ised)%mass(i,Nbed-ks)*tmp3
        end do
      end do

     end if !thck_to_add > 0
  end do !loop over nodes

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Update_Stratigraphy "
    
  End Subroutine Update_Stratigraphy

!==========================================================================
! Add New Layer                 
!   If first deposit time step:  
!       a.) combine bottom layers 
!       b.) shift upper layers down
!       c.) create new top layer
!==========================================================================

  Subroutine Add_New_Layer 
  use lims,     only: m
  implicit none
  integer :: i,k,ksed,ised
  real(sp) :: eps,t_mass,cnt
 
  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Add_New_Layer "

  !initialize
  eps = epsilon(eps) 

  do i=1,m
  
    cnt = 0.0
    do ised=1,nsed
      cnt = cnt + sed(ised)%depm(i)
    end do

    if(cnt > 0)then !!need new surface layer

      if(nbed > 1)then!--->
      !combine bottom two layers
      bed(i,nbed,ithck) =     bed(i,nbed-1,ithck) + bed(i,nbed,ithck)
      bed(i,nbed,iporo) = .5*(bed(i,nbed-1,iporo) + bed(i,nbed,iporo))
      bed(i,nbed,iaged) = .5*(bed(i,nbed-1,iaged) + bed(i,nbed,iaged))

      t_mass =0.0
      do ised=1,nsed
        sed(ised)%mass(i,nbed) = sed(ised)%mass(i,nbed-1)+sed(ised)%mass(i,nbed) 
        t_mass = t_mass +sed(ised)%mass(i,nbed)
      end do
      do ised=1,nsed
        sed(ised)%frac(i,nbed) = sed(ised)%mass(i,nbed)/MAX(t_mass,eps)
      end do

      !push layers down
      do k=Nbed-1,2,-1
        bed(i,k,ithck) = bed(i,k-1,ithck)
        bed(i,k,iporo) = bed(i,k-1,iporo)
        bed(i,k,iaged) = bed(i,k-1,iaged)
        do ised =1,nsed
          sed(ised)%frac(i,k) = sed(ised)%frac(i,k-1)
          sed(ised)%mass(i,k) = sed(ised)%mass(i,k-1)
        end do
      end do

      !refresh top layer in multi layer bed
      bed(i,1,ithck)=0.0
      do ised=1,nsed
        sed(ised)%mass(i,1)=0.0
      end do
      end if ! <---

     

      !update surface layer properties
      t_mass = 0.0
      do ised=1,nsed
        sed(ised)%mass(i,1) = sed(ised)%mass(i,1)+ sed(ised)%depm(i)
        t_mass = t_mass + sed(ised)%mass(i,1)
        bed(i,1,ithck)=bed(i,1,ithck)+sed(ised)%depm(i)/    &
                         (sed(ised)%Srho*(1.0-bed(i,1,iporo)))
      end do
      bed(i,1,iaged) = T_model  
      do ised=1,nsed
        sed(ised)%frac(i,1) = sed(ised)%mass(i,1)/MAX(t_mass,eps) 
        sed(ised)%depm(i  ) = 0.0
      end do
 
    end if  !need new surface layer
  end do !node loop

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Add_New_Layer "

  End Subroutine Add_New_Layer

!==========================================================================
! Update Total Sediment Thickness for Reporting
!==========================================================================

  Subroutine Update_Thickness_Delta
  use lims,     only: m
  implicit none
  integer :: i,k,ised
 
  !shift last thickness to lthck
  bottom(:,lthck) = bottom(:,nthck)
  bottom(:,nthck) = 0.0

  !calculate new thickness
  do i=1,m
    do k=1,nbed
      bottom(i,nthck) = bottom(i,nthck) + bed(i,k,ithck)
    end do
  end do

  !calculate new delta
  do i=1,m
    bottom(i,dthck) = bottom(i,dthck) + bottom(i,nthck) - bottom(i,lthck)
  end do

  !calculate new total mass
  do i=1,m
    bottom(i,tmass) = 0.0 
    do ised=1,nsed
      do k=1,nbed
        bottom(i,tmass) = bottom(i,tmass) + sed(ised)%mass(i,k)
      end do
    end do
  end do
  
  End Subroutine Update_Thickness_Delta

!==========================================================================
! Report Layer Stats for Debug   
!==========================================================================

  Subroutine Layer_Report(i)
  implicit none
  integer, intent(in) :: i
  integer :: j,ised

  write(*,*)'----------bottom props'
  write(*,*)'active thickness',bottom(i,iactv)
  write(*,*)'total  thickness',bottom(i,nthck)
  write(*,*)'mean grain size ',bottom(i,isd50)
  write(*,*)'mean density    ',bottom(i,idens)
  write(*,*)'mean settle V   ',bottom(i,iwset)
  write(*,*)'mean E stress   ',bottom(i,itauc)

  write(*,*)'--------- bed props'
  write(*,*)'total  thickness',bottom(i,nthck)
  write(*,*)'total  mass'     ,bottom(i,tmass)
  if(bottom(i,nthck) <= 0.)stop

  End Subroutine Layer_Report

!==========================================================================
!  Setup open boundary nudging for sediment obc
!   - check for existence of sed nudging file
!   - read data (concentration of each sed type at each obc)
!==========================================================================
  Subroutine Setup_Sed_OBC
  Use Control,  only : casename,input_dir,msr,serial,par
!  Use Mod_OBCs, only : iobcn,i_obc_gl,iobcn_gl
  Use Mod_Utils, only: pstop
# if defined (MULTIPROCESSOR)
  Use Mod_Par,  only : nlid
# endif

  integer :: i,j,tmp,ncnt,i1
  integer,  allocatable :: vtmp(:),atemp(:,:)
  real(sp), allocatable :: sed_obc_gl(:,:)
  integer, parameter    :: sedobc = 81
  character(len=120)    :: fname
  logical :: fexist
  

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Setup_Sed_OBC "

  !------------------------------------------------------
  ! make sure sediment nudging data file exists 
  !------------------------------------------------------
  fname = "./"//trim(input_dir)//"/"//trim(casename)//"_sed_nudge.dat"

  inquire(file=trim(fname),exist=fexist)
  if(msr .and. .not.fexist)then
    write(*,*)'sediment nudging data: ',fname,' does not exist'
    write(*,*)'halting.....'
    call pstop
  end if

  open(unit=sedobc,file=fname,form='formatted')
  !------------------------------------------------------
  ! allocate global data
  !------------------------------------------------------
  allocate(sed_obc_gl(Nobc_gl,nsed)) ; sed_obc_gl = 0.

  !header
  read(sedobc,*)
  read(sedobc,*)
  read(sedobc,*)

  !# obc nodes
  read(sedobc,*) tmp 

  if(tmp /= Nobc_gl .and. msr)then
    write(*,*)'number of open boundary nodes for sediment nudging: ',tmp
    write(*,*)'this is not consistent with number of open boundary nodes in'
    write(*,*)'the model: ',Nobc_gl
    write(*,*)'halting'
    call pstop
  endif

  allocate(vtmp(Nobc_gl))
  do i=1,Nobc_gl  
    read(sedobc,*) i1,vtmp(i),(sed_obc_gl(i,j),j=1,nsed)
  end do
  close(sedobc)

  ! !------------------------------------------------------
  ! ! make sure global obc node list is consistent
  ! !------------------------------------------------------
  ! 
  ! if(msr)then
  ! do i=1,iobcn_gl
  !   if(i_obc_gl(i) /= vtmp(i))then
  !     write(*,*)'obc node list for sediment nudging not consistent with model'
  !     write(*,*)'obc node list' 
  !     write(*,*)'halting'
  !     call pstop
  !   end if
  ! end do
  ! endif
    
  !---------------------------------------------------------
  ! shift sediment nudging data to local open boundary nodes
  !---------------------------------------------------------

  if(serial)then
    do i=1,Nobc_gl
      do j=1,nsed
        sed(j)%cobc(i) = sed_obc_gl(i,j)
      end do
    end do
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
    allocate(atemp(Nobc_gl,nsed))
    ncnt = 0
    !!set up local open boundary nodes
    do i=1,Nobc_gl
      i1 = nlid( vtmp(i) )
      if(i1 /= 0)then
        ncnt = ncnt + 1
        do j=1,nsed
          atemp(ncnt,j) = sed_obc_gl(i,j)
        end do
      end if
    end do

    !transfer sed nudging data global --> local
    if(ncnt > 0)then
      do j=1,nsed
        sed(j)%cobc(1:ncnt) = atemp(1:ncnt,j)
      end do
    end if

    deallocate(atemp) 
  end if
# endif

  if(msr)write(*,*)'obc setup for sediment model: complete'

  !---------------------------------------------------------
  ! cleanup 
  !---------------------------------------------------------
  deallocate(sed_obc_gl)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Setup_Sed_OBC "

  End Subroutine Setup_Sed_OBC



   SUBROUTINE UPDATE_DEPTH
!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: TEMP, DX, DY
   INTEGER  :: I,I1,K,J1,J2,CC
   LOGICAL  :: ISONB2
!------------------------------------------------------------------------------|
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: UPDATE_DEPTH"

# if defined (MULTIPROCESSOR)
    !IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,bottom )
    !IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,H )
# endif

   DO I = 1,M
!
!   DEPTHS ON OPEN BOUNADY NODES AND THE NODES NEXT TO THOSE NODES ARE KEPT INVARIANT
!
    ISONB2=.FALSE.
     DO I1 = 1,IOBCN
       J1 = I_OBC_N(I1)
       J2 = NEXT_OBC(I1)
       IF(I==J1 .OR. I==J2) ISONB2=.TRUE.
     ENDDO
    IF(.NOT. ISONB2) H(I) = H(I) - bottom(I,morph)
    bottom(I,morph) = 0.0
   ENDDO
# if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,H )
# endif
   CALL N2E2D(H,H1)
   !DO I = 1, N
   ! CC = 0
   ! H1(I) = 0.
   ! DO K = 1,3
   !  DX = VX(NV(I,K)) - XC(I)
   !  DY = VY(NV(I,K)) - YC(I)
   !  IF( (DX*U(I,KBM1) - DY*V(I,KBM1)) <  0.) THEN !node is upstream
   !   CC = CC + 1
   !   H1(I) = H1(I) + H(NV(I,K))
   !  ENDIF
   ! ENDDO 
   !  H1(I) = H1(I)/CC
   !ENDDO
# if defined (MULTIPROCESSOR)
    !IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,H1 )
    !IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,EL1,ET1)
    !IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,EL ,ET )
# endif
   D   = H + EL
   DT  = H + ET
   D1  = H1+EL1
   DT1 = H1+ET1
# if defined (MULTIPROCESSOR)
    !IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,D1,DT1 )
    !IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,D ,DT )
# endif


   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: UPDATE_DEPTH"

   END SUBROUTINE UPDATE_DEPTH

# endif
   

End Module Mod_Sed 
