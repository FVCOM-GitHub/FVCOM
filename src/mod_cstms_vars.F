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
Module Mod_CSTMS_vars

#if defined (SEDIMENT) && (CSTMS_SED)
Use Mod_Par
Use Mod_Prec 
Use Mod_Types

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
! morph_fac    => morphological factor
!-----------------------------------------------------------------------
!       Horizontal and vertical constant mixing coefficients.
!
! Akt_bak      => Background vertical mixing coefficient (m2/s) for
!                 tracers.                  
! nl_tnu2      => NLM lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
! nl_tnu4      => NLM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for tracers.
!-----------------------------------------------------------------------
! LtracerSponge=> Switches to increase/decrease horizontal viscosity
!                 and/or diffusion  in specific areas of the
!                 application domain (like sponge areas).
! LtracerSrc   => Switches to activate point Source/Sinks in an
!                 application: Tracer transport
! LtracerCLM   => Switches to read and process climatology fields.
! LnudgeTCLM   => Switched to nudge to climatology fields.
! Tnudg        => Inverse time-scales (1/s) for nudging at open
!                 boundaries and sponge areas.
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
! Jianzhong Ge 03/05/2013
! wset         => store settling velocity
! Jianzhong Ge 03/05/2013
! arraysize    => spatial size of sediment arrays
!--------------------------------------------------
type sed_type
  character(len=20) :: sname    
  character(len=20) :: sname2    
  character(len=20) :: stype 
  real(sp)          :: csed_initial ! initial concentration
  real(sp)          :: Erate        ! erosion rate
  real(sp)          :: Sd50         ! mediam grain diameter
  real(sp)          :: Srho         ! grain density
  real(sp)          :: Wsed         ! settling velocity
  real(sp)          :: poros        ! porosity
  real(sp)          :: tau_ce       ! shear for erosion
  real(sp)          :: tau_cd       ! shear for deposition
  real(sp)          :: morph_fac    ! morphological factor
  real(sp)          :: Akt_bak        ! m2/s
  real(sp)          :: nl_tnu2        ! m2/s
  real(sp)          :: nl_tnu4        ! m2 s-1/2
  real(sp)          :: Tnudg          ! Tracers
  logical           :: LtracerSponge   ! diffusion
  logical           :: LtracerSrc      ! tracers
  logical           :: LtracerCLM      ! tracers
  logical           :: LnudgeTCLM      ! tracers
  real(sp)          :: Chin       ! concentration for hindered settling 
  real(sp)          :: Wrdc       ! reduction scale for hindered settling velocity
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
  real(sp), allocatable :: wset0(:,:) ! vertical settling velocity
  real(sp), allocatable :: t_cd(:)    ! =tau_cd
  real(sp), allocatable :: t_ce(:)    ! =tau_ce
  real(sp), allocatable :: rate(:)    ! =erate
  integer           :: arraysize
end type sed_type






TYPE TYPE_SEDBED
  real(sp), allocatable :: avgbedldu(:,:)
  real(sp), allocatable :: avgbedldv(:,:)

  real(sp), allocatable :: bed(:,:,:)
  real(sp), allocatable :: bed_frac(:,:,:)
  real(sp), allocatable :: bed_mass(:,:,:,:)

  real(sp), allocatable :: bed_thick0(:)
  real(sp), allocatable :: bed_thick(:,:)

  real(sp), allocatable :: bedldu(:,:)
  real(sp), allocatable :: bedldv(:,:)

  real(sp), allocatable :: bottom(:,:)

  real(sp), allocatable :: ero_flux(:,:)
  real(sp), allocatable :: settling_flux(:,:)

  real(sp), allocatable :: Dstp_max(:,:)

END TYPE TYPE_SEDBED


TYPE TYPE_SEDFLOCS
  real(sp), allocatable :: f_diam(:)
  real(sp), allocatable :: f_vol(:)
  real(sp), allocatable :: f_rho(:)
  real(sp), allocatable :: f_cv(:)
  real(sp), allocatable :: f_l3(:)
  real(sp), allocatable :: f_mass(:)
  real(sp), allocatable :: f_coll_prob_sh(:,:)
  real(sp), allocatable :: f_coll_prob_ds(:,:)
  real(sp), allocatable :: f_l1_sh(:,:)
  real(sp), allocatable :: f_l1_ds(:,:)
  real(sp), allocatable :: f_g3(:,:)
  real(sp), allocatable :: f_l4(:,:)
  real(sp), allocatable :: f_g1_sh(:,:,:)
  real(sp), allocatable :: f_g1_ds(:,:,:)
  real(sp), allocatable :: f_g4(:,:,:)
END TYPE TYPE_SEDFLOCS



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


!
!--------------------------------------------------
! Physical Dynamics Control Parameters and Logicals
!--------------------------------------------------
integer :: ncs
integer :: nns
integer :: nst
logical :: oned_model
logical :: susload
logical :: bedload
logical :: sed_morph
logical :: sed_flocs
logical :: sed_biomass
logical :: cohesive_bed
logical :: mixed_bed
logical :: noncohesive_bed2
logical :: sed_biodiff
logical :: sed_defloc

!--------------------------------------------------
! seagrass biomass options
!--------------------------------------------------
logical :: seagrass_bottom
logical :: seagrass_sink

!--------------------------------------------------
! bedload transport algortihm options
!--------------------------------------------------
logical :: BEDLOAD_MPM
logical :: BEDLOAD_SOULSBY
logical :: SLOPE_NEMETH
logical :: SLOPE_LESSER

!--------------------------------------------------
! vertical settling options
!--------------------------------------------------
logical :: LINEAR_CONTINUATION
logical :: NEUMANN

!--------------------------------------------------
! bottom condition for critical deposition shear stress
!--------------------------------------------------
logical :: SED_TAU_CD_CONST
logical :: SED_TAU_CD_LIN

!--------------------------------------------------
! bed biodiffusitivity options
!--------------------------------------------------
logical :: DB_PROFILE

!--------------------------------------------------
!Initial condition parameters
!  Read in sediment input file
!  Can be overwritten in init_sed.F or through 
!  restart
!
!  init_bed_porosity:  initial bed porosity [0,1]
!                      [1:NBED] 
!  init_bed_thickness: initial bed thickness [m] [> 0]
!                      [1:NBED] 
!  init_bed_fraction:  initial bed Sediment fraction of each size
!                      class in each bed layer(nondimensional: 0-1.0)
!                      Sum of bed_frac = 1.0. 
!                      [1:NST,1:NBED] 
!  init_bed_age     :  initial bed biodiffusity [s] [> 0]
!                      [1:NBED] 
!  init_bed_tau_crit:  initial bed sediment critical stress for erosion [N/m2] [> 0]
!                      [1:NBED] 
!  init_bed_biodiff :  initial bed sediment layer bio-diffusivity [m2/s] [> 0]
!                      [1:NBED] 
!  init_bed_mass    :  initial bed sediment mass of each size class in each bed layer   !
!                      (kg/m2).
!                      [1:NST,1:NBED] 
!--------------------------------------------------
 
real(sp), allocatable :: init_bed_porosity(:)
real(sp), allocatable :: init_bed_thickness(:) 
real(sp), allocatable :: init_bed_age(:)
real(sp), allocatable :: init_bed_biodiff(:)
real(sp), allocatable :: init_bed_tau_crit(:)
real(sp), allocatable :: init_bed_mass(:,:)
real(sp), allocatable :: init_bed_fraction(:,:) 

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

!
!-----------------------------------------------------------------------
!  Tracer identification indices.
!-----------------------------------------------------------------------
!
integer, allocatable :: idsed(:)    ! Cohesive and non-cohesive
integer, allocatable :: idmud(:)    ! Cohesive sediment
integer, allocatable :: isand(:)    ! Non-cohesive sediment



!  BED properties indices:                                             !
!  ======================                                              !
!                                                                      !
!   MBEDP           Number of bed properties (array dimension).        !
!   idBmas(:)       Sediment mass index.                               !
!   idSbed(:)       IO indices for bed properties variables.           !
!   idfrac(:)       Sediment class fraction (non-dimensional).         !
!   ithck           Sediment layer thickness (m).                      !
!   iaged           Sediment layer age (s).                            !
!   iporo           Sediment layer porosity (non-dimensional).         !
!   idiff           Sediment layer bio-diffusivity (m2/s).             !
!   ibtcr           Sediment critical stress for erosion (Pa).         !
!                                                                      !
!  BOTTOM properties indices:                                          !
!  =========================                                           !
!                                                                      !
!   MBOTP           Number of bottom properties (array dimension).     !
!   idBott(:)       IO indices for bottom properties variables.        !
!   isd50           Median sediment grain diameter (m).                !
!   idens           Median sediment grain density (kg/m3).             !
!   iwsed           Mean settling velocity (m/s).                      !

!   nthck           new total thickness of sediment layer              !
!   lthck           last thickness of sediment layer                   !
!   dthck           accumulated delta of layer thickness [m]           !
!   tmass           total mass in sediment layer [kg/m^2]              !
!   morph           differential in thickness of sediment layer        !

!   Keep it for be compatibel with original version
!   iwset           Mean settling velocity (m/s).


!   itauc           Mean critical erosion stress (m2/s2).              !
!   irlen           Sediment ripple length (m).                        !
!   irhgt           Sediment ripple height (m).                        !
!   ibwav           Bed wave excursion amplitude (m).                  !
!   izdef           Default bottom roughness (m).                      !
!   izapp           Apparent bottom roughness (m).                     !
!   izNik           Nikuradse bottom roughness (m).                    !
!   izbio           Biological bottom roughness (m).                   !
!   izbfm           Bed form bottom roughness (m).                     !
!   izbld           Bed load bottom roughness (m).                     !
!   izwbl           Bottom roughness used wave BBL (m).                !
!   iactv           Active layer thickness for erosive potential (m).  !
!   ishgt           Sediment saltation height (m).                     !
!   imaxD           Maximum inundation depth.                          !
!   isgrH           Seagrass height.                                   !
!   isgrD           Seagrass shoot density.                            !
!   idnet           Erosion or deposition.                             !
!   idoff           Offset for calculation of dmix erodibility         !
!                     profile (m).                                     !
!   idslp           Slope  for calculation of dmix or erodibility      !
!                     profile.                                         !
!   idtim           Time scale for restoring erodibility profile (s).  !
!   idbmx           Bed biodifusivity maximum.                         !
!   idbmm           Bed biodifusivity minimum.                         !
!   idbzs           Bed biodifusivity zs.                              !
!   idbzm           Bed biodifusivity zm.                              !
!   idbzp           Bed biodifusivity phi.                             !
!   idprp           Cohesive behavior.                                 !
!                                                                      !
!   nTbiom          Number of hours for depth integration              !
!=======================================================================
!-----------------------------------------------------------------------
!  bed            Sediment properties in each bed layer:               !
!                   bed(:,:,:,ithck) => layer thickness                !
!                   bed(:,:,:,iaged) => layer age                      !
!                   bed(:,:,:,iporo) => layer porosity                 !
!                   bed(:,:,:,idiff) => layer bio-diffusivity          !
!  bed_frac       Sediment fraction of each size class in each bed     !
!                   layer(nondimensional: 0-1.0).  Sum of              !
!                   bed_frac = 1.0.                                    !
!  bed_mass       Sediment mass of each size class in each bed layer   !
!                   (kg/m2).
!  bed_thick0     Sum all initial bed layer thicknesses (m).           !
!  bed_thick      Instantaneous total bed thickness (m).               !

!  bedldu         Bed load u-transport (kg/m/s).                       !
!  bedldv         Bed load v-transport (kg/m/s).                       !

!  bottom         Exposed sediment layer properties:                   !
!                   bottom(:,:,isd50) => mean grain diameter           !
!                   bottom(:,:,idens) => mean grain density            !
!                   bottom(:,:,iwsed) => mean settling velocity        !
!                   bottom(:,:,nthck) => new total thickness of layer  !  
!                   bottom(:,:,lthck) => last thickness of  layer      !
!                   bottom(:,:,dthck) => accumulated delta of layer [m]!
!                   bottom(:,:,tmass) => total mass in  layer [kg/m^2] !
!                   bottom(:,:,morph) => differential in thickness of layer
!                   bottom(:,:,itauc) => mean critical erosion stress  !
!                   bottom(:,:,irlen) => ripple length                 !
!                   bottom(:,:,irhgt) => ripple height                 !
!                   bottom(:,:,ibwav) => bed wave excursion amplitude  !
!                   bottom(:,:,izNik) => Nikuradse bottom roughness    !
!                   bottom(:,:,izbio) => biological bottom roughness   !
!                   bottom(:,:,izbfm) => bed form bottom roughness     !
!                   bottom(:,:,izbld) => bed load bottom roughness     !
!                   bottom(:,:,izapp) => apparent bottom roughness     !
!                   bottom(:,:,izwbl) => wave bottom roughness         !
!                   bottom(:,:,izdef) => default bottom roughness      !
!                   bottom(:,:,iactv) => active layer thickness        !
!                   bottom(:,:,ishgt) => saltation height              !
!                   bottom(:,:,idoff) => tau critical offset           !
!                   bottom(:,:,idslp) => tau critical slope            !
!                   bottom(:,:,idtim) => erodibility time scale        !
!                   bottom(:,:,idbmx) => diffusivity db_max            !
!                   bottom(:,:,idbmm) => diffusivity db_m              !
!                   bottom(:,:,idbzs) => diffusivity db_zs             !
!                   bottom(:,:,idbzm) => diffusivity db_zm             !
!                   bottom(:,:,idbzp) => diffusivity db_zphi           !
!                   bottom(:,:,idprp) => cohesive behavior             !
!                   bottom(:,:,isgrH) => Seagrass height               !
!                   bottom(:,:,isgrD) => Seagrass shoot density        !
!  ero_flux       Flux from erosion.                                   !
!  settling_flux  Flux from settling.                                  !
!  tcr_min         minimum shear for erosion
!  tcr_max         maximum shear for erosion
!  tcr_slp         Tau_crit profile slope
!  tcr_off         Tau_crit profile offset
!  tcr_tim         Tau_crit consolidation rate
!  transC          cohesive transition
!  transN          noncohesive transition
!
!-----------------------------------------------------------------------
!  Bed and bottom properties indices.
!-----------------------------------------------------------------------
!
!  Set size of properties arrays.
!
      integer,  pointer  :: n_bot_vars 
      integer,  pointer  :: n_bed_chars
      character(len=80), dimension(34) :: bot_snames = &
       (/       "bot_sd50"                                     &
            &  ,"bot_dens"                                     &
            &  ,"bot_wsed"                                     &
            &  ,"bot_nthck"                                    &
            &  ,"bot_lthck"                                    &
            &  ,"bot_dthck"                                    &
            &  ,"bot_tmass"                                    &
            &  ,"bot_morph"                                    &
            &  ,"bot_tauc"                                     &
            &  ,"bot_rlen"                                     &
            &  ,"bot_rhgt"                                     &
            &  ,"bot_bwav"                                     &
            &  ,"bot_zdef"                                     &
            &  ,"bot_zapp"                                     &
            &  ,"bot_zNik"                                     &
            &  ,"bot_zbio"                                     &
            &  ,"bot_zbfm"                                     &
            &  ,"bot_zbld"                                     &
            &  ,"bot_zwbl"                                     &
            &  ,"bot_actv"                                     &
            &  ,"bot_shgt"                                     &
            &  ,"bot_maxD"                                     &
            &  ,"bot_dnet"                                     &
            &  ,"bot_doff"                                     &
            &  ,"bot_dslp"                                     &
            &  ,"bot_dtim"                                     &
            &  ,"bot_dbmx"                                     &
            &  ,"bot_dbmm"                                     &
            &  ,"bot_dbzs"                                     &
            &  ,"bot_dbzm"                                     &
            &  ,"bot_dbzp"                                     &
            &  ,"bot_dprp"                                     &
            &  ,"bot_sgrD"                                     &
            &  ,"bot_sgrH"/)
      character(len=80), dimension(34) :: bot_lnames = &
      (/          "mean grain diameter"                  &
            &  ,  "mean grain density"                   &
            &  ,  "mean settle velocity"                 &
            &  ,  "new total thickness of bed layer"     &
            &  ,  "last thickness of bed layer"          &
            &  ,  "accumulated delta of lay thickness"   &
            &  ,  "total mass in sediment layer"         &
            &  ,  "differential in bed layer thickness"  &
            &  ,  "critical erosion stress"              &
            &  ,  "ripple length"                        &
            &  ,  "ripple height"                        &
            &  ,  "wave excursion amplitude"             &
            &  ,  "default bottom roughness"             &
            &  ,  "apparent bottom roughness"            &
            &  ,  "Nikuradse bottom roughness"           &
            &  ,  "biological bottom roughness"          &
            &  ,  "bed form bottom roughness"            &
            &  ,  "bed load bottom roughness"            &
            &  ,   "wave bottom roughness"               &
            &  ,  "active layer thickness"               &
            &  ,  "saltation height"                     &
            &  ,  "maximum inundation depth"             &
            &  ,  "erosion or deposition"                &
            &  ,   "tau critical offset"                 &
            &  ,   "tau critical slope"                  &
            &  ,   "erodibility time scale"              &
            &  ,   "diffusivity db_max"                  &
            &  ,   "diffusivity db_m"                    &
            &  ,   "diffusivity db_zs"                   &
            &  ,   "diffusivity db_zm"                   &
            &  ,   "diffusivity db_zphi"                 &
            &  ,   "cohesive behavior"                   &
            &  ,  "seagrass shoot density"               &
            &  ,  "seagrass height" /)
      character(len=80), dimension(34) :: bot_units = (/"m","kg&
     &/m^3","m/s","m","m","m","kg","m","N/m^2","m","m","m","m","m","m","m","m"&
     &,"m","m","m","m","m","m","N/m^2","-","s","m","m","m","m","m","&
     &-","-","m"/)

      character(len=80), dimension(5) :: bed_snames = (/"bed_thick","bed_age","bed_por","bed_diff","bed_btcr"/)
      character(len=80), dimension(5) :: bed_lnames = (/"bed layer thickness","bed layer age","bed layer porosity","bed layer bio-diffusivity","bed critical stress"/)
      character(len=80), dimension(5) :: bed_units = (/"m","days","-","-","N/m2"/)
!----------------------------------------------------------------

      integer,  target  :: MBEDP          ! Bed properties
      integer,  target  :: MBOTP          ! Bottom properties
      integer, allocatable :: idSbed(:)    ! bed properties IDs
      integer, allocatable :: idBott(:)    ! bottom properties IDs
!
!  Set properties indices.
!
      integer, parameter :: ithck = 1      ! layer thickness
      integer, parameter :: iaged = 2      ! layer age
      integer, parameter :: iporo = 3      ! layer porosity
      integer, parameter :: idiff = 4      ! layer bio-diffusivity
!for COHESIVE_BED or SED_BIODIFF or MIXED_BED
      integer, parameter :: ibtcr = 5      ! layer critical stress
!----------------------------------------------------------------

      integer, parameter :: isd50 = 1      ! mean grain diameter
      integer, parameter :: idens = 2      ! mean grain density
      integer, parameter :: iwsed = 3      ! mean settle velocity
      integer, parameter :: iwset = 3      ! mean settle velocity
      integer, parameter :: nthck = 4      ! new total thickness of sediment layer
      integer, parameter :: lthck = 5      ! last thickness of sediment layer
      integer, parameter :: dthck = 6      ! accumulated delta of layer thickness [m]
      integer, parameter :: tmass = 7      ! total mass in sediment layer [kg/m^2] 
      integer, parameter :: morph = 8      ! differential in thickness of sediment layer
      integer, parameter :: itauc = 9      ! critical erosion stress
      integer, parameter :: irlen = 10     ! ripple length
      integer, parameter :: irhgt = 11     ! ripple height
      integer, parameter :: ibwav = 12     ! wave excursion amplitude
      integer, parameter :: izdef = 13     ! default bottom roughness
      integer, parameter :: izapp = 14     ! apparent bottom roughness
      integer, parameter :: izNik = 15     ! Nikuradse bottom roughness
      integer, parameter :: izbio = 16     ! biological bottom roughness
      integer, parameter :: izbfm = 17     ! bed form bottom roughness
      integer, parameter :: izbld = 18     ! bed load bottom roughness
      integer, parameter :: izwbl = 19     ! wave bottom roughness
      integer, parameter :: iactv = 20     ! active layer thickness
      integer, parameter :: ishgt = 21     ! saltation height
      integer, parameter :: imaxD = 22     ! maximum inundation depth
      integer, parameter :: idnet = 23     ! erosion or deposition
!for COHESIVE_BED or SED_BIODIFF or MIXED_BED
      integer, parameter :: idoff = 24     ! tau critical offset
      integer, parameter :: idslp = 25     ! tau critical slope
      integer, parameter :: idtim = 26     ! erodibility time scale
      integer, parameter :: idbmx = 27     ! diffusivity db_max
      integer, parameter :: idbmm = 28     ! diffusivity db_m
      integer, parameter :: idbzs = 29     ! diffusivity db_zs
      integer, parameter :: idbzm = 30     ! diffusivity db_zm
      integer, parameter :: idbzp = 31     ! diffusivity db_zphi
!----------------------------------------------------------------
!for MIXED_BED
      integer, parameter :: idprp = 32     ! cohesive behavior
!----------------------------------------------------------------
      integer            :: isgrD          ! seagrass shoot density
      integer            :: isgrH          ! seagrass height
!----------------------------------------------------------------

      real(sp),dimension(:,:), pointer :: bottom

      real(sp), pointer :: bed(:,:,:)

!
!  Sediment metadata indices vectors.
!
      integer, allocatable :: idBmas(:)    ! class mass indices
      integer, allocatable :: idfrac(:)    ! class fraction indices
      integer, allocatable :: idUbld(:)    ! bed load u-points
      integer, allocatable :: idVbld(:)    ! bed load v-points
!
!-----------------------------------------------------------------------
!  Input sediment parameters.
!-----------------------------------------------------------------------
!
      real(sp)  :: newlayer_thick   ! deposit thickness criteria
      real(sp)  :: minlayer_thick   ! 2nd layer thickness criteria
      real(sp)  :: bedload_coeff    ! bedload rate coefficient

!      real(sp), allocatable :: Csed(:)       ! initial concentration
!      real(sp), allocatable :: Erate(:)      ! erosion rate
!      real(sp), allocatable :: Sd50(:)       ! mediam grain diameter
!      real(sp), allocatable :: Srho(:)       ! grain density
!      real(sp), allocatable :: Wsed(:)       ! settling velocity
!      real(sp), allocatable :: poros(:)      ! porosity
!      real(sp), allocatable :: tau_ce(:)     ! shear for erosion
!      real(sp), allocatable :: tau_cd(:)     ! shear for deposition
!      real(sp), allocatable :: morph_fac(:)  ! morphological factor

!for  COHESIVE_BED or  MIXED_BED
      real(sp)  :: tcr_min      ! minimum shear for erosion
      real(sp)  :: tcr_max      ! maximum shear for erosion
      real(sp)  :: tcr_slp      ! Tau_crit profile slope
      real(sp)  :: tcr_off      ! Tau_crit profile offset
      real(sp)  :: tcr_tim      ! Tau_crit consolidation rate
!-----------------------------------------------------------------------

!for  MIXED_BED
      real(sp)  :: transC       ! cohesive transition
      real(sp)  :: transN       ! noncohesive transition
!-----------------------------------------------------------------------

!for SED_BIODIFF
      real(sp)  :: Dbmx         ! Dbmax  Maximum biodiffusivity
      real(sp)  :: Dbmm         ! Dbmin  Minimum biodiffusivity
      real(sp)  :: Dbzs         ! Dbzs   Depth of maximum biodiff
      real(sp)  :: Dbzm         ! Dbzm   Depth end exp biodiff
      real(sp)  :: Dbzp         ! Dbzp   Depth of minimum biodiff
!-----------------------------------------------------------------------

!for SED_FLOCS and SED_DEFLOC
      real(sp), allocatable :: mud_frac_eq(:) ! Equilibrium fractional class distribution
      real(sp)  :: t_dfloc      ! Time scale of bed deflocculation
!-----------------------------------------------------------------------
!for  SED_BIOMASS
      integer :: nTbiom         ! Number of hours for depth integration
      real(sp)  :: sgr_diam     ! Seagrass diameter
      real(sp)  :: sgr_density  ! Seagrass density (no shoot density)
      real(sp)  :: sgr_Hthres   ! Seagrass height threshold 

      real(sp), allocatable,dimension(:) :: SgrN !  Seagrass in mols of N 
!-----------------------------------------------------------------------


     TYPE (TYPE_SEDBED) ,         target :: SEDBED

     TYPE (TYPE_SEDFLOCS),        target :: SEDFLOCS

     type(sed_type), allocatable, target :: sed(:)

!  nnew      Baroclinic (slow) time-step index corresponding to the    !
!              newest values for 3D primitive equation variables.      !
!  nstp      Baroclinic (slow) time-step index to which the current    !
!              changes are added to compute new 3D primitive equation  !
!              variables.                                              !

    integer :: nnew,nstp

    !J. Ge for tracer advection
     type(sed_type), allocatable, target :: sed0(:),sed2(:)
    !J. Ge for tracer advection

#endif

End Module Mod_CSTMS_vars
