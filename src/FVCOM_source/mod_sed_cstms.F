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
Module Mod_Sed_CSTMS
#if defined (SEDIMENT) && (CSTMS_SED)
Use Mod_Par
Use Mod_Prec 
Use Mod_Types
Use Mod_wd
Use Control, only : seddis
Use all_vars,only : CNSTNT,UNIFORM,SEDIMENT_PARAMETER_TYPE
# if defined (WAVE_CURRENT_INTERACTION)
USE MOD_WAVE_CURRENT_INTERACTION
# endif
# if defined (MULTIPROCESSOR)
Use Mod_Par
# endif
# if defined (FLUID_MUD)
USE MOD_FLUID_MUD
# endif

Use Mod_CSTMS_vars

Use Mod_FlocMod

implicit none 






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

  call alloc_sedbed

  if(SED_FLOCS)then

    CALL allocate_sedflocs

    CALL initialize_sedflocs_param (SEDFLOCS % f_mass,              &
     &                     SEDFLOCS % f_diam,                       &
     &                     SEDFLOCS % f_g1_sh,                      &
     &                     SEDFLOCS % f_g1_ds,                      &
     &                     SEDFLOCS % f_g3,                         &
     &                     SEDFLOCS % f_l1_sh,                      &
     &                     SEDFLOCS % f_l1_ds,                      &
     &                     SEDFLOCS % f_coll_prob_sh,               &
     &                     SEDFLOCS % f_coll_prob_ds,               &
     &                     SEDFLOCS % f_l3)

  end if

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
! update bottom properties                           
!----------------------------------------------------------------
  call update_bottom_properties

!----------------------------------------------------------------
! set arguments for sed names and N_sed for FVCOM to use 
! in setting up river forcing
!----------------------------------------------------------------
  N_sed = Nst
  do i=1,Nst
    sed_names(i) = sed(i)%sname
  end do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: setup_sed" 

End Subroutine Setup_Sed



!=======================================================================
!                                                                      !
!  This routine it is the main driver for the sediment-transport       !
!  model. Currently, it includes calls to the following routines:      !
!                                                                      !
!  * Vertical settling of sediment in the water column.                !
!  * Erosive and depositional flux interactions of sediment            !
!    between water column and the bed.                                 !
!  * Transport of multiple grain sizes.                                !
!  * Bed layer stratigraphy.                                           !
!  * Bed morphology.                                                   !
!  * Bedload based on Meyer Peter Mueller.                             !
!  * Bedload based on Soulsby combined waves + currents                !
!    (p166 Soulsby 1997)                                               !
!  * Bedload slope term options: Nemeth et al, 2006, Coastal           !
!    Engineering, v 53, p 265-275; Lesser et al, 2004, Coastal         !
!    Engineering, v 51, p 883-915.                                     !
!                                                                      !
!  * Seawater/sediment vertical level distribution:                    !
!                                                                      !
!         W-level  RHO-level                                           !
!                                                                      !
!            N     _________                                           !
!                 |         |                                          !
!                 |    N    |                                          !
!          N-1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |   N-1   |  A                                       !
!            2    |_________|  W                                       !
!                 |         |  A                                       !
!                 |    2    |  T                                       !
!            1    |_________|  E                                       !
!                 |         |  R                                       !
!                 |    1    |                                          !
!            0    |_________|_____ bathymetry                          !
!                 |/////////|                                          !
!                 |    1    |                                          !
!            1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |    2    |  D                                       !
!            2    |_________|  I                                       !
!                 |         |  M                                       !
!                 |  Nbed-1 |  E                                       !
!        Nbed-1   |_________|  N                                       !
!                 |         |  T                                       !
!                 |  Nbed   |                                          !
!         Nbed    |_________|                                          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
!***********************************************************************



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
  !
  !  Imported variable declarations.
  !
  Implicit None
  real(sp), intent(in ) :: DTin,Tin
  real(sp), intent(in ) :: taub_in(m) 
  integer :: ierr

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Advance_Sed_CSTMS "

  !------------------------------------------------------
  !Check for Initialization
  !------------------------------------------------------
  if(iint < sed_start) return
  nstp = 1+MOD(iint-sed_start,2)
  nnew = 3-nstp

  !------------------------------------------------------
  !Increment Sed Model Counter and Report Sed Model Init 
  !------------------------------------------------------
  sed_its = sed_its + 1
  if(msr .and.  sed_its == 1)then
     write(*,*)'========Sed Model Initiated======='
  endif

  !------------------------------------------------------
  !Set up Model Forcing  -> time step and bottom shear 
  !Convert taux/tauy to node-based values
  !Convert Model tau/rho to tau [N/m^2]
  !------------------------------------------------------
  DTsed = DTin
  T_model = Tin
# if defined (WAVE_CURRENT_INTERACTION)
  if(MB_BBL_USE.or.SG_BBL_USE.or.SSW_BBL_USE)then
     taub = taucwmax*rho_water
  else
     taub = taub_in*rho_water
  endif
# else
  taub = taub_in*rho_water
# endif

  !------------------------------------------------------
  !set bed age to current model time if first sed call 
  !------------------------------------------------------
  if(sed_its ==1)then
     bed(:,:,iaged) = T_model
  endif

  !------------------------------------------------------
  !Calculate Active Layer Thickness --> bottom(:,iactv)
  !------------------------------------------------------
  call calc_active_layer

  !-----------------------------------------------------------------------
  !  Compute sediment bedload transport.
  !-----------------------------------------------------------------------
  !
  if(BEDLOAD)then
     CALL calc_sed_bedload
  endif
  !
  !-----------------------------------------------------------------------
  !  Compute sediment suspended transport.
  !-----------------------------------------------------------------------
  if(SUSLOAD)then
     !
     !-----------------------------------------------------------------------
     !  Compute sediment flocculation
     !-----------------------------------------------------------------------
     !
     if(SED_FLOCS)then
        CALL calc_sed_flocmod
     endif
     !
     !-----------------------------------------------------------------------
     !  Compute sediment vertical settling.
     !-----------------------------------------------------------------------
     !
     CALL calc_sed_settling

     !
     !-----------------------------------------------------------------------
     !  Compute bed-water column exchanges: erosion and deposition.
     !-----------------------------------------------------------------------
     !
     CALL calc_sed_fluxes
     !

     !-----------------------------------------------------------------------
     !  Compute sediment due to biomass.
     !-----------------------------------------------------------------------
     !
     if(SED_BIOMASS)then
        CALL calc_sed_biomass
     endif

     !-----------------------------------------------------------------------
     !  Compute fluid mud dynamics.
     !-----------------------------------------------------------------------
     !
# if defined (FLUID_MID)
     call calc_fluid_mud
# endif

  endif  ! end of suspended load

  !
  !-----------------------------------------------------------------------
  !  Compute sediment bed stratigraphy.
  !-----------------------------------------------------------------------
  !
  if (COHESIVE_BED .or. MIXED_BED)then
     CALL calc_sed_bed_cohesive
  elseif(NONCOHESIVE_BED2)then
     CALL calc_sed_bed2
  else
     CALL calc_sed_bed
  endif
  !
  !-----------------------------------------------------------------------
  !  Compute sediment bed biodiffusivity.
  !-----------------------------------------------------------------------
  !
  if(SED_BIODIFF)then
     CALL calc_sed_biodiff
  endif
  !
  !-----------------------------------------------------------------------
  !  Compute sediment surface layer properties.
  !-----------------------------------------------------------------------
  !
  CALL Update_Bottom_Properties

  !
  !-----------------------------------------------------------------------
  !  Horizontal Advection of Sediment Concentration 
  !-----------------------------------------------------------------------
  !
  CALL calc_sed_advection

  !
  !------------------------------------------------------
  !Vertical Diffusion of Sediment Concentrations
  !------------------------------------------------------
  !
  CALL calc_sed_diffusion

  !------------------------------------------------------
  ! Set Open Boundary Conditions 
  !------------------------------------------------------

  CALL calc_sed_bcond

  !------------------------------------------------------
  ! Update Depth Delta and Morphology Array
  ! this array is positive if material has accumulated,
  ! negative if eroded
  !------------------------------------------------------
  call update_thickness_delta

  !------------------------------------------------------
  ! Update Concentration Variables to next time level
  !------------------------------------------------------
  !
  CALL finalize_sed_step

  !------------------------------------------------------
  ! Report the status of sediment
  !------------------------------------------------------
  !

  if(mod(sed_its,n_report)==0)call sed_report   

  RETURN

END SUBROUTINE Advance_Sed



!=======================================================================
!Read Sediment Parameters From Sediment Input File
!=======================================================================
Subroutine Read_Sed_Params
  Use Input_Util
  Implicit None
  integer linenum,i,k1,iscan
  real(sp)           :: ftemp
  character(len=120) :: stemp

  real(sp),allocatable:: ncstmp(:)
  real(sp),allocatable:: nnstmp(:)
  real(sp),allocatable:: nsttmp(:)
  logical, allocatable:: nswitch(:)
  character(len=80),allocatable:: strtmp(:)

  !
  !  Imported variable declarations
  !

  !
  !  Local variable declarations.
  !
  integer :: iTrcStr, iTrcEnd
  integer :: ifield, igrid, itracer, itrc, nline


  !-----------------------------------------------------------------------
  !  Initialize.
  !-----------------------------------------------------------------------
  !
  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Read_Sed_Params " 
  !
  !-----------------------------------------------------------------------
  !  Read in cohesive and non-cohesive model parameters.
  !-----------------------------------------------------------------------
  !

  !read in interation interval for reporting 
  Call Get_Val(n_report,sedfile,'N_REPORT',line=linenum,echo=.true.) 

  !read in Depositional bed layer thickness criteria to create a new layer
  Call Get_Val(newlayer_thick,sedfile,'NEWLAYER_THICK',line=linenum,echo=.true.)

  !read in number of cohesive sediment class
  Call Get_Val(ncs,sedfile,'NCS',line=linenum,echo=.true.)

  !read in number of non-cohesive sediment class
  Call Get_Val(nns,sedfile,'NNS',line=linenum,echo=.true.)

  !total number of classes for cohesive and non-cohesive sediment
  NST = NCS+NNS

  if(nst < 1)then
    write(*,*)'sediment input file must have at least one sediment class'
    write(*,*)'currently nst = ',nst
    stop
  endif

  NSED = NST  ! to be consistent with original output code.

  allocate(sed(nst))

  !J. Ge for tracer advection
  Allocate(sed0(nst))
  Allocate(sed2(nst))
  !J. Ge for tracer advection



  !read in switches for bedload calculation
  Call Get_Val(BEDLOAD,sedfile,'BEDLOAD',line=linenum,echo=.true.)

  !read in switches for suspended sediment calculation
  Call Get_Val(SUSLOAD,sedfile,'SUSLOAD',line=linenum,echo=.true.)

  !read in minumum sediment density
  Call Get_Val(min_Srho,sedfile,'MIN_SRHO',line=linenum,echo=.false.)

  !read in switches for 1-D sediment calculation
  Call Get_Val(oned_model,sedfile,'SED_ONED',line=linenum,echo=.true.)

  !read in switches for cohesive bed layer stratigraphy
  Call Get_Val(COHESIVE_BED,sedfile,'COHESIVE_BED',line=linenum,echo=.true.)

  if(COHESIVE_BED .and. NNS>0 )then
     write(*,*)'cohesive bed dynamics has been specified, however,the'
     write(*,*)'number of non-cohesive sediment classes are not-zero'
     write(*,*)'currently nns = ',nns
     stop     
  end if


  !read in switches for non-cohesive bed layer stratigraphy
  Call Get_Val(NONCOHESIVE_BED2,sedfile,'NONCOHESIVE_BED2',line=linenum,echo=.true.)

  !read in switches for mixed-bed layer stratigraphy
  Call Get_Val(MIXED_BED,sedfile,'MIXED_BED',line=linenum,echo=.true.)

  !read in switches for morphology calculation
  Call Get_Val(SED_MORPH,sedfile,'SED_MORPH',line=linenum,echo=.true.)

  !read in switches for floc bed layer stratigraphy
  Call Get_Val(SED_FLOCS,sedfile,'SED_FLOCS',line=linenum,echo=.true.)

  !read in switches for sediment biomass due to vegation growth 
  Call Get_Val(SED_BIOMASS,sedfile,'SED_BIOMASS',line=linenum,echo=.true.)

  !read in switches for sediment bed layer mixing (biodiffusion)
  Call Get_Val(SED_BIODIFF,sedfile,'SED_BIODIFF',line=linenum,echo=.true.)

  !read in switches for sediment defloccuation
  Call Get_Val(SED_DEFLOC,sedfile,'SED_DEFLOC',line=linenum,echo=.true.)


  if(SED_BIOMASS)then
     !read in switches for bottom seagrass influences
     Call Get_Val(SEAGRASS_BOTTOM ,sedfile,'SEAGRASS_BOTTOM ',line=linenum,echo=.true.)

     !read in switches for bottom seagrass sink
     Call Get_Val(SEAGRASS_SINK,sedfile,'SEAGRASS_SINK ',line=linenum,echo=.true.)
  endif

  if(MIXED_BED)then
     !read in cohesive transition
     Call Get_Val(transC,sedfile,'TRANSC',line=linenum,echo=.true.)

     !read in noncohesive transition
     Call Get_Val(transN,sedfile,'TRANSN',line=linenum,echo=.true.)
  endif

  if(BEDLOAD)then

     !read in switches for Meyer-Peter Mueller formulation
     Call Get_Val(BEDLOAD_MPM,sedfile,'BEDLOAD_MPM',line=linenum,echo=.true.)

     !read in switches for Soulsby formulation
     Call Get_Val(BEDLOAD_SOULSBY,sedfile,'BEDLOAD_SOULSBY',line=linenum,echo=.true.)

     !read in switches for memeth slope calculation
     Call Get_Val(SLOPE_NEMETH,sedfile,'SLOPE_NEMETH',line=linenum,echo=.true.)

     !read in switches for lesser slope calculation
     Call Get_Val(SLOPE_LESSER,sedfile,'SLOPE_LESSER',line=linenum,echo=.true.)

  endif


  !read in switch for constant tau cd
  Call Get_Val(sed_tau_cd_const,sedfile,'SED_TAU_CD_CONST',line=linenum,echo=.true.)

  !read in switch for linear tau cd
  Call Get_Val(sed_tau_cd_lin,sedfile,'SED_TAU_CD_LIN',line=linenum,echo=.true.)

  !read in Bed load transport rate coefficient
  Call Get_Val(bedload_coeff,sedfile,'BEDLOAD_COEFF',line=linenum,echo=.true.)

  !read in linear continuation for vertical settling
  Call Get_Val(linear_continuation,sedfile,'LINEAR_CONTINUATION',line=linenum,echo=.true.)

  !read in neumann continuation for vertical settling
  Call Get_Val(neumann,sedfile,'NEUMANN',line=linenum,echo=.true.)

  if(SED_BIODIFF)then

     !read in setting biodiffusivity profile
     Call Get_Val(DB_PROFILE,sedfile,'DB_PROFILE',line=linenum,echo=.true.)    

     !read in Maximum biodiffusivity
     Call Get_Val(Dbmx,sedfile,'DBMAX',line=linenum,echo=.true.)

     !read in Minimum biodiffusivity
     Call Get_Val(Dbmm,sedfile,'DBMIN',line=linenum,echo=.true.)

     !read in Depth of maximum biodiffusivity
     Call Get_Val(Dbzs,sedfile,'DBZS',line=linenum,echo=.true.)

     !read in  Depth of end of exponential biodiffusivity
     Call Get_Val(Dbzm,sedfile,'DBZM',line=linenum,echo=.true.)

     !read in Depth of minimum biodiffusivity
     Call Get_Val(Dbzp,sedfile,'DBZP',line=linenum,echo=.true.)

  endif

  if(ncs>0)then

     sed(1:ncs)%stype = 'cohesive'
     sed(ncs+1:nst)%stype = 'non-cohesive'

     allocate(ncstmp(ncs))

     !read in sediment names for Cohesive sediment
     allocate(strtmp(ncs))
     Call Get_Val_Array(strtmp,sedfile,'MUD_NAME',ncs,echo=.true.)
     sed(1:ncs)%sname=strtmp
     do i=1,ncs
        sed(i)%sname2=trim(strtmp(i))//'_bedload'
     end do
     deallocate(strtmp)


     !read in Median grain diameter (mm) for Cohesive sediment
     Call Get_Val_Array(ncstmp,sedfile,'MUD_SD50',ncs,echo=.true.)
     sed(1:ncs)%Sd50=ncstmp*0.001  !convert settling velocity to m/s from input mm/s

     !read in  Sediment concentration (kg/m3).
     Call Get_Val_Array(ncstmp,sedfile,'MUD_CSED',ncs,echo=.true.)
     sed(1:ncs)%Csed_initial=ncstmp

     !read in Sediment grain density (kg/m3).
     Call Get_Val_Array(ncstmp,sedfile,'MUD_SRHO',ncs,echo=.true.)
     sed(1:ncs)%Srho=ncstmp

     !read in Particle settling velocity (mm/s).
     Call Get_Val_Array(ncstmp,sedfile,'MUD_WSED',ncs,echo=.true.)
     sed(1:ncs)%Wsed=ncstmp*0.001 !convert settling velocity to m/s from input mm/s

     !read in Surface erosion rate (kg/m2/s).
     Call Get_Val_Array(ncstmp,sedfile,'MUD_ERATE',ncs,echo=.true.)
     sed(1:ncs)%Erate=ncstmp

     !read in Critical shear for erosion and deposition (N/m2).
     Call Get_Val_Array(ncstmp,sedfile,'MUD_TAU_CE',ncs,echo=.true.)
     sed(1:ncs)%tau_ce=ncstmp

     Call Get_Val_Array(ncstmp,sedfile,'MUD_TAU_CD',ncs,echo=.true.)
     sed(1:ncs)%tau_cd=ncstmp

     !read in  Porosity (nondimensional: 0.0-1.0):  Vwater/(Vwater+Vsed).
     Call Get_Val_Array(ncstmp,sedfile,'MUD_POROS',ncs,echo=.true.)
     sed(1:ncs)%poros=ncstmp

     !read in Harmonic/biharmonic horizontal diffusion of tracer for nonlinear model
     Call Get_Val_Array(ncstmp,sedfile,'MUD_TNU2',ncs,echo=.true.)
     sed(1:ncs)%nl_tnu2=ncstmp

     Call Get_Val_Array(ncstmp,sedfile,'MUD_TNU4',ncs,echo=.true.)
     sed(1:ncs)%nl_tnu4=ncstmp

     ! read in Logical switches (TRUE/FALSE) to increase horizontal diffusivity
     ! of cohesive sediment trace in specific areas of the application
     ! domain (like sponge areas)
     allocate(nswitch(ncs)) 
     Call Get_Val_Array(nswitch,sedfile,'MUD_Sponge',ncs,echo=.true.)
     sed(1:ncs)%LtracerSponge=nswitch

     !read in Vertical mixing coefficients for tracers in nonlinear model and
     ! basic state scale factor in adjoint-based algorithms.
     Call Get_Val_Array(ncstmp,sedfile,'MUD_AKT_BAK',ncs,echo=.true.)
     sed(1:ncs)%Akt_bak=ncstmp


     !read in  Nudging/relaxation time scales,
     Call Get_Val_Array(ncstmp,sedfile,'MUD_TNUDG',ncs,echo=.true.)
     sed(1:ncs)%Tnudg=ncstmp

     !read in Morphological time scale factor (greater than or equal to 1.0)
     Call Get_Val_Array(ncstmp,sedfile,'MUD_MORPH_FAC',ncs,echo=.true.)
     sed(1:ncs)%morph_fac=ncstmp

     !read in Logical switches (TRUE/FALSE) to specify which variables to consider on
     ! tracers point Sources/Sinks (like river runoff). 
     Call Get_Val_Array(nswitch,sedfile,'MUD_Ltracer',ncs,echo=.true.)
     sed(1:ncs)%LtracerSrc=nswitch


     !read in Logical switches (T/F) to process cohesive sediment tracer climatology. 
     Call Get_Val_Array(nswitch,sedfile,'MUD_Ltclm',ncs,echo=.true.)
     sed(1:ncs)%LtracerCLM=nswitch

     !read in  Logical switches (T/F) to activate nugding of cohesive
     !sediment tracer climatology. 
     Call Get_Val_Array(nswitch,sedfile,'MUD_Tnudge',ncs,echo=.true.)
     sed(1:ncs)%LnudgeTCLM=nswitch

     deallocate(nswitch)

  endif


  if(COHESIVE_BED.or.MIXED_BED)then

     !read in  minimum shear for erosion
     Call Get_Val(tcr_min,sedfile,'MUD_TAUCR_MIN',line=linenum,echo=.true.) 

     !read in  maximum shear for erosion
     Call Get_Val(tcr_max,sedfile,'MUD_TAUCR_MAX',line=linenum,echo=.true.)

     !read in Tau_crit profile slope
     Call Get_Val(tcr_slp,sedfile,'MUD_TAUCR_SLOPE',line=linenum,echo=.true.)

     !read in Tau_crit profile offset
     Call Get_Val(tcr_off,sedfile,'MUD_TAUCR_OFF',line=linenum,echo=.true.)

     !read in  Tau_crit consolidation rate
     Call Get_Val(tcr_tim,sedfile,'MUD_TAUCR_TIME',line=linenum,echo=.true.)

  endif


  if(nns>0)then

     allocate(nswitch(nns))

     allocate(nnstmp(nns))

     !read in sediment names for non-Cohesive sediment
     allocate(strtmp(nns))
     Call Get_Val_Array(strtmp,sedfile,'SAND_NAME',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%sname=strtmp
     do i=1,nns
        sed(ncs+i)%sname2=trim(strtmp(i))//'_bedload'
     end do

     deallocate(strtmp)

     !read in Median grain diameter (mm) for non-cohesive sediment
     Call Get_Val_Array(nnstmp,sedfile,'SAND_SD50',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%Sd50=nnstmp*.001  !convert mean grain diameter from mm to m 

     !read in Sediment concentration (kg/m3) for non-cohesive sediment
     Call Get_Val_Array(nnstmp,sedfile,'SAND_CSED',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%Csed_initial=nnstmp

     !read in grain density (kg/m3) for non-cohesive sediment
     Call Get_Val_Array(nnstmp,sedfile,'SAND_SRHO',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%Srho=nnstmp

     !read in Particle settling velocity (mm/s) for non-cohesive sediment
     Call Get_Val_Array(nnstmp,sedfile,'SAND_WSED',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%Wsed=nnstmp*.001 !convert settling velocity to m/s from input mm/s

     !read in  Surface erosion rate (kg/m2/s) for non-cohesive sediment
     Call Get_Val_Array(nnstmp,sedfile,'SAND_ERATE',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%Erate=nnstmp

     !read in Critical shear for erosion and deposition (N/m2).
     Call Get_Val_Array(nnstmp,sedfile,'SAND_TAU_CE',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%tau_ce=nnstmp

     Call Get_Val_Array(nnstmp,sedfile,'SAND_TAU_CD',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%tau_cd=nnstmp

     !read in Porosity (nondimensional: 0.0-1.0):  Vwater/(Vwater+Vsed).
     Call Get_Val_Array(nnstmp,sedfile,'SAND_POROS',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%poros=nnstmp

     !read in Harmonic/biharmonic horizontal diffusion of tracer
     Call Get_Val_Array(nnstmp,sedfile,'SAND_TNU2',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%nl_tnu2=nnstmp

     Call Get_Val_Array(nnstmp,sedfile,'SAND_TNU4',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%nl_tnu4=nnstmp

     !read in Logical switches (TRUE/FALSE) to increase horizontal diffusivity
     ! of noncohesive sediment trace
     Call Get_Val_Array(nswitch,sedfile,'SAND_Sponge',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%LtracerSponge=nswitch

     !read in Vertical mixing coefficients for tracers in nonlinear model and
     ! basic state scale factor
     Call Get_Val_Array(nnstmp,sedfile,'SAND_AKT_BAK',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%Akt_bak=nnstmp 

     !read in Nudging/relaxation time scales,
     Call Get_Val_Array(nnstmp,sedfile,'SAND_TNUDG',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%Tnudg=nnstmp

     !read in Morphological time scale factor (greater than or equal to 1.0)
     Call Get_Val_Array(nnstmp,sedfile,'SAND_MORPH_FAC',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%morph_fac=nnstmp

     !read in  Logical switches (TRUE/FALSE) to activate tracers point Sources/Sinks
     ! (like river runoff) and to specify which tracer variables to consider
     Call Get_Val_Array(nswitch,sedfile,'SAND_Ltsrc',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%LtracerSrc=nswitch

     !read in Logical switches (TRUE/FALSE) to read and process
     ! noncohesive sediment tracers climatology
     Call Get_Val_Array(nswitch,sedfile,'SAND_Ltclm',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%LtracerCLM=nswitch

     !read in Logical switches (TRUE/FALSE) to nudge the desired noncohesive sediment
     ! tracer climatology field.
     Call Get_Val_Array(nswitch,sedfile,'SAND_Tnudge',nns,echo=.true.)
     sed(ncs+1:ncs+nns)%LnudgeTCLM=nswitch 

  endif

  !------------------------------------------------------------------------------
  ! Flocculation Sediment Parameters, [1] values expected.
  !------------------------------------------------------------------------------

  if(SED_FLOCS)then

     !read in l_ADS  : boolean set to .true. if differential settling aggregation 
     Call Get_Val(l_ADS,sedfile,'l_ADS',line=linenum,echo=.true.) 

     !read in l_ASH : boolean set to .true. if shear aggregation
     Call Get_Val(l_ASH,sedfile,'l_ASH',line=linenum,echo=.true.) 

     !read in l_COLLFRAG  : boolean set to .true. if collision-induced fragmentation enable
     Call Get_Val(l_COLLFRAG,sedfile,'l_COLLFRAG',line=linenum,echo=.true.) 

     !read in f_dp0 : Primary particle size (m), typically 4e-6 m
     Call Get_Val(f_dp0,sedfile,'f_dp0',line=linenum,echo=.true.) 

     !read in f_nf : Floc fractal dimension, typically ranging from 1.6 to 2.6
     Call Get_Val(f_nf,sedfile,'f_nf',line=linenum,echo=.true.) 

     !read in  f_dmax : (unused) Maximum diameter (m)
     Call Get_Val(f_dmax,sedfile,'f_dmax',line=linenum,echo=.true.) 

     !read in  f_nb_frag : number of fragments by shear erosion. If binary/ternary : 2.
     Call Get_Val(f_nb_frag,sedfile,'f_nb_frag',line=linenum,echo=.true.) 

     !read in f_alpha :  flocculation efficiency, ranging from 0 .to 1.
     Call Get_Val(f_alpha,sedfile,'f_alpha',line=linenum,echo=.true.) 

     !read in f_beta : shear fragmentation rate
     Call Get_Val(f_beta,sedfile,'f_beta',line=linenum,echo=.true.) 

     !read in  f_ater : for ternary breakup, use 0.5, for binary : 0
     Call Get_Val(f_ater,sedfile,'f_ater',line=linenum,echo=.true.) 

     !read in  f_ero_frac :  fraction of the shear fragmentation term
     !transfered to shear erosion. Ranging from 0. (no erosion) to 1.
     !(all erosion)
     Call Get_Val(f_ero_frac,sedfile,'f_ero_frac',line=linenum,echo=.true.) 

     !read in  f_ero_nbfrag : Number of fragments induced by shear erosion.   
     Call Get_Val(f_ero_nbfrag,sedfile,'f_ero_nbfrag',line=linenum,echo=.true.) 

     !read in f_ero_iv : fragment size class (could be changed to a
     ! particle  size or a particle distribution (INTEGER)
     Call Get_Val(f_ero_iv,sedfile,'f_ero_iv',line=linenum,echo=.true.) 

     !read in f_collfragparam : Fragmentation rate for collision
     !-induced breakup
     Call Get_Val(f_collfragparam,sedfile,'f_collfragparam',line=linenum,echo=.true.) 

     !read in  f_clim : min concentration below which flocculation
     ! processes are not calculated
     Call Get_Val(f_clim,sedfile,'f_clim',line=linenum,echo=.true.) 

     !read in l_testcase - if .TRUE. sets G(t) to values from Verney
     ! et al., 2011 lab experiment
     Call Get_Val(l_testcase,sedfile,'l_testcase',line=linenum,echo=.true.)

     !read in gls_p
     Call Get_Val(gls_p,sedfile,'GLS_P',line=linenum,echo=.true.)

     !read in gls_m
     Call Get_Val(gls_m,sedfile,'GLS_M',line=linenum,echo=.true.)

     !read in gls_n
     Call Get_Val(gls_n,sedfile,'GLS_N',line=linenum,echo=.true.)

     !read in gls_cmu0
     Call Get_Val(gls_cmu0,sedfile,'GLS_CMU0',line=linenum,echo=.true.)

     !read in FLOC_TURB_DISS
     Call Get_Val(FLOC_TURB_DISS,sedfile,'FLOC_TURB_DISS',line=linenum,echo=.true.)

     !read in FLOC_BBL_DISS
     Call Get_Val(FLOC_BBL_DISS,sedfile,'FLOC_BBL_DISS',line=linenum,echo=.true.)

  endif

  if(SED_FLOCS.and.SED_DEFLOC)then
     allocate(mud_frac_eq(ncs))
     !read in Equilibrium fractional class distribution (they should add
     ! up to 1).
     Call Get_Val_Array(mud_frac_eq,sedfile,'MUD_FRAC_EQ',ncs,echo=.true.)
     mud_frac_eq(1:ncs)=mud_frac_eq

     !read in Time scale of flocculation descomposition in bed
     Call Get_Val(t_dfloc,sedfile,'t_dfloc',line=linenum,echo=.true.)

  endif

  if(SED_BIOMASS)then


     !read in Number of hours for depth integration
     Call Get_Val(nTbiom,sedfile,'nTbiom',line=linenum,echo=.true.)

     !read in 
     Call Get_Val(sgr_diam,sedfile,'SGR_DIAM',line=linenum,echo=.true.)

     !read in 
     Call Get_Val(sgr_density,sedfile,'SGR_DENS',line=linenum,echo=.true.)

     !read in 
     Call Get_Val(sgr_Hthres,sedfile,'SGR_HTHRES',line=linenum,echo=.true.)

  endif

  !
  !-----------------------------------------------------------------------
  !  Report input parameters.
  !-----------------------------------------------------------------------
  !

  !read in start interval for sed model 
  Call Get_Val(sed_start,sedfile,'SED_START',line=linenum,echo=.true.)  

  !hot start option
  Call Get_Val(sed_hot_start,sedfile,'SED_HOT_START',line=linenum,echo=.true.)  

  !read in morphology parameters
  Call Get_Val(morpho_model,sedfile,'MORPHO_MODEL',line=linenum,echo=.true.)
  Call Get_Val(morpho_factor,sedfile,'MORPHO_FACTOR',line=linenum,echo=.true.)  
  if(morpho_model)then
     Call Get_Val(morpho_incr,sedfile,'MORPHO_INCR',line=linenum,echo=.true.)  
     Call Get_Val(morpho_strt,sedfile,'MORPHO_STRT',line=linenum,echo=.true.)  
  else
     morpho_incr   = 1
  endif

  !read in number of bed layers 
  Call Get_Val(nbed,sedfile,'NBED',line=linenum,echo=.true.)
  if(nbed < 1)then
    write(*,*)'sediment input file must have at least one bed layer'
    write(*,*)'currently nbed = ',nbed
    stop
  endif

  !read in logical for infinite sediment supply 
  Call Get_Val(inf_bed,sedfile,'INF_BED',line=linenum,echo=.true.)

  !read in initial bed porosity  
  allocate(init_bed_porosity(nbed))
  if(nbed==1)then
     Call Get_Val(ftemp,sedfile,'INIT_BED_POROSITY',line=linenum,echo=.true.)
     init_bed_porosity(1) = ftemp
  else
     open(unit=999,file=trim(sedfile),form='formatted')
     iscan = scan_file(999,"INIT_BED_POROSITY",fvec = init_bed_porosity,nsze = nbed)
     if(iscan /= 0)then
        write(*,*)'problem reading INIT_BED_POROSITY from sediment param file'
        stop
     endif
     close(999)
  endif

  if(minval(init_bed_porosity) < 0. .or. maxval(init_bed_porosity) > 1)then
     write(*,*)'error in init_bed_porosity in sed param file'
     write(*,*)'you entered: ',init_bed_porosity
     write(*,*)'must be in range [0,1]'
     stop
  endif

  !read in initial bed thickness
  allocate(init_bed_thickness(nbed))
  if(nbed==1)then
     Call Get_Val(ftemp,sedfile,'INIT_BED_THICKNESS',line=linenum,echo=.true.)
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

  !read in initial bed age
  allocate(init_bed_age(nbed))
  if(nbed==1)then
     Call Get_Val(ftemp,sedfile,'INIT_BED_AGE',line=linenum,echo=.true.)
     init_bed_age(1) = ftemp
  else
     open(unit=999,file=trim(sedfile),form='formatted')
     iscan = scan_file(999,"INIT_BED_AGE",fvec = init_bed_age,nsze = nbed)
     if(iscan /= 0)then
        write(*,*)'problem reading INIT_BED_AGE from sediment param file'
        stop
     endif
     close(999)
  endif
  if(minval(init_bed_age) <= 0.)then
     write(*,*)'error in init_bed_age in sed param file'
     write(*,*)'you entered: ',init_bed_age
     write(*,*)'must be > 0'
  endif

  !read in initial bed bio-diffusivity
  if(SED_BIODIFF)then
     allocate(init_bed_biodiff(nbed))
     if(nbed==1)then
        Call Get_Val(ftemp,sedfile,'INIT_BED_BIODIFF',line=linenum,echo=.true.)
        init_bed_biodiff(1) = ftemp
     else
        open(unit=999,file=trim(sedfile),form='formatted')
        iscan = scan_file(999,"INIT_BED_BIODIFF",fvec = init_bed_biodiff,nsze = nbed)
        if(iscan /= 0)then
           write(*,*)'problem reading INIT_BED_BIODIFF from sediment param file'
           stop
        endif
        close(999)
     endif
     if(minval(init_bed_age) <= 0.)then
        write(*,*)'error in init_bed_bioddiff in sed param file'
        write(*,*)'you entered: ',init_bed_biodiff
        write(*,*)'must be > 0'
     endif
  end if


  !read Sediment critical stress for erosion
  allocate(init_bed_tau_crit(nbed))
  if(nbed==1)then
     Call Get_Val(ftemp,sedfile,'INIT_BED_TAU_CRIT',line=linenum,echo=.true.)
     init_bed_tau_crit(1) = ftemp
  else
     open(unit=999,file=trim(sedfile),form='formatted')
     iscan = scan_file(999,"INIT_BED_TAU_CRIT",fvec = init_bed_tau_crit,nsze = nbed)
     if(iscan /= 0)then
        write(*,*)'problem reading INIT_BED_TAU_CRIT from sediment param file'
        stop
     endif
     close(999)
  endif
  if(minval(init_bed_tau_crit) <= 0.)then
     write(*,*)'error in init_bed_tau_crit in sed param file'
     write(*,*)'you entered: ',init_bed_tau_crit
     write(*,*)'must be > 0'
  endif

  !read in initial bed fraction   
  allocate(init_bed_fraction(nst,nbed))
  allocate(nsttmp(nst*nbed))
  Call Get_Val_Array(nsttmp,sedfile,'INIT_BED_FRACTION',nst*nbed,echo=.true.)
  do i=1,nbed
     init_bed_fraction(1:nst,i)=nsttmp((i-1)*nst+1:(i-1)*nst+nst)
     if(sum(init_bed_fraction(1:nst,i))/=1.0)then
        write(*,*)'error in init_bed_fraction in sed param file'
        write(*,*)'also error in your chosen career'
        write(*,*)'in bed layer:',i 
        write(*,*)'you entered: ',init_bed_fraction(1:nst,i)
        write(*,*)'must have the summary = 1'
        stop
     end if
  end do
  if(minval(init_bed_fraction) < 0. .or. maxval(init_bed_fraction) > 1.)then
     write(*,*)'error in init_bed_fraction in sed param file'
     write(*,*)'also error in your chosen career'
     write(*,*)'you entered: ',init_bed_fraction  
     write(*,*)'must be >= 0 and <= 1'
     stop
  endif
  deallocate(nsttmp)

  !read in nudging switch 
  Call Get_Val(sed_nudge,sedfile,'SED_NUDGE',line=linenum,echo=.true.)

  !read in nudging relaxation factor 
  if(sed_nudge) then
     Call Get_Val(sed_alpha,sedfile,'SED_ALPHA',line=linenum,echo=.false.)
     Call Get_Val(sed_ramp ,sedfile,'SED_RAMP ',line=linenum,echo=.false.)
  endif

  !read in params controlling output 
  Call Get_Val(sed_dumpbed,sedfile,'SED_DUMPBED',line=linenum,echo=.true.)
  Call Get_Val(sed_dumpbot,sedfile,'SED_DUMPBOT',line=linenum,echo=.true.)

  !hindered effect on vertical diffusion
  Call Get_Val(VERT_HINDERED,sedfile,'VERT_HINDERED',line=linenum,echo=.true.)

  !read in point source switch 
  Call Get_Val(sed_source,sedfile,'SED_PTSOURCE',line=linenum,echo=.true.)

  !check values
  if(nst < 1)then
     write(*,*)'sediment input file must have at least one sediment class'
     write(*,*)'currently nst = ',nst
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

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Read_Sed_Params " 

End Subroutine Read_Sed_Params




!=======================================================================
! Allocate Sediment Variables                      
!=======================================================================
  Subroutine Alloc_Sed_Vars
  use lims,     only: nt,mt,kb,numqbc,kbm1
!  use mod_obcs, only: iobcn
  implicit none
  integer i,ic,tmp1,tmp2,tmp3

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Alloc_Sed_Vars "

  if(COHESIVE_BED.or.SED_BIODIFF.or. MIXED_BED)then
     MBEDP = 5      ! Bed properties
  else
     MBEDP = 4      ! Bed properties
  endif
  allocate(idSbed(MBEDP))


  if(SEAGRASS_BOTTOM .and. MIXED_BED)then
      MBOTP = 34     ! Bottom properties
  elseif(SEAGRASS_BOTTOM .and. .not.(MIXED_BED))then
      MBOTP = 25     ! Bottom properties
  elseif(COHESIVE_BED.or.SED_BIODIFF)then
      MBOTP = 31     ! Bottom properties
  elseif(MIXED_BED)then
      MBOTP = 32     ! Bottom properties
  else
      MBOTP = 23     ! Bottom properties
  endif
  allocate(idBott(MBOTP))

  if(SED_BIOMASS.and.SEAGRASS_BOTTOM.and.MIXED_BED)then
     isgrD = 33     ! seagrass shoot density
     isgrH = 34     ! seagrass height
  elseif(SED_BIOMASS.and.SEAGRASS_BOTTOM.and.(.not.(MIXED_BED)))then 
     isgrD = 24     ! seagrass shoot density
     isgrH = 25     ! seagrass height
  endif

  n_bot_vars  => MBOTP

  n_bed_chars => MBEDP    

!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
!  Allocate various nested grid depended parameters
!

  if(COHESIVE_BED.or.MIXED_BED)then
    tcr_min = 0.0_sp
    tcr_max = 0.0_sp
    tcr_slp = 0.0_sp
    tcr_off = 0.0_sp
    tcr_tim = 0.0_sp
  endif

  if(MIXED_BED)then
    transC = 0.0_sp
    transN = 0.0_sp
  endif

  if(SED_BIOMASS)then
    sgr_diam = 0.0_sp
    sgr_density = 0.0_sp
    sgr_Hthres = 0.0_sp
  endif
!
!  Allocate sediment tracers indices vectors.
!

  allocate ( idsed(MAX(1,NST)) )
  allocate ( idmud(MAX(1,NCS)) )
  allocate ( isand(MAX(1,NNS)) )
  allocate ( idBmas(NST) )
  allocate ( idfrac(NST) )
  allocate ( idUbld(NST) )
  allocate ( idVbld(NST) )

!
!  Set cohesive and noncohesive suspended sediment tracers
!  identification indices.
!
      ic=0
      DO i=1,NCS
        ic=ic+1
        idmud(i)=ic
        idsed(i)=idmud(i)
      END DO
      DO i=1,NNS
        ic=ic+1
        isand(i)=ic
        idsed(NCS+i)=isand(i)
      END DO

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Alloc_Sed_Vars "

  End Subroutine Alloc_Sed_Vars



  SUBROUTINE Alloc_sedbed
  Use Scalar
  Use Control, only : ireport,iint,msr,par
  Use Lims,    only : m,mt,nt,kbm1,numqbc,kb,nprocs,myid

# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif

  implicit none
  integer i,tmp1,tmp2,tmp3


!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Alloc_sedbed"
!

  !ensure arrays have nonzero dimension
  tmp1 = max(numqbc,1)
  tmp2 = max(Nobc+1,1)
  tmp3 = max(Nobc  ,1)

  if(BEDLOAD.and.SED_DUMPBED)then
     allocate ( SEDBED% avgbedldu(0:mt,NST) )      ;SEDBED% avgbedldu = 0.0_sp
     allocate ( SEDBED% avgbedldv(0:mt,NST) )      ;SEDBED% avgbedldv = 0.0_sp
  endif

  allocate ( SEDBED % bed(0:mt,Nbed,MBEDP) )       ;SEDBED % bed = 0.0_sp
  allocate ( SEDBED % bed_frac(0:mt,Nbed,NST) )    ;SEDBED % bed_frac = 0.0_sp
  allocate ( SEDBED % bed_mass(0:mt,Nbed,2,NST) )  ;SEDBED % bed_mass = 0.0_sp

  if(SED_MORPH)then
    allocate ( SEDBED % bed_thick0(0:mt) )         ;SEDBED % bed_thick0 = 0.0_sp
    allocate ( SEDBED % bed_thick(0:mt,1:2)  )         ;SEDBED % bed_thick = 0.0_sp
  endif

  if(BEDLOAD)then
    allocate ( SEDBED % bedldu(0:mt,NST) )         ;SEDBED % bedldu = 0.0_sp
    allocate ( SEDBED % bedldv(0:mt,NST) )         ;SEDBED % bedldv = 0.0_sp
  end if

  allocate ( SEDBED % bottom(0:mt,MBOTP) )         ;SEDBED % bottom = 0.0_sp
 
  bottom => SEDBED % bottom
  bed    => SEDBED % bed

  do i=1,nst
    allocate(sed(i)%conc(0:mt,kb  ))      ; sed(i)%conc = 0.0_sp
    allocate(sed(i)%cnew(0:mt,kb  ))      ; sed(i)%cnew = 0.0_sp
    allocate(sed(i)%mass(0:mt,nbed))      ; sed(i)%mass = 0.0_sp
    allocate(sed(i)%frac(0:mt,nbed))      ; sed(i)%frac = 0.0_sp
    allocate(sed(i)%bflx(0:mt     ))      ; sed(i)%bflx = 0.0_sp
    allocate(sed(i)%eflx(0:mt     ))      ; sed(i)%eflx = 0.0_sp
    allocate(sed(i)%dflx(0:mt     ))      ; sed(i)%dflx = 0.0_sp
    allocate(sed(i)%cdis(500      ))      ; sed(i)%cdis = 0.0_sp
    allocate(sed(i)%cflx(tmp2,kbm1))      ; sed(i)%cflx = 0.0_sp
    allocate(sed(i)%cobc(tmp3     ))      ; sed(i)%cobc = 0.0_sp
    allocate(sed(i)%depm(0:mt     ))      ; sed(i)%depm = 0.0_sp
    allocate(sed(i)%wset0(0:mt,kb  ))     ; sed(i)%wset0 = 0.0_sp
    sed(i)%wset0 = sed(i)%Wsed

    allocate(sed(i)%t_cd(0:mt     ))      ; sed(i)%t_cd = 0.0_sp
    allocate(sed(i)%t_ce(0:mt     ))      ; sed(i)%t_ce = 0.0_sp
    allocate(sed(i)%rate(0:mt     ))      ; sed(i)%rate = 0.0_sp

!J. Ge for tracer advection
    allocate(sed0(i)%conc(0:mt,kb  ))     ; sed0(i)%conc = 0.0
    allocate(sed2(i)%conc(0:mt,kb  ))     ; sed2(i)%conc = 0.0
!J. Ge for tracer advection

  end do

  !allocate suspended sediment arrays
  allocate(csed(0:mt,kb)) ; csed = 0.0

  !allocate bottom shear stress array 
  allocate(taub(0:mt)) ; taub = 0.0_sp

  if(SUSLOAD)then
    allocate ( SEDBED % ero_flux(0:mt,NST) )       ;SEDBED % ero_flux = 0.0_sp
    allocate ( SEDBED % settling_flux(0:mt,NST) )  ;SEDBED % settling_flux = 0.0_sp
  endif

  if(SED_BIOMASS)then
    allocate ( SEDBED % Dstp_max(0:mt,24) )        ;SEDBED % Dstp_max = 0.1_sp
    if(SEAGRASS_SINK)then
      allocate (SgrN (0:mt))                       ;SgrN = 0.0_sp
    end if
  endif

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Alloc_sedbed"
  RETURN

  END SUBROUTINE Alloc_sedbed



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
    do ised=1,nst
      if(SEDIMENT_PARAMETER_TYPE/=UNIFORM)then
        sum1 = sum1*(sed(ised)%t_ce(i))**sedbed%bed_frac(i,1,ised) 
      else
        sum1 = sum1*(sed(ised)%tau_ce)**sedbed%bed_frac(i,1,ised) 
      end if  
      sum2 = sum2*(sed(ised)%Sd50   )**sedbed%bed_frac(i,1,ised) 
      sum3 = sum3*(sed(ised)%wsed   )**sedbed%bed_frac(i,1,ised) 
      sum4 = sum4*(sed(ised)%Srho   )**sedbed%bed_frac(i,1,ised) 
    end do
    bottom(i,itauc) = sum1
    bottom(i,isd50) = sum2
    bottom(i,iwsed) = sum3
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
  write(*,*)'!  nbed                  :',nbed    
  write(*,*)'!  nst                   :',nst    
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
  do i=1,nst
    write(*,'(1X,A1,1X,A12,1X,A6,6(F7.4,1X),F8.2)')'!', &
              trim(sed(i)%sname),trim(sed(i)%stype(1:6)), &
              sed(i)%Sd50,sed(i)%Wsed, &
              sed(i)%tau_ce,sed(i)%tau_cd,sed(i)%erate,sed(i)%poros, &
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
  do i=1,nst
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
  do i=1,nst
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
  do i=1,nst
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

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Setup_Sed_OBC "

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Setup_Sed_OBC "

  End Subroutine Setup_Sed_OBC




!=======================================================================
!                                                                      !
!  This routine computes sediment bedload transport using the Meyer-   !
!  Peter and Muller (1948) formulation  for unidirectional flow and    !
!  Souksby and Damgaard (2005) algorithm that accounts for combined    !
!  effect of currents and waves.                                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Meyer-Peter, E. and R. Muller, 1948: Formulas for bedload transport !
!    In: Report on the 2nd Meeting International Association Hydraulic !
!    Research, Stockholm, Sweden, pp 39-64.                            !
!                                                                      !
!  Soulsby, R.L. and J.S. Damgaard, 2005: Bedload sediment transport   !
!    in coastal waters, Coastal Engineering, 52 (8), 673-689.          !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
  Subroutine calc_sed_bedload
    use control, only: msr ,par, mpi_fvcom_group
    use lims   , only: m,mt,kbm1, msrid
    use all_vars
# if defined (MULTIPROCESSOR)
    Use Mod_Par
# endif
    implicit none
    !
    !  Local variable declarations.
    !
    integer :: i, ised, j, k,cnt,jj,cc,ibed
    integer :: n1,n2,n3
    integer :: ia,ib,j1,j2

    integer, allocatable ::IAB_LIMIT(:)

    real(sp), parameter :: eps = 1.0E-14_sp

    real(sp) :: cff1, cff2, cff3, cff4, cff5
    real(sp) :: Dstp, bed_change
    real(sp), parameter :: one = 1.0_sp
    ! for BEDLOAD_MPM
    real(sp), dimension(0:mt) :: tau_w

    real(sp) :: sed_angle,max_slope
    real(sp) :: bedld, bedld_mass,poro
    real(sp) :: smgd, smgdr, osmgd, Umag
    real(sp) :: rhs_bed, Ra, phi, Clim,sed_slope
    real(sp) :: tau_x,tau_y,flowdir,bldx,bldy

    real(sp), allocatable :: cff(:)

    real(sp), dimension(0:mt) :: bld_dr
    real(sp), dimension(0:mt) :: dzdx
    real(sp), dimension(0:mt) :: dzdy
    real(sp), dimension(0:mt) :: a_slopex
    real(sp), dimension(0:mt) :: a_slopey

    real(sp), dimension(:),allocatable :: bflxe ! bedload scalar at cell center
    real(sp), dimension(:),allocatable :: bflxn ! bedload scalar at node point

    real(sp), dimension(0:nt) :: bustrcwmax_e
    real(sp), dimension(0:nt) :: bvstrcwmax_e

    real(sp), dimension(:), allocatable :: btx   !x-component of bedload transport
    real(sp), dimension(:), allocatable :: bty   !y-component of bedload transport

    real(sp), dimension(:), allocatable :: btxn  !x-component of bedload transport
    real(sp), dimension(:), allocatable :: btyn  !y-component of bedload transport

    real(sp), dimension(0:nt) :: dbdx_e !bed slope at cell face
    real(sp), dimension(0:nt) :: dbdy_e !bed slope at cell face

    real(sp), dimension(0:mt) :: dbdx_n !bed slope at node 
    real(sp), dimension(0:mt) :: dbdy_n !bed slope at node 

    ! for BEDLOAD_MPM 
    real(sp), dimension(0:mt) :: angleu
    real(sp), dimension(0:mt) :: anglev

    ! for BEDLOAD_SOULSBY
    real(sp) :: theta_mean, theta_wav, w_asym
    real(sp) :: theta_max, theta_max1, theta_max2
    real(sp) :: phi_x1, phi_x2, phi_x, phi_y
    real(sp) :: bedld_x, bedld_y, tau_cur, waven, wavec

    real(sp) :: dbdx,dbdy,beta,slopex,slopey

    real(sp), dimension(0:mt) :: phic
    real(sp), dimension(0:mt) :: phicw
    real(sp), dimension(0:mt) :: tau_wav
    real(sp), dimension(0:mt) :: tau_mean
    real(sp), parameter :: kdmax = 100.0_sp

    real(sp),allocatable,dimension(:,:)::temp
    real(sp),allocatable,dimension(:)::temp1

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start: calc_sed_bedload "

    allocate(btxn(0:mt));
    allocate(btyn(0:mt));
    allocate(btx(0:nt));
    allocate(bty(0:nt));
    allocate(bflxe(0:nt));
    allocate(bflxn(0:mt));

    !
    !-----------------------------------------------------------------------
    ! Compute maximum bottom stress for MPM bedload or suspended load.
    !-----------------------------------------------------------------------
    !
    if(BEDLOAD_MPM)then
       tau_w = taub /rho_water
     endif

    !
    !-----------------------------------------------------------------------
    !  Compute bedload sediment transport.
    !-----------------------------------------------------------------------
    !
    ! Compute some constant bed slope parameters.
    !
    sed_angle=DTAN(33.0_sp*pi/180.0_sp)
    !
    !  Compute angle between currents and waves (radians).
    !
    do i=1,m
# if defined (WAVE_CURRENT_INTERACTION)
       if(BEDLOAD_SOULSBY)then
          if(MB_BBL_USE.or.SG_BBL_USE.or.SSW_BBL_USE)then
             !
             ! Compute angle between currents and waves, measure CCW from current
             ! direction toward wave vector.
             !
             IF (bustrc(i).eq.0.0_sp) THEN
                phic(i)=0.5_sp*pi*SIGN(1.0_sp,bvstrc(i))
             ELSE
                phic(i)=ATAN2(bvstrc(i),bustrc(i))
             ENDIF
             phicw(i)=1.5_sp*pi-Dwave(i)-phic(i) !-angler(i)
             !
             ! Compute stress components at rho points.
             !
             tau_cur=SQRT(bustrc(i)**2+bvstrc(i)**2)
             tau_wav(i)=SQRT(bustrw(i)**2+bvstrw(j)**2)
             tau_mean(i)=tau_cur*(1.0_sp+1.2_sp*((tau_wav(i)/(tau_cur&
                  &+tau_wav(i)+eps))**3.2_sp))

          endif
       endif
# endif
       !
       if(BEDLOAD_MPM)then  ! && !defined BSTRESS_UPWIND
# if defined (WAVE_CURRENT_INTERACTION)
          if(MB_BBL_USE.or.SG_BBL_USE.or.SSW_BBL_USE)then
             cff1=bustrcwmax(i)
             cff2=bvstrcwmax(i)
          else
             cff1=wubot_n(i)
             cff2=wvbot_n(i)
          end if
# else
          cff1=wubot_n(i)
          cff2=wvbot_n(i)
# endif
          Umag=SQRT(cff1*cff1+cff2*cff2)+eps
          angleu(i)=cff1/Umag
          anglev(i)=cff2/Umag
       endif
    end do
    !
    DO ised=NCS+1,NST

       btxn  = 0.0_sp
       btyn  = 0.0_sp
       btx   = 0.0_sp
       bty   = 0.0_sp
       bflxe = 0.0_sp
       bflxn = 0.0_sp

       smgd=(sed(ised)%Srho/rho0-1.0_sp)*grav*sed(ised)%Sd50
       osmgd=1.0_sp/smgd
       smgdr=SQRT(smgd)*sed(ised)%Sd50*sed(ised)%Srho
       !
       do i=1,m
# if defined (WAVE_CURRENT_INTERACTION)
          if(BEDLOAD_SOULSBY)then
             if(MB_BBL_USE.or.SG_BBL_USE.or.SSW_BBL_USE)then
                !
                ! Compute wave asymmetry factor, based on Fredosoe and Deigaard.
                !
                Dstp=d(i)
                waven=2.0_sp*pi/(wlen(i)+eps)
                wavec=SQRT(grav/waven*tanh(waven*Dstp))
                cff4=MIN(waven*Dstp,kdmax)
                cff1=-0.1875_sp*wavec*(waven*Dstp)**2/(SINH(cff4))**4
                cff2=0.125_sp*grav*HSC1(i)**2/(wavec*Dstp+eps)
                cff3=pi*HSC1(i)/(Pwave_bot(i)*SINH(cff4)+eps)
                !
                ! Compute wave asymmetry factor, based on the note of Soulsby
                !
                cff1=MIN(0.375_sp*(HSC1(i)/Dstp)*                        &
                     &                    ((waven*Dstp)/(SINH(cff4))**3),0.15_sp)
                w_asym=2.0_sp*cff1/(1.0_sp+cff1**2.0_sp)
                !
                ! Compute nondimensional stresses.
                !
                theta_wav=tau_wav(i)*osmgd+eps
                theta_mean=tau_mean(i)*osmgd
                !
                cff1=theta_wav*(1.0_sp+w_asym)
                cff2=theta_wav*(1.0_sp-w_asym)
                theta_max1=SQRT((theta_mean+                              &
                     &                     cff1*COS(phicw(i)))**2+                  &
                     &                    (cff1*SIN(phicw(i)))**2)
                theta_max2=SQRT((theta_mean+                              &
                     &                     cff2*COS(phicw(i)+pi))**2+               &
                     &                    (cff2*SIN(phicw(i)+pi))**2)
                theta_max=MAX(theta_max1,theta_max2)
                !
                ! Motion initiation factor.
                !
                cff3=0.5_sp*(1.0_sp+SIGN(1.0_sp,                           &
                     &                               theta_max/sed(ised)%tau_ce-1.0_sp))
                !
                ! Calculate bed loads in direction of current and perpendicular
                ! direction.
                !
                phi_x1=12.0_sp*SQRT(theta_mean)*                           &
                     &           MAX((theta_mean-sed(ised)%tau_ce),0.0_sp)
                phi_x2=12.0_sp*(0.9534_sp+0.1907*COS(2.0_sp*phicw(i)))*    &
                     &           SQRT(theta_wav)*theta_mean+                         &
                     &           12.0_sp*(0.229_sp*w_asym*theta_wav**1.5_sp*         &
                     &                    COS(phicw(i)))
                !         phi_x=MAX(phi_x1,phi_x2) !  <- original
                IF (ABS(phi_x2).gt.phi_x1) THEN
                   phi_x=phi_x2
                ELSE
                   phi_x=phi_x1
                END IF
                bedld_x=phi_x*smgdr*cff3
                !
                cff5=theta_wav**1.5_sp+1.5_sp*(theta_mean**1.5_sp)
                phi_y=12.0_sp*0.1907_sp*theta_wav*theta_wav*               &
                     &          (theta_mean*SIN(2.0_sp*phicw(i))+1.2_sp*w_asym*      &
                     &          theta_wav*SIN(phicw(i)))/cff5*cff3
                bedld_y=phi_y*smgdr
                !
                ! Partition bedld into x- and y- directions.
                ! (btxn and btyn have dimensions of kg/m).
                !
                btxn(i)=(bedld_x*COS(phic(i))-bedld_y*SIN(phic(i)))*DTsed
                btyn(i)=(bedld_x*SIN(phic(i))+bedld_y*COS(phic(i)))*DTsed

             end if
          end if
# endif
          if(BEDLOAD_MPM)then
             !
             ! Magnitude of bed load at rho points. Meyer-Peter Muller formulation.
             ! (BEDLD has dimensions of kg m-1 s-1).
             !
             bedld=8.0_sp*(MAX((tau_w(i)*osmgd-0.047_sp),               &
                  &              0.0_sp)**1.5_sp)*smgdr
             !
             ! Partition bedld into xi and eta directions, still at rho points.
             ! (btxn and btyn have dimensions of kg/m ).
             !
             !         btxn(i)=angleu(i)*bedld*on_r(i)*DTsed ! original roms formular
             !         btyn(i)=anglev(i)*bedld*om_r(i)*DTsed ! original roms formular
             btxn(i) = angleu(i)*bedld*DTsed
             btyn(i) = anglev(i)*bedld*DTsed
          endif
          !
          ! Correct for along-direction slope. Limit slope to 0.9*sed angle.
          !
          !      cff1=0.5_sp*(1.0_sp+SIGN(1.0_sp,FX_r(i)))
          !      cff2=0.5_sp*(1.0_sp-SIGN(1.0_sp,FX_r(i)))
          !      cff3=0.5_sp*(1.0_sp+SIGN(1.0_sp,FE_r(i)))
          !      cff4=0.5_sp*(1.0_sp-SIGN(1.0_sp,FE_r(i)))
          !      if(SLOPE_NEMETH)then
          !            dzdx=(h(i+1,j)-h(i,j))/om_u(i+1,j)*cff1+                    &
          !     &           (h(i-1,j)-h(i,j))/om_u(i  ,j)*cff2
          !            dzdy=(h(i,j+1)-h(i,j))/on_v(i,j+1)*cff3+                    &
          !     &           (h(i,j-1)-h(i,j))/on_v(i  ,j)*cff4
          !            a_slopex=1.7_sp*dzdx
          !            a_slopey=1.7_sp*dzdy
          !
          ! Add contriubiton of bed slope to bed load transport fluxes.
          !
          !            FX_r(i,j)=FX_r(i,j)*(1.0_sp+a_slopex)
          !            FE_r(i,j)=FE_r(i,j)*(1.0_sp+a_slopey)
          !
          !   elseif(SLOPE_LESSER)then
          !            dzdx=MIN(((h(i+1,j)-h(i  ,j))/om_u(i+1,j)*cff1+             &
          !     &                (h(i  ,j)-h(i-1,j))/om_u(i  ,j)*cff2),0.52_sp)*   &
          !     &                SIGN(1.0_sp,FX_r(i,j))
          !            dzdy=MIN(((h(i,j+1)-h(i,j  ))/on_v(i,j+1)*cff3+             &
          !     &                (h(i,j  )-h(i,j-1))/on_v(i  ,j)*cff4),0.52_sp)*   &
          !     &                SIGN(1.0_sp,FE_r(i,j))
          !            cff=DATAN(dzdx)
          !            a_slopex=sed_angle/(COS(cff)*(sed_angle-dzdx))
          !            cff=DATAN(dzdy)
          !            a_slopey=sed_angle/(COS(cff)*(sed_angle-dzdy))
          !
          ! Add contriubiton of bed slope to bed load transport fluxes.
          !
          !            FX_r(i,j)=FX_r(i,j)*a_slopex
          !            FE_r(i,j)=FE_r(i,j)*a_slopey
          !    endif

       end do

       ! calculate the bed slope at cell center

       do i=1,n
          n1 = nv(i,1) ; n2 = nv(i,2) ; n3 = nv(i,3)

          dbdx_e(i) = awx(i,1)*h(n1)+awx(i,2)*h(n2)+awx(i,3)*h(n3)

          dbdy_e(i) = awy(i,1)*h(n1)+awy(i,2)*h(n2)+awy(i,3)*h(n3)

       end do

       ! interpolate slope into node 
       CALL E2N2D(dbdx_e,dbdx_n)
       CALL E2N2D(dbdy_e,dbdy_n)



       !----------CALCULATE DERIVATIVES OF DEPTH WITH X AND Y AT NODES----------------!
       !CALL DEPTH_GRADIENT

       if(SLOPE_NEMETH)then
          !
          ! Correct for along-direction slope. Limit slope to 0.9*sed angle.
          !
          bld_dr = sign(1.0_sp,btxn)
          dzdx=dbdx_n*bld_dr
          bld_dr = sign(1.0_sp,btyn)
          dzdy=dbdy_n*bld_dr
          a_slopex=1.7_sp*dzdx
          a_slopey=1.7_sp*dzdy
       
          ! Add contriubiton of bed slope to bed load transport fluxes.
          !
          btxn=btxn*(1.0_sp+a_slopex)
          btyn=btyn*(1.0_sp+a_slopey)
       elseif(SLOPE_LESSER)then
          sed_slope = tan(sed_angle*3.14159/180.)
          max_slope = 0.8*sed_slope
          bld_dr = sign(1.0_sp,btxn)
          dzdx=MIN(dbdx_n*bld_dr,max_slope)
          bld_dr = sign(1.0_sp,btyn)
          dzdy=MIN(dbdy_n*bld_dr,max_slope)
          allocate(cff(0:mt))
          cff=ATAN(dzdx)
          a_slopex=sed_angle/(COS(cff)*(sed_angle-dzdx))
          cff=ATAN(dzdy)
          a_slopey=sed_angle/(COS(cff)*(sed_angle-dzdy))
          deallocate(cff)
          !
          ! Add contriubiton of bed slope to bed load transport fluxes.
          !
          btxn=btxn*a_slopex
          btyn=btyn*a_slopey
       end if


       if(SED_MORPH)then
          !
          ! Apply morphology factor.
          !
          btxn=btxn*sed(ised)%morph_fac
          btyn=btyn*sed(ised)%morph_fac
       endif
       !
       ! Apply bedload transport rate coefficient. Also limit
       ! bedload to the fraction of each sediment class.
       !

       btxn=btxn*bedload_coeff*sedbed%bed_frac(:,1,ised)
       btyn=btyn*bedload_coeff*sedbed%bed_frac(:,1,ised)

       !
       ! Limit bed load to available bed mass.
       !
       ! ROMS does this, but doesn't make sense, mass is linked to divergence of  
       ! transport, not horizontal fluxes , check with J. Warner

       !     bedld_mass=ABS(FX_r(i,j))+ABS(btyn(i,j))+eps
       !     FX_r=MIN(ABS(FX_r),sedbed%bed_mass(1:m,1,nstp,ised)*om_r*on_r*ABS(FX_r)/bedld_mass)*SIGN(1.0_sp,FX_r)
       !     btyn=MIN(ABS(btyn),sedbed%bed_mass(1:m,1,nstp,ised)*om_r*on_r*ABS(btyn)/bedld_mass)*SIGN(1.0_sp,btyn)


       !exchange btx/bty across interprocessor boundaries 
# if defined (MULTIPROCESSOR)
       if(par)then
          call aexchange(nc,myid,nprocs,btxn,btyn)
       endif
# endif

       call n2e2d(btxn,btx)
       call n2e2d(btyn,bty)

       !add bed slope contribution to our horizontal fluxes
       !Eq 35, Warner etal, Comp & Geo, 2008 => qbl_slope = slopex,slopey 
       !minimize the maximum downslope, but do not truncate upslope
       !upslope is limited by Lesser formula
       !range of slopefactor (slopex)  = [.55,5.64] for ..
       !..  upslope(dhdx<0)=-infty -> downslope(dhdx>0) = infty
       !dbdx > 0 is a downslope (see bld_dr to enforce this)

       !sed_slope = tan(sed_angle*3.14159/180.)
       !max_slope = 0.8*sed_slope
       !do i=1,n
       !   n1 = nv(i,1) ; n2 = nv(i,2) ; n3 = nv(i,3)

       !   bld_dr = sign(1.0_sp,btx(i))
       !   dbdx   = min( (awx(i,1)*h(n1)+awx(i,2)*h(n2)+awx(i,3)*h(n3))*bld_dr ,max_slope)
       !   beta   = atan(dbdx)
       !   slopex = sed_slope/( (sed_slope-dbdx)*cos(beta)) 
       !   btx(i) = btx(i)*slopex

       !   bld_dr = sign(1.0_sp,bty(i))
       !   dbdy   = min( (awy(i,1)*h(n1)+awy(i,2)*h(n2)+awy(i,3)*h(n3))*bld_dr ,max_slope)
       !   beta   = atan(dbdy)
       !   slopey = sed_slope/( (sed_slope-dbdy)*cos(beta)) 
       !   bty(i) = bty(i)*slopey

       !end do



# if defined (WAVE_CURRENT_INTERACTION)
       if(MB_BBL_USE.or.SG_BBL_USE.or.SSW_BBL_USE)then
          call n2e2d(bustrcwmax,bustrcwmax_e)
          call n2e2d(bvstrcwmax,bvstrcwmax_e)
       end if
# endif

       !calculate the flux on the elements using first order upwind formulation 
       do i=1,ne
          ia=iec(i,1)
          ib=iec(i,2)
          j1=ienode(i,1)
          j2=ienode(i,2)
          !bed stress in x and y averaged to an edge , used for upwind selection
# if defined (WAVE_CURRENT_INTERACTION)
          if(MB_BBL_USE.or.SG_BBL_USE.or.SSW_BBL_USE)then
             tau_x = -0.5*(bustrcwmax_e(ia)+bustrcwmax_e(ib))
             tau_y = -0.5*(bvstrcwmax_e(ia)+bvstrcwmax_e(ib))
          else
             tau_x = 0.5*(wubot(ia)+wubot(ib))  
             tau_y = 0.5*(wvbot(ia)+wvbot(ib))
          endif
# else
          tau_x = 0.5*(wubot(ia)+wubot(ib))  
          tau_y = 0.5*(wvbot(ia)+wvbot(ib))  
# endif
          !dot the bed stress vector with length-weighted edge normal 
          !note bed stress opposes flow dir so we need neg sign
          !the vector (-dltyc,dltxc) points from ia to ib
          flowdir  = -(-tau_x*dltyc(i) + tau_y*dltxc(i)) ! -(tau dot dL)
          !upwind select the x-dir and y-dir bedload transport 
          bldx  = ((one-sign(one,flowdir))*btx(ib)+(one+sign(one,flowdir))*btx(ia))*0.5_sp 
          bldy  = ((one-sign(one,flowdir))*bty(ib)+(one+sign(one,flowdir))*bty(ia))*0.5_sp 
          !accum contributions to bedload transport divergence calculations at elements ia/ib
          !on either side of edge i 
          ! => bedload in kg
          bflxe(ia) = bflxe(ia) - bldx*dltyc(i) + bldy*dltxc(i)
          bflxe(ib) = bflxe(ib) + bldx*dltyc(i) - bldy*dltxc(i)
       end do


       !limit fluxes to prevent bottom from breaking thru water surface.
       allocate(cff(0:nt)); cff = 0.0_sp
       do i=1,n
          ! Total thickness available.
          Dstp=MAX(d1(i),0.0_sp)
          ! Thickness change that wants to occur.
          n1 = nv(i,1) ; n2 = nv(i,2) ; n3 = nv(i,3)
          poro=(bed(n1,1,iporo)+bed(n2,1,iporo)+bed(n3,1,iporo))/3.0
          bed_change= bflxe(i)/(sed(ised)%Srho*(1.0_sp-poro))
          ! calculate the coefficients to make sure the change breaking
          ! the water surface
          cff(i)=MAX(bed_change-Dstp,0.0_sp)
          cff(i)=cff(i)/ABS(bed_change+eps)
          ! Limit that change to be less than available.
          bflxe(i) = bflxe(i)*(1-cff(i))
       end do

       ! re-Limit the UV transport components that accumulated to bedload 
       allocate(iab_limit(0:nt)); iab_limit = 0 ! flag to avoid multiple cutting
       do i=1,ne
          ia=iec(i,1)
          ib=iec(i,2)
          if(iab_limit(ia)==0)then
             btx(ia) = btx(ia)*(1-cff(ia))
             bty(ia) = bty(ia)*(1-cff(ia))
             iab_limit(ia)=1
          endif
          if(iab_limit(ib)==0)then
             btx(ib) = btx(ib)*(1-cff(ib))
             bty(ib) = bty(ib)*(1-cff(ib))
             iab_limit(ib)=1
          endif
       end do
       deallocate(iab_limit)
       deallocate(cff)


       !divide by element area to complete divergence => bedload in kg/m^2
       !if there was a divergence of bedload transport, this should be a positive number
       do i=1,n
          bflxe(i) = bflxe(i)/art(i)
       end do

       !exchange bflxe across interprocessor boundaries
#  if defined (MULTIPROCESSOR)
       if(par)then
          call aexchange(ec,myid,nprocs,bflxe)
          call aexchange(ec,myid,nprocs,btx,bty)
       endif
#  endif


       ! interpolate element-based blfx to the nodes
       call e2n2d(bflxe,bflxn)
       call e2n2d(btx,btxn)   
       call e2n2d(bty,btyn)  


#  if defined (MULTIPROCESSOR)
       if(par)then
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,1,myid,nprocs,bflxn)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,1,myid,nprocs,btxn)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,1,myid,nprocs,btyn)
          call aexchange(nc,myid,nprocs,bflxn)
          call aexchange(nc,myid,nprocs,btxn,btyn)
       endif
#  endif


       ! evaluate change in bed properties.
       sedbed%bed_mass(1:m,1,nnew,ised)=MAX(sedbed%bed_mass(1:m,1,nstp,ised)&
            &-bflxn(1:m),0.0_sp)
       if(.not.SUSLOAD)then
          DO k=2,Nbed
             sedbed%bed_mass(1:m,k,nnew,ised)=sedbed%bed_mass(1:m,k,nstp,ised)
          END DO
       endif


       bed(1:m,1,ithck)=MAX(bed(1:m,1,ithck)- bflxn(1:m)/(sed(ised)%Srho*(1.0_sp&
            &-bed(1:m,1,iporo))),0.0_sp)

       !
       !-----------------------------------------------------------------------
       !  Output bedload fluxes.
       !-----------------------------------------------------------------------
       !

       sedbed%bedldu(1:m,ised)=btxn(1:m)/Dtsed
       sedbed%bedldv(1:m,ised)=btyn(1:m)/Dtsed

    END DO


    !
    !  Need to update bed mass for the cohesive sediment types, becasue 
    !  they did not partake in the bedload transport.
    !

    DO ised=1,NCS
       sedbed%bed_mass(1:m,1,nnew,ised)=sedbed%bed_mass(1:m,1,nstp,ised)
       if(.not.SUSLOAD)then
          DO k=2,Nbed
             sedbed%bed_mass(1:m,k,nnew,ised)=sedbed%bed_mass(1:m,k,nstp,ised)
          END DO
       endif
    END DO

    !
    !  Update mean surface properties.
    !  Sd50 must be positive definite, due to BBL routines.
    !  Srho must be >1000, due to (s-1) in BBL routines.
    !
    DO i=1,m
       cff3=0.0_sp
       DO ised=1,NST
          cff3=cff3+sedbed%bed_mass(i,1,nnew,ised)
       END DO
       IF (cff3.eq.0.0_sp) THEN
          cff3=eps
       END IF
       DO ised=1,NST
          sedbed%bed_frac(i,1,ised)=sedbed%bed_mass(i,1,nnew,ised)/cff3
       END DO
    End DO



    Call Update_Bottom_Properties

    deallocate(btx,bty,btxn,btyn,bflxe,bflxn)

    RETURN


    if(dbg_set(dbg_sbr)) write(ipt,*) "End: calc_sed_bedload "

  End Subroutine calc_sed_bedload






!=================================================================== ===
!                                                                      !
!  This routine computes the vertical settling (sinking) of suspended  !
!  sediment via a semi-Lagrangian advective flux algorithm. It uses a  !
!  parabolic,  vertical reconstructuion of the suspended  sediment in  !
!  the water column with PPT/WENO constraints to avoid oscillations.   !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Colella, P. and P. Woodward, 1984: The piecewise parabolic method   !
!    (PPM) for gas-dynamical simulations, J. Comp. Phys., 54, 174-201. !
!                                                                      !
!  Liu, X.D., S. Osher, and T. Chan, 1994: Weighted essentially        !
!    nonoscillatory shemes, J. Comp. Phys., 115, 200-212.              !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
  Subroutine calc_sed_settling

  implicit none
!
!  Local variable declarations.
!
  integer :: i, indx, ised, j, k, ks

  real(sp) :: cff, cu, cffL, cffR, dltL, dltR

  integer, dimension(0:kb) :: ksource

  real(sp), dimension(0:kb) :: FC

  real(sp), dimension(0:kb) :: Hz_inv
  real(sp), dimension(0:kb) :: Hz_inv2
  real(sp), dimension(0:kb) :: Hz_inv3
  real(sp), dimension(0:kb) :: qc
  real(sp), dimension(0:kb) :: qR
  real(sp), dimension(0:kb) :: qL
  real(sp), dimension(0:kb) :: WR
  real(sp), dimension(0:kb) :: WL

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: calc_sed_settling "

!
!-----------------------------------------------------------------------
!  Add sediment vertical sinking (settling) term.
!-----------------------------------------------------------------------
!
!  Compute inverse thicknesses to avoid repeated divisions.
!

!  DO ised=1,NST
!    print*,'before settling:',ised,maxval(sed(ised)%conc),minval(sed(ised)%conc)
!  end do

  do i=1,m
    DO k=1,kbm1
       Hz_inv(k)=1.0_sp/(dt(i)*dz(i,k))
    END DO
    DO k=2,kbm1
       Hz_inv2(k)=1.0_sp/(dt(i)*dz(i,k)+dt(i)*dz(i,k-1))
    END DO
    DO k=2,kbm1-1
       Hz_inv3(k)=1.0_sp/(dt(i)*(dz(i,k-1)+dz(i,k)+dz(i,k+1))) 
    END DO
!
!  Copy concentration of suspended sediment into scratch array "qc"
!  (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for sediment
!  concentration.
!

    SED_LOOP: DO ised=1,NST
    DO k=1,kbm1
       qc(k)=sed(ised)%conc(i,k)
    END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking of suspended sediment.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!
    DO k=2,kbm1
        FC(k)=(qc(k-1)-qc(k))*Hz_inv2(k)
    END DO

    DO k=kbm1-1,2,-1
       dltR=(dt(i)*dz(i,k))*FC(k)
       dltL=(dt(i)*dz(i,k))*FC(k+1)
       cff=dt(i)*(dz(i,k-1)+2.0_sp*dz(i,k)+dz(i,k+1))
       cffR=cff*FC(k)
       cffL=cff*FC(k+1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
       IF ((dltR*dltL).le.0.0_sp) THEN
         dltR=0.0_sp
         dltL=0.0_sp
       ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
         dltR=cffL
       ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
         dltL=cffR
       END IF
!
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
       cff=(dltR-dltL)*Hz_inv3(k)
       dltR=dltR-cff*(dt(i)*dz(i,k-1)) 
       dltL=dltL+cff*(dt(i)*dz(i,k+1)) 
       qR(k)=qc(k)+dltR
       qL(k)=qc(k)-dltL
       WR(k)=(2.0_sp*dltR-dltL)**2
       WL(k)=(dltR-2.0_sp*dltL)**2
    END DO
    cff=1.0E-14_sp
    DO k=kbm1-1,3,-1
       dltL=MAX(cff,WL(k  ))
       dltR=MAX(cff,WR(k-1))
       qR(k)=(dltR*qR(k)+dltL*qL(k-1))/(dltR+dltL)
       qL(k-1)=qR(k)
    END DO
    FC(1)=0.0_sp              ! no-flux boundary condition
    if(LINEAR_CONTINUATION)then
      qL(1)=qR(2)
      qR(1)=2.0_sp*qc(1)-qL(1)
    elseif(NEUMANN)then
      qL(1)=qR(2)
      qR(1)=1.5_sp*qc(1)-0.5_sp*qL(1)
    else
      qR(1)=qc(1)         ! default strictly monotonic
      qL(1)=qc(1)         ! conditions
      qR(2)=qc(1)
    endif
    if(LINEAR_CONTINUATION)then
      qR(kbm1)=qL(kbm1-1)
      qL(kbm1)=2.0_sp*qc(kbm1)-qR(kbm1)
    elseif(NEUMANN)then
      qR(kbm1)=qL(kbm1-1)
      qL(kbm1)=1.5_sp*qc(kbm1)-0.5_sp*qR(kbm1)
    else
      qL(kbm1-1)=qc(kbm1)                 ! bottom grid boxes are
      qR(kbm1)=qc(kbm1)                 ! re-assumed to be
      qL(kbm1)=qc(kbm1)                 ! piecewise constant.
    endif
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
    DO k=kbm1,1,-1
       dltR=qR(k)-qc(k)
       dltL=qc(k)-qL(k)
       cffR=2.0_sp*dltR
       cffL=2.0_sp*dltL
       IF ((dltR*dltL).lt.0.0_sp) THEN
         dltR=0.0_sp
         dltL=0.0_sp
       ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
         dltR=cffL
       ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
         dltL=cffR
       END IF
       qR(k)=qc(k)+dltR
       qL(k)=qc(k)-dltL
    END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by kbm1).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
    cff=DTsed*ABS(sed(ised)%Wsed)
    DO k=kbm1,1,-1
       FC(k+1)=0.0_sp
       WL(k)=dt(i)*zz(i,k+1)+cff
       WR(k)=dt(i)*dz(i,k)*qc(k)
       ksource(k)=k
    END DO
    DO k=kbm1,1,-1
       DO ks=k,2,-1
         IF (WL(k).gt.dt(i)*zz(i,ks)) THEN
            ksource(k)=ks-1
            FC(k+1)=FC(k+1)+WR(ks)
         END IF
       END DO
    END DO
!
!  Finalize computation of flux: add fractional part.
!
    DO k=kbm1,1,-1
       ks=ksource(k)
       cu=MIN(1.0_sp,(WL(k)-dt(i)*zz(i,ks+1))*Hz_inv(ks))
       FC(k+1)=FC(k+1)+                                      &
     &           dt(i)*dz(i,ks)*cu*                                  &
     &           (qL(ks)+                                      &
     &            cu*(0.5_sp*(qR(ks)-qL(ks))-                &
     &                (1.5_sp-cu)*                               &
     &              (qR(ks)+qL(ks)-2.0_sp*qc(ks))))
     END DO
     DO k=kbm1,1,-1
        sed(ised)%conc(i,k)=sed(ised)%conc(i,k)+(FC(k)-FC(k+1))
     END DO
     sedbed%settling_flux(i,ised)=FC(kb)

    END DO SED_LOOP

  END DO

!  DO ised=1,NST
!    print*,'after settling:',ised,maxval(sed(ised)%conc),minval(sed(ised)%conc)
!  end do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: calc_sed_settling "

  RETURN

  End Subroutine calc_sed_settling



!=======================================================================
!                                                                      !
!  This computes sediment bed and water column exchanges: deposition,  !
!  resuspension, and erosion.                                          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
  Subroutine calc_sed_fluxes

    implicit none
    !
    !  Local variable declarations.
    !
    integer :: Ksed, i, indx, ised, j, k, ks
    integer :: bnew, cnt,jj,cc

    real(sp) :: cff, cff1, cff2, cff3, cff4, cff5
    real(sp), dimension(0:mt) :: tau_w

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start:  calc_sed_fluxes"

    if(BEDLOAD)then
      bnew=nnew
    else
      bnew=nstp
    endif

    !
    !-----------------------------------------------------------------------
    !  Compute sediment deposition, resuspension, and erosion.
    !-----------------------------------------------------------------------
    !
    if(BEDLOAD_MPM.or.SUSLOAD)then
      tau_w = taub
    endif

# if !defined (FLUID_MUD)
    !
    !-----------------------------------------------------------------------
    !  Sediment deposition and resuspension near the bottom.
    !-----------------------------------------------------------------------
    !
    !  The deposition and resuspension of sediment on the bottom "bed"
    !  is due to precepitation settling_flux, already computed, and the
    !  resuspension (erosion, hence called ero_flux). The resuspension is
    !  applied to the bottom-most grid box value qc(:,1) so the total mass
    !  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
    !  after both fluxes are applied.
    !
    NODE_LOOP : DO i=1,m
# if defined (WET_DRY)
       if(iswetn(i)/=1)cycle  ! no calculation when dry node
# endif
       SED_LOOP: DO ised=1,NST
          !
          !  Calculate critical shear stress in Pa
          !
          if(COHESIVE_BED)then
!             cff =1.0_sp/bed(i,1,ibtcr)
            cff=1.0_sp/sed(ised)%tau_ce
          elseif(MIXED_BED)then
             cff = MAX(bottom(i,idprp)*bed(i,1,ibtcr)+          &
                  &        (1.0_sp-bottom(i,idprp))*sed(ised)%tau_ce,           &
                  &            sed(ised)%tau_ce)
             cff=1.0_sp/cff
          else
             cff=1.0_sp/sed(ised)%tau_ce
          endif
          !
          !  Compute erosion, ero_flux (kg/m2).
          !
          cff1=(1.0_sp-bed(i,1,iporo))*sedbed%bed_frac(i,1,ised)
          cff2=DTsed*sed(ised)%Erate*cff1
          cff3=sed(ised)%Srho*cff1
          cff4=sedbed%bed_mass(i,1,bnew,ised)
          cff5=sedbed%settling_flux(i,ised) !CRS
          if(SED_TAU_CD_CONST.or.SED_TAU_CD_LIN)then
             if (tau_w(i).GT.sed(ised)%tau_cd) then
                cff5=0.0_sp
             endif
          endif
          if(SED_TAU_CD_LIN)then
             if (tau_w(i).LE.sed(ised)%tau_cd.AND.                      &
                  &                       sed(ised)%tau_cd.GT.0.0_sp) then
                cff5=(1.0_sp-tau_w(i)/sed(ised)%tau_cd)*cff5
             endif
          endif

          sedbed%ero_flux(i,ised)=                                         &
               &         MIN(MAX(0.0_sp,cff2*(cff*tau_w(i)-1.0_sp)),      &
               &             MIN(cff3*bottom(i,iactv),cff4)+              &
               &                cff5)

          !
          !  Update global tracer variables (mT units) for erosive flux.
          !
          sed(ised)%conc(i,kbm1)=sed(ised)%conc(i,kbm1)+sedbed%ero_flux(i,ised)
          if (SED_TAU_CD_CONST .or. SED_TAU_CD_LIN)then
             sed(ised)%conc(i,kbm1)=sed(ised)%conc(i,kbm1) +               &
                  &                       (sedbed%settling_flux(i,ised)-cff5)
             sedbed%settling_flux(i,ised)=cff5
          endif

       END DO SED_LOOP
    END DO NODE_LOOP

# endif

    ! J. Ge for fluid mud layer
# if defined (FLUID_MUD)
    Settling = 0.0_SP
    do ised=1,nst
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

    if(dbg_set(dbg_sbr)) write(ipt,*) "End: calc_sed_fluxes"

    return

  End Subroutine calc_sed_fluxes


!=======================================================================
!                                                                      !
!  This computes sediment biomass due to vegation growth               !
!                                                                      !
!  References:                                                         !
!                                                                      !
!                                                                      !
!=======================================================================
  Subroutine calc_sed_biomass

  implicit none
!
!  Local variable declarations.
!
  integer :: i, j, k, ised
  integer :: sstp, nbio_steps
  real(sp) :: cff, Dstp
  real(sp) :: sgr_kgmmol

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: calc_sed_biomass "

!
!-----------------------------------------------------------------------
!  Compute 
!-----------------------------------------------------------------------
!
!  Compute number of model steps for each hour.
!
  nbio_steps=MAX(1,INT(3600.0_sp/DTsed))
!
!  Compute number of hourly values we need to save.
!  If we want a 1 day avg, then need 24 values.
!
  NODE_LOOP : DO i=1,m
# if defined (WET_DRY)
    if(iswetn(i)/=1)cycle  ! no calculation when dry node
# endif
!
!  Only update the max depth once per hour.
!
    IF (MOD(iint,nbio_steps).eq.0) THEN
!
!  Determine the index for placement of new value.
!
!      sstp=1+MOD(iint-ntstart,24)
       sstp=1+MOD(iint,24)
!
!  Save instantaneous depth at this instance and recompute max daily depth.
!
       Dstp=d(i)
       sedbed%Dstp_max(i,sstp)=Dstp
       cff=0.0_sp
       DO k=1,nTbiom                        
         cff=MAX(cff,sedbed%Dstp_max(i,k))
       END DO
       bottom(i,imaxD)=cff
    END IF
!
!  Seagrass as a bottom property
!
    if(SEAGRASS_BOTTOM)then
      if(SEAGRASS_SINK)then
!         sgr_diam=0.01_sp
!         SgrN has units of millimole_nitrogen meter-3
!         change moles to kg (2.8e-5 kg/millimole N2) switch to right formula
!         cylinder height is mass over (density*pi*r*r)
!         shoot height is cylinder height / shoot density
          sgr_kgmmol = 2.8e-5_sp
!         sgr_density = 500.0_sp
!         sgr_Hthres = 1.25_sp
          cff = SgrN(i)*sgr_kgmmol/                                   &
     &                    (0.25_sp*sgr_diam*sgr_diam*pi*sgr_density)
!
          bottom(i,isgrH) = cff/bottom(i,isgrD)
!
          IF (bottom(i,isgrH).gt.sgr_Hthres) THEN
               bottom(i,isgrD)=bottom(i,isgrD)*                     &
     &                    bottom(i,isgrH)/sgr_Hthres
               bottom(i,isgrH) = sgr_Hthres
          END IF
      else
          bottom(i,isgrH)=1.25_sp
          bottom(i,isgrD)=400.0_sp
      endif
    endif
!
!  Update settling flux for depositing bio mass.
!
!
!  Require (for now) that the first sed class be the new biomass.
!
    ised=1
    cff=0.0_sp                               ! remove this line
!    cff= funct( bottom(i,imaxD), DTsed)   ! need real eq. in kg/m^2


    sedbed%settling_flux(i,ised)=sedbed%settling_flux(i,ised)+cff
  END DO NODE_LOOP


  if(dbg_set(dbg_sbr)) write(ipt,*) "End: calc_sed_biomass "

  return

  End Subroutine calc_sed_biomass



!=======================================================================
!                                                                      !
!  This routine computes sediment bed layer stratigraphy.              !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
  Subroutine calc_sed_bed_cohesive

    implicit none

    !
    !  Local variable declarations.
    !
    integer :: Ksed, i, ised,ibed, j, k, ks
    integer :: bnew, nnn, cnt,cc,jj

    real(sp), parameter :: eps = 1.0E-14_sp

    real(sp) :: cff, cff1, cff2, cff3,cff4
    real(sp) :: thck_avail, thck_to_add

    real(sp), dimension(NST) :: nlysm

    real(sp), dimension(NST) :: dep_mass
    real(sp), dimension(0:mt) :: tau_w

    real(sp) :: pcoh 

    real(sp) :: alpha, bzactv, frt, tcb_top, tcb_bot
    real(sp), dimension(Nbed) :: tcr

    real(sp), dimension(Nbed) :: bmz

    real(sp), dimension(NCS) :: masseq

    real(sp),allocatable,dimension(:,:)::temp

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start:  calc_sed_bed_cohesive"

    if(BEDLOAD)then
       bnew=nnew
    else
       bnew=nstp
    endif

    ! KLUDGE alert
    ! minlayer_thick = 0.0005
    minlayer_thick = newlayer_thick
    !
    !-----------------------------------------------------------------------
    ! Compute sediment bed layer stratigraphy.
    !-----------------------------------------------------------------------
    !
    if(BEDLOAD_MPM.or.SUSLOAD)then
       tau_w = taub
    endif

    !
    !-----------------------------------------------------------------------
    !  Update bed properties according to ero_flux and dep_flux.
    !-----------------------------------------------------------------------
    !
    if(SUSLOAD)then
       I_LOOP : DO i=1,m
          !
          !  The deposition and resuspension of sediment on the bottom "bed"
          !  is due to precipitation flux FC(:,0), already computed, and the
          !  resuspension (erosion, hence called ero_flux). The resuspension is
          !  applied to the bottom-most grid box value qc(:,1) so the total mass
          !  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
          !  after both fluxes are applied.
          !
          SED_LOOP: DO ised=1,NST  
             dep_mass(ised)=0.0_sp
             if(SED_MORPH)then
                ! Apply morphology factor.
                sedbed%ero_flux(i,ised)=sedbed%ero_flux(i,ised)*sed(ised)%morph_fac
                sedbed%settling_flux(i,ised)=sedbed%settling_flux(i,ised)*            &
                     &  sed(ised)%morph_fac
             endif

             !  Update bed mass arrays.
             sedbed%bed_mass(i,1,nnew,ised)=MAX(sedbed%bed_mass(i,1,bnew,ised)-     &
                  &    (sedbed%ero_flux(i,ised)-             &
                  &    sedbed%settling_flux(i,ised)),        &
                  &        0.0_sp)
             DO k=2,Nbed
                sedbed%bed_mass(i,k,nnew,ised)=sedbed%bed_mass(i,k,nstp,ised)
             END DO
          END DO SED_LOOP

          cff3=0.0_sp
          DO ised=1,NST
             cff3=cff3+sedbed%bed_mass(i,1,nnew,ised)
          END DO

          IF (cff3.eq.0.0_sp) THEN 
             cff3=eps 
          END IF
          bed(i,1,ithck)=0.0_sp
          DO ised=1,NST
             sedbed%bed_frac(i,1,ised)=sedbed%bed_mass(i,1,nnew,ised)/cff3
             bed(i,1,ithck)=MAX(bed(i,1,ithck)+        &
                  &   sedbed%bed_mass(i,1,nnew,ised)/       &
                  &   (sed(ised)%Srho*                      &
                  &   (1.0_sp-bed(i,1,iporo))),0.0_sp)
          END DO

       END DO I_LOOP
    endif !/* SUSPLOAD section */

    !
    !-----------------------------------------------------------------------
    !  At this point, all deposition or erosion is complete, and
    !  has been added/subtracted to top layer. Thickness has NOT been corrected.
    !-----------------------------------------------------------------------
    !
    I_LOOP_CB : DO i=1,m
       !       Calculate active layer thickness, bottom(i,iactv).
       !       (trunk version allows this to be zero...this has minimum of 6*D50)
       bottom(i,iactv)=MAX(0.0_sp,                                 &
            &   0.007_sp*                               &
            &   (tau_w(i)-bottom(i,itauc)))+   &
            &   6.0_sp*bottom(i,isd50)
       if(COHESIVE_BED)then 
          bottom(i,iactv)=MAX(0.0_sp,                                 &
               &     0.007_sp*                               &
               &     (tau_w(i)-bed(i,1,ibtcr)))+    &
               &     6.0_sp*bottom(i,isd50)
       endif
       if(MIXED_BED)then
          cff1=  MAX(0.0_sp,                                 &
               &     0.007_sp*                               &
               &     (tau_w(i)-bed(i,1,ibtcr)))+    &
               &     6.0_sp*bottom(i,isd50)
          cff2=  MAX(0.0_sp,                                 &
               &     0.007_sp*                               &
               &     (tau_w(i)-bottom(i,itauc)))+   &
               &     6.0_sp*bottom(i,isd50)
          bottom(i,iactv)=MAX(cff1,cff2)
       endif
       if(SED_MORPH)then
          ! Apply morphology factor.
          bottom(i,iactv)=MAX(bottom(i,iactv)*sed(ised)%morph_fac,      &
               &          bottom(i,iactv))
       endif
       if (COHESIVE_BED .or. MIXED_BED)then
          !  Find first layer with tc > tb
          !  Remember, the critical stresses apply to the TOP of each layer
          Ksed = 1
          bzactv = 0.0_sp
          frt = 0.0_sp ! CRS
          cff1 = tau_w(i)
          tcb_top = bed(i,1,ibtcr)
          tcb_bot = bed(i,2,ibtcr)
          if (MIXED_BED)then
             !         Calc. tau crit for bottom of next layer
             !         Update cohesive property of seds in top layer
             cff3 = 0.0_sp
             cff4 = 1.0_sp
             DO ised=1,NCS
                cff3=cff3+sedbed%bed_frac(i,1,ised)
                cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,1,ised)
             END DO
             DO ised=NCS+1,NST
                cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,1,ised)
             ENDDO
             pcoh=min(max((cff3-transN)/(transC-transN)        &
                  &      ,0.0_sp),1.0_sp)
             tcb_top = pcoh*bed(i,1,ibtcr)+(pcoh-1.0_sp)*cff4   !BF
          endif
          IF(cff1 .GT. tcb_top)THEN
             !            Calculate tcb_temp for next layer

             tcb_bot = bed(i,Ksed+1,ibtcr)
             if (MIXED_BED)then
                ! Recalculate cohesive fraction and mean grain tau crit
                ! Note that Ksed is used for grain props, and Ksed+1 for
                ! bed tau crit at bottom of layer Ksed
                cff3 = 0.0_sp
                cff4 = 1.0_sp
                DO ised=1,NCS
                   cff3=cff3+sedbed%bed_frac(i,Ksed,ised)
                   cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,Ksed,ised)
                END DO
                DO ised=NCS+1,NST
                   cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,Ksed,ised)
                ENDDO
                ! Calculate cohesive behavior and blended tau crit
                pcoh=min(max((cff3-transN)/(transC-transN),    &
                     &            0.0_sp),1.0_sp)
                tcb_bot =pcoh*bed(i,Ksed+1,ibtcr)+(pcoh-1.0_sp)*cff4  !BF
             endif
             DO WHILE ( (Ksed.LE.(Nbed-1)) .AND.                        &
                  &            (cff1 .GT. tcb_bot))
                !ALA   Instead of adding entire 2nd layer, add just what is needed
                !ALA!  Add entire layer
                !ALA   bzactv = bzactv + bed(i,Ksed,ithck)
                !      Add thickness equivalent to difference between cff1 and tcb_bot
                IF (cff1.GT.bed(i,Ksed+1,ibtcr)) THEN
                   bzactv = bzactv + bed(i,Ksed,ithck)
                   tcb_top = tcb_bot
                   tcb_bot = bed(i,Ksed+1,ibtcr)
                ELSE
                   bzactv = bzactv + MIN(bed(i,Ksed,ithck),           &
                        &    (MAX(0.0_sp,0.007_sp*(cff1-tcb_bot))))
                   tcb_top = cff1
                   tcb_bot = bed(i,Ksed+1,ibtcr)   
                ENDIF
                if(MIXED_BED)then
                   ! Recalculate cohesive fraction and mean grain tau crit
                   cff3 = 0.0_sp
                   cff4 = 1.0_sp
                   DO ised=1,NCS
                      cff3=cff3+sedbed%bed_frac(i,Ksed,ised)
                      cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,Ksed,ised)
                   END DO
                   DO ised=NCS+1,NST
                      cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,Ksed,ised)
                   ENDDO
                   !               Calculate cohesive behavior and blended tau crit
                   pcoh=min(max((cff3-transN)/(transC-transN), &
                        &  0.0_sp),1.0_sp)
                   tcb_bot = pcoh*bed(i,Ksed+1,ibtcr)+(pcoh-1.0_sp)*cff4    !BF
                endif
                Ksed = Ksed+1
             ENDDO
             frt = MAX(0.0_sp,(cff1-tcb_top) ) /                        &
                  &   MAX( eps,                                             &
                  &   ( tcb_bot-tcb_top ))
          ENDIF
          bzactv = bzactv+frt*bed(i,Ksed,ithck)
          bzactv = MAX( bzactv, 6.0_sp*bottom(i,isd50) ) !CRS
          bottom(i,iactv)=min(bottom(i,iactv),bzactv) !?CRS
          !          bottom(i,iactv)=max(bottom(i,iactv),bzactv) !?CRS
       endif !/* defined COHESIVE_BED || defined MIXED_BED */

    END DO I_LOOP_CB
 

    J_LOOP2 : DO i=1,m

       !
       !          Calculate net deposition and erosion
       cff=0.0_sp
       cff2=0.0_sp
       DO ised=1,NST
          cff=cff+sedbed%settling_flux(i,ised)
          cff2=cff2+sedbed%ero_flux(i,ised)
          dep_mass(ised)=0.0_sp
          IF ((sedbed%ero_flux(i,ised)-sedbed%settling_flux(i,ised)).lt.      &
               &             0.0_sp) THEN
             dep_mass(ised)=sedbed%settling_flux(i,ised)-              &
                  &                sedbed%ero_flux(i,ised)
           END IF
       END DO
       bottom(i,idnet)=cff-cff2

       IF ( cff-cff2.GT.0.0_sp) THEN ! NET deposition
          !  Deposition. Determine if we need to create a new bed layer 

          if( COHESIVE_BED .or. MIXED_BED)then
             !
             !  Calculate tau_crit of deposited bed
             !
             !  Calculate new mass in top layer
             bmz(1) = 0.0_sp
             DO ised=1,NST
                bmz(1) = bmz(1)+sedbed%bed_mass(i,1,nnew,ised)
             END DO
             IF (Nbed.GT.1) THEN
                !  Average of cff1 and cff2, where
                !    cff1 = linear extension of previous tcr slope to new surface
                !    cff2 = minimum deposition
                cff1 = bed(i,1,ibtcr) -                               &
                     &             bottom(i,idnet) *                                  &
                     &             (bed(i,2,ibtcr)-bed(i,1,ibtcr)) /                &
                     &             (bmz(1)-bottom(i,idnet))
                cff2 = MAX( tau_w(i) , tcr_min )
                bed(i,1,ibtcr) = MIN(bed(i,1,ibtcr),                &
                     &                                 MAX( 0.5_sp*(cff1+cff2), cff2 ))
                !                cff1 = bottom(i,idnet)/MAX(bmz(1),eps)
                !                cff2 = bed(i,1,ibtcr)+bed(i,2,ibtcr)-2*tcr_min
                !                bed(i,1,ibtcr) =bed(i,1,ibtcr) - cff1*cff2 
             ELSE
                !  TODO: Not sure what mud tau_crit should be for single-layer bed.
                !  Try weighted average of dep and old value
                cff1 = (bed(i,1,ibtcr)*(bmz(1)-bottom(i,idnet))     &
                     &                  + cff2*bottom(i,idnet)) / bmz(1)
                bed(i,1,ibtcr) = MAX( cff1, cff2 )
             END IF
          endif
          ! (no test for age here)
          bed(i,1,iaged)=T_model
          IF(bed(i,1,ithck).gt.                                   &
               &             MAX(bottom(i,iactv),newlayer_thick)) THEN
             ! Top layer is too thick
              IF (Nbed.gt.2) THEN
                IF(bed(i,2,ithck).lt.minlayer_thick) THEN
                   ! Layer 2 is smaller than minimum size
                   ! Instead of pushing down all layers, just combine top 2 layers
                   cff=0.0_sp
                   cff1=0.0_sp
                   cff2=0.0_sp
                   DO ised=1,NST
                      cff =cff +dep_mass(ised)
                      cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                      cff2=cff2+sedbed%bed_mass(i,2,nnew,ised)
                   END DO
                   if(COHESIVE_BED .or. MIXED_BED)then
                      !
                      !   Assign new tau_crit at 2nd layer
                      !
                      bed(i,2,ibtcr)= bed(i,2,ibtcr) -                &
                           &   (cff1-cff)/(cff1+cff2)*(bed(i,3,ibtcr)-bed(i,1,ibtcr))
                      !         bed(i,2,ibtcr) = 0.5_sp*( bed(i,1,ibtcr)+         &
                      !     &      bed(i,2,ibtcr) )
                   endif
                   !

                   !  Update bed mass
                   DO ised=1,NST
                      sedbed%bed_mass(i,2,nnew,ised)=                    &
                           &    MAX(sedbed%bed_mass(i,2,nnew,ised)+           &
                           &        sedbed%bed_mass(i,1,nnew,ised)-               &
                           &               dep_mass(ised),0.0_sp)
                      sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                   END DO
                   ! ALA - average time and porosity
                   ! ALA CHECK WITH CRS cff1 or cff1-cff for first layer
                   bed(i,2,iaged)=(bed(i,1,iaged)*cff1+         &
                        &      bed(i,2,iaged)*cff2)/(cff1+cff2)
                   bed(i,1,iaged)=T_model
                   bed(i,2,iporo)=(bed(i,1,iporo)*cff1+         &
                        &      bed(i,2,iporo)*cff2)/(cff1+cff2)
                   ! ALA CHECK WITH CRS POROSITY OF 1ST LAYER
                   bed(i,1,iporo)=bed(i,1,iporo)
                ELSE
                   ! Layer 2 is > minlayer thick, need another layer
                   !  Combine bottom layers.
                   cff1=0.0_sp
                   cff2=0.0_sp
                   DO ised=1,NST
                      cff1=cff1+sedbed%bed_mass(i,Nbed-1,nnew,ised)
                      cff2=cff2+sedbed%bed_mass(i,Nbed,nnew,ised)
                   END DO
                   bed(i,Nbed,iporo)=  (bed(i,Nbed-1,iporo)*cff1+                &
                        &        bed(i,Nbed,iporo)*cff2)/(cff1+cff2)
                   bed(i,Nbed,iaged)=  (bed(i,Nbed-1,iaged)*cff1+                &
                        &        bed(i,Nbed,iaged)*cff2)/(cff1+cff2)
                   if(COHESIVE_BED .or. MIXED_BED)then
                      !
                      !   Assign tcrit at top of new bottom bed tcrit for Nbed-1
                      bed(i,Nbed,ibtcr)= bed(i,Nbed-1,ibtcr)
                      !   bed(i,Nbed,ibtcr)=bed(i,Nbed-1,ibtcr) +                &
                      !     &         cff2*(bed(i,Nbed,ibtcr)-bed(i,Nbed-1,ibtcr))/cff1
                   endif

                   DO ised=1,NST
                      sedbed%bed_mass(i,Nbed,nnew,ised)=                 &
                           &   sedbed%bed_mass(i,Nbed-1,nnew,ised)+          &
                           &   sedbed%bed_mass(i,Nbed  ,nnew,ised)
                   END DO
                   !
                   !  Push layers down.
                   DO k=Nbed-1,2,-1
                      bed(i,k,iporo)=bed(i,k-1,iporo)
                      bed(i,k,iaged)=bed(i,k-1,iaged)
                      DO ised =1,NST
                         sedbed%bed_mass(i,k,nnew,ised)=                 &
                              &    sedbed%bed_mass(i,k-1,nnew,ised)
                      END DO
                      if(COHESIVE_BED .or. MIXED_BED)then
                         bed(i,k,ibtcr)=bed(i,k-1,ibtcr)
                      endif
                   END DO
                   !  Set new top parameters for top 2 layers
                   if(COHESIVE_BED .or. MIXED_BED)then
                      !  Tau_crit at top already determined
                      !  Set tau_crit at 2nd layer
                      cff=0.0_sp
                      cff1=0.0_sp
                      DO ised=1,NST
                         cff =cff +dep_mass(ised)
                         cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                      END DO

                      cff2 = (bed(i,2,ibtcr)-bed(i,1,ibtcr))/cff1
                      bed(i,2,ibtcr) = bed(i,1,ibtcr)+cff*cff2
                   endif

                   DO ised=1,NST
                      sedbed%bed_mass(i,2,nnew,ised)=                    &
                           &       MAX(sedbed%bed_mass(i,2,nnew,ised)-           &
                           &         dep_mass(ised),0.0_sp)
                      sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                   END DO
                END IF
             ELSEIF (Nbed.eq.2) THEN 
                ! NBED=2
                !
                if(COHESIVE_BED .or. MIXED_BED)then
                   !  Tau_crit at top already determined 
                   !  Set tau_crit at 2nd layer
                   cff=0.0_sp
                   cff1=0.0_sp
                   DO ised=1,NST
                      cff =cff +dep_mass(ised)
                      cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                   END DO
                   cff2 = (bed(i,2,ibtcr)-bed(i,1,ibtcr))/cff1
                   bed(i,2,ibtcr) = bed(i,1,ibtcr)+cff*cff2
                endif
                cff1=0.0_sp
                cff2=0.0_sp
                DO ised=1,NST
                   cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                   cff2=cff2+sedbed%bed_mass(i,2,nnew,ised)
                END DO
                DO ised=1,NST
                   sedbed%bed_mass(i,2,nnew,ised)=                       &
                        &     MAX(sedbed%bed_mass(i,2,nnew,ised)+              &
                        &         sedbed%bed_mass(i,1,nnew,ised)-                  &
                        &              dep_mass(ised),0.0_sp)
                   sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                END DO
                ! ALA - average time and porosity
                bed(i,2,iaged)=(bed(i,1,iaged)*cff1+            &
                     &                      bed(i,2,iaged)*cff2)/(cff1+cff2)
                bed(i,1,iaged)=T_model
                bed(i,2,iporo)=(bed(i,1,iporo)*cff1+            &
                     &                      bed(i,2,iporo)*cff2)/(cff1+cff2)
                ! ALA CHECK WITH CRS POROSITY OF 1ST LAYER
                bed(i,1,iporo)=bed(i,1,iporo)                    
             ELSE
                ! NBED=1
             END IF
          ELSE
             ! Net deposition has occurred, but no new bed layer was created 
          END IF
       ELSE
          ! Net erosion occurred
          if(COHESIVE_BED .or. MIXED_BED)then
             bmz(1) = 0.0_sp
             DO ised=1,NST
                bmz(1) = bmz(1)+sedbed%bed_mass(i,1,nnew,ised)
             END DO
             !  recalc tc for top of bed, based on linear
             !  interpolation and mass removed / orig. mass in top layer
             bed(i,1,ibtcr)=bed(i,1,ibtcr)+                        &
                  &             MIN(1.0_sp,-bottom(i,idnet)/MAX(eps,bmz(1)))*      &
                  &             MAX(0.0_sp,(bed(i,2,ibtcr)-bed(i,1,ibtcr)))     
          endif
          bed(i,1,iaged)=T_model
          IF (Nbed.eq.2) THEN 
             ! NBED=2
             DO ised=1,NST
                sedbed%bed_mass(i,2,nnew,ised)=                       &
                     &      MAX(sedbed%bed_mass(i,2,nnew,ised)+              &
                     &          sedbed%bed_mass(i,1,nnew,ised)-                  &
                     &           dep_mass(ised),0.0_sp)
                sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
             END DO
          ELSEIF (Nbed.eq.1) THEN
             ! ALF NO NEED TO DO ANYTHING
          ELSE
          END IF
       END IF

       ! Recalculate thickness and fractions for all layers.
       DO k=1,Nbed
          cff3=0.0_sp
          DO ised=1,NST
             cff3=cff3+sedbed%bed_mass(i,k,nnew,ised)
          END DO
          IF (cff3.eq.0.0_sp) THEN 
             cff3=eps 
          END IF
          bed(i,k,ithck)=0.0_sp
          DO ised=1,NST
             sedbed%bed_frac(i,k,ised)=sedbed%bed_mass(i,k,nnew,ised)/cff3
             bed(i,k,ithck)=MAX(bed(i,k,ithck)+                 &
                  &       sedbed%bed_mass(i,k,nnew,ised)/               &
                  &            (sed(ised)%Srho*                          &
                  &         (1.0_sp-bed(i,k,iporo))),0.0_sp)
          END DO
       END DO
    END DO J_LOOP2

    J_LOOP3 : DO i=1,m
       IF (bottom(i,iactv).gt.bed(i,1,ithck)) THEN

          IF (Nbed.eq.1) THEN
             bottom(i,iactv)=bed(i,1,ithck)
          ELSE
             thck_to_add=bottom(i,iactv)-bed(i,1,ithck)
             thck_avail=0.0_sp
             Ksed=1                                        ! initialize
             DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                   thck_avail=thck_avail+bed(i,k,ithck)
                   Ksed=k
                END IF
             END DO
             !
             !  Catch here if there was not enough bed material.
             !

             IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,iactv)=bed(i,1,ithck)+thck_avail
                thck_to_add=thck_avail
             END IF
             !
             !  Update bed mass of top layer and fractional layer.
             !
             cff2=MAX(thck_avail-thck_to_add,0.0_sp)/                  &
                  &             MAX(bed(i,Ksed,ithck),eps)	        
             DO ised=1,NST
                cff1=0.0_sp
                DO k=1,Ksed
                   cff1=cff1+sedbed%bed_mass(i,k,nnew,ised)
                END DO
                cff3=cff2*sedbed%bed_mass(i,Ksed,nnew,ised)
                sedbed%bed_mass(i,1   ,nnew,ised)=cff1-cff3
                sedbed%bed_mass(i,Ksed,nnew,ised)=cff3
             END DO

             !
             !  Update thickness of fractional layer ksource_sed.
             !
             bed(i,Ksed,ithck)=MAX(thck_avail-thck_to_add,0.0_sp)

             if(COHESIVE_BED .or. MIXED_BED)then
                !
                ! Update tau_cr of fractional layer
                if(Nbed>2)then
                    bed(i,Ksed,ibtcr) = bed(i,Ksed+1,ibtcr)-              &
                     &                 cff2*(bed(i,Ksed+1,ibtcr)-bed(i,Ksed,ibtcr))
                end if
             endif

             !
             !  Update bed fraction of top layer.
             !
             cff3=0.0_sp
             DO ised=1,NST
                cff3=cff3+sedbed%bed_mass(i,1,nnew,ised)
             END DO
             IF (cff3.eq.0.0_sp) THEN
                cff3=eps
             END IF
             DO ised=1,NST
                sedbed%bed_frac(i,1,ised)=sedbed%bed_mass(i,1,nnew,ised)/cff3
             END DO
             !
             !  Upate bed thickness of top layer.
             !
             bed(i,1,ithck)=bottom(i,iactv)
             !
             !  Pull all layers closer to the surface.
             !
             DO k=Ksed,Nbed
                ks=Ksed-2
                bed(i,k-ks,ithck)=bed(i,k,ithck)
                bed(i,k-ks,iporo)=bed(i,k,iporo)
                bed(i,k-ks,iaged)=bed(i,k,iaged)

                if(COHESIVE_BED .or. MIXED_BED)then
                   bed(i,k-ks,ibtcr)=bed(i,k,ibtcr)
                endif

                DO ised=1,NST
                   sedbed%bed_frac(i,k-ks,ised)=sedbed%bed_frac(i,k,ised)
                   sedbed%bed_mass(i,k-ks,nnew,ised)=sedbed%bed_mass(i,k,nnew,ised)
                END DO
             END DO
             !
             !  Add new layers onto the bottom. Split what was in the bottom layer to
             !  fill these new empty cells. ("ks" is the number of new layers).
             !
             ks=Ksed-2
             cff=1.0_sp/REAL(ks+1,sp)
             if (COHESIVE_BED .or. MIXED_BED)then
                !alf              cff1 = cff*(tcr_max-bed(i,Nbed-ks,ibtcr))
                !alf          Use the value from the bottom layer as max.
                cff1 = cff*(bed(i,Nbed,ibtcr)-bed(i,Nbed-ks,ibtcr))

                ! Interpolate bottom layer to tau_crit_max
                ! TODO: should be the reference profile and not tau_crit_max
                DO k=Nbed-1,Nbed-ks,-1

                   bed(i,k,ibtcr)=bed(i,Nbed-ks,ibtcr)+           &
                        &                       REAL(k-Nbed+ks,sp) * cff1
                ENDDO
             endif

             ! ALA CHECK WITH CRS about bed_frac 
             nnn=0
             DO ised=1,NST
                nlysm(ised)=newlayer_thick*REAL(ks+1,sp)*          &
                     &                  (sed(ised)%Srho*                                 &
                     &                  (1.0_sp-bed(i,Nbed-ks,iporo)))*               &
                     &                  sedbed%bed_frac(i,Nbed-ks,ised)
                IF (ks.gt.0) THEN
                   IF (sedbed%bed_mass(i,Nbed-ks,nnew,ised).gt.         &
                        &          nlysm(ised)) THEN
                      nnn=nnn+1
                      nlysm(ised)=  newlayer_thick*REAL(ks,sp)*        &
                           &      (sed(ised)%Srho*                         &
                           &      (1.0_sp-bed(i,Nbed-ks,iporo)))*       &
                           &      sedbed%bed_frac(i,Nbed-ks,ised)
                   END IF
                END IF
             END DO
             IF (nnn.eq.NST) THEN
                bed(i,Nbed,ithck)=bed(i,Nbed-ks,ithck)-           &
                     &         newlayer_thick*REAL(ks,sp)
                DO ised=1,NST
                   sedbed%bed_mass(i,Nbed,nnew,ised)=                      &
                        &       sedbed%bed_mass(i,Nbed-ks,nnew,ised)-nlysm(ised)
                END DO
                DO k=Nbed-1,Nbed-ks,-1
                   bed(i,k,ithck)=newlayer_thick
                   bed(i,k,iaged)=bed(i,Nbed-ks,iaged)
                   DO ised=1,NST
                      sedbed%bed_frac(i,k,ised)=sedbed%bed_frac(i,Nbed-ks,ised)
                      sedbed%bed_mass(i,k,nnew,ised)=                      &
                           &    nlysm(ised)/REAL(ks,sp)
                   END DO
                END DO
             ELSE
                cff=1.0_sp/REAL(ks+1,sp)
                DO k=Nbed,Nbed-ks,-1
                   bed(i,k,ithck)=bed(i,Nbed-ks,ithck)*cff
                   bed(i,k,iaged)=bed(i,Nbed-ks,iaged)
                   DO ised=1,NST
                      sedbed%bed_frac(i,k,ised)=sedbed%bed_frac(i,Nbed-ks,ised)
                      sedbed%bed_mass(i,k,nnew,ised)=                      &
                           &     sedbed%bed_mass(i,Nbed-ks,nnew,ised)*cff
                   END DO
                END DO
             END IF
          END IF  ! Nbed > 1
       END IF  ! increase top bed layer
       if(MIXED_BED)then
          !
          !         Update cohesive property of seds in top layer
          cff1 = 0.0_sp
          DO ised=1,NCS
             cff1=cff1+sedbed%bed_frac(i,1,ised)
          END DO
          bottom(i,idprp)=min(max((cff1-transN)/                  &
               &     (transC-transN),0.0_sp),1.0_sp)
       endif

       if(SED_FLOCS .and. SED_DEFLOC)then
          alpha = MIN(DTsed/t_dfloc,1.0_sp)
          DO k=2,Nbed
             cff3=0.0_sp
             DO ised=1,NST
                cff3=cff3+sedbed%bed_mass(i,k,nnew,ised)
             END DO
             cff1 = 0.0_sp
             DO ised=1,NCS
                cff1 = cff1+sedbed%bed_mass(i,k,nnew,ised)
             END DO
             DO ised=1,NCS
                masseq(ised)=mud_frac_eq(ised)*cff1
                sedbed%bed_mass(i,k,nnew,ised)=alpha*masseq(ised)+           &
                     &                  sedbed%bed_mass(i,k,nnew,ised)*(1.0_sp-alpha)
                sedbed%bed_frac(i,k,ised)=sedbed%bed_mass(i,k,nnew,ised)/cff3
             END DO
          END DO
       endif

       if(COHESIVE_BED .or. MIXED_BED)then
          !
          !  Key to this algorithm: "mass available depth" (mad) [kg/m2]
          !  present tau crit profile (tcp)
          !  representative tau crit profile (tcr)          
          !  Compute depth and mass depth of sediment column
          !  Compute cumulative depth
          !  Compute mass depth
          bmz(1) = 0.0_sp
          DO ised=1,NST
             bmz(1) = bmz(1)+sedbed%bed_mass(i,1,nnew,ised)
          ENDDO
          DO k=2,Nbed
             bmz(k) = bmz(k-1)
             DO ised=1,NST
                bmz(k)=bmz(k)+sedbed%bed_mass(i,k,nnew,ised)
             ENDDO
          ENDDO

          !  Calculate representative critical shear stress profile
          !  Note that the values are for the TOP of the layer...we
          !  assume the bottom of the bottom layer has tcr = tcr_max
          tcr(1) = tcr_min
          DO k=2,Nbed
             tcr(k) = tcr_min
             IF (bmz(k-1).GT.eps) THEN
                !                tcr(k) = exp((log(bmz(k-1))-                            &
                !     &               bottom(i,idoff))/                                &
                !     &               bottom(i,idslp))
                tcr(k) = exp((log(bmz(k-1))-                            &
                     &                    tcr_off)/                                 &
                     &                    tcr_slp)
             ENDIF
             tcr(k) = MIN( MAX( tcr(k), tcr_min), tcr_max )
          ENDDO
          !          ks=Ksed-2
          !          DO k=Nbed,Nbed-ks,-1
          !             bed(i,k,ibtcr)=tcr(k)
          !          ENDDO
          !  Relax tau crit profile bottom(i,k,ibtcr) 
          !  towards representative profile tcr...100 x slower for "swelling"
          !  than consolidation
          !  TODO - make the factor an input parameter
          IF( bed(i,1,ibtcr).LE.tcr(1)) THEN
             !             alpha = MIN(DTsed/bottom(i,idtim),1.0_sp)
             alpha = MIN(DTsed/tcr_tim,1.0_sp)
          ELSE
             !             alpha = MIN(DTsed/(100.0_sp*bottom(i,idtim)),1.0_sp)
             alpha = MIN(DTsed/(100.0_sp*tcr_tim),1.0_sp)
          ENDIF
          DO k=1,Nbed
             bed(i,k,ibtcr)=bed(i,k,ibtcr)+                         &
                  &            alpha*(tcr(k)-bed(i,k,ibtcr))
             bed(i,k,ibtcr)=                                          &
                  &          MIN( MAX( bed(i,k,ibtcr), tcr_min), tcr_max )
          ENDDO
          !ALA    Enforce last layer = tcr_max
          !ALA          bed(i,Nbed,ibtcr)= tcr_max
       endif !/* cohesive bed */

    END DO J_LOOP3

    !
    !-----------------------------------------------------------------------
    ! Store old bed thickness.
    !-----------------------------------------------------------------------
    !
    if (SED_MORPH)then
       do i=1,m
          sedbed%bed_thick(i,nnew)=0.0_sp
          DO k=1,Nbed
             sedbed%bed_thick(i,nnew)=sedbed%bed_thick(i,nnew)+                  &
                  &                       bed(i,k,ithck)
          END DO
       end do

# if defined (MULTIPROCESSOR)
       if(par)then
          call aexchange(nc,myid,nprocs,sedbed%bed_thick)
       endif
# endif
       bottom(i,morph) = sedbed%bed_thick(i,nnew)
    endif

# if defined (MULTIPROCESSOR)
    if(par)then
       allocate(temp(0:mt,1:nbed))
       DO ised=1,NST
          temp = 0.0_sp
          temp = sedbed%bed_frac(:,:,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_frac(:,:,ised) = temp
          temp = 0.0_sp
          temp = sedbed%bed_mass(:,:,nnew,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_mass(:,:,nnew,ised) = temp
       end do

       do ibed=1,MBEDP
          temp = 0.0_sp      
          temp = bed(:,:,ibed)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          bed(:,:,ibed)=temp
       end do
       deallocate(temp)
    endif
# endif

    if(dbg_set(dbg_sbr)) write(ipt,*) "End:  calc_sed_bed_cohesive"

    return



  end Subroutine calc_sed_bed_cohesive


!=======================================================================
!                                                                      !
!  This routine computes sediment bed layer stratigraphy.              !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
  Subroutine calc_sed_bed_cohesive2

    implicit none

    !
    !  Local variable declarations.
    !
    integer :: Kbed, i, ised,ibed, j, k, ks
    integer :: bnew, nnn, cnt,cc,jj

    real(sp), parameter :: eps = 1.0E-14_sp

    real(sp) :: cff, cff1, cff2, cff3,cff4
    real(sp) :: thck_avail, thck_to_add

    real(sp), dimension(NST) :: nlysm

    real(sp), dimension(NST) :: dep_mass
    real(sp), dimension(0:mt) :: tau_w

    real(sp) :: pcoh 

    real(sp) :: alpha, bzactv, frt, tcb_top, tcb_bot
    real(sp), dimension(Nbed) :: tcr

    real(sp), dimension(Nbed) :: bmz

    real(sp), dimension(NCS) :: masseq

    real(sp),allocatable,dimension(:,:)::temp

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start:  calc_sed_bed_cohesive"

    if(BEDLOAD)then
       bnew=nnew
    else
       bnew=nstp
    endif

    ! KLUDGE alert
    ! minlayer_thick = 0.0005
    minlayer_thick = newlayer_thick
    !
    !-----------------------------------------------------------------------
    ! Compute sediment bed layer stratigraphy.
    !-----------------------------------------------------------------------
    !
    if(BEDLOAD_MPM.or.SUSLOAD)then
       tau_w = taub
    endif

    !
    !-----------------------------------------------------------------------
    !  Update bed properties according to ero_flux and dep_flux.
    !-----------------------------------------------------------------------
    !
    if(SUSLOAD)then
       NODE_LOOP : DO i=1,m
          !
          !  The deposition and resuspension of sediment on the bottom "bed"
          !  is due to precipitation flux FC(:,0), already computed, and the
          !  resuspension (erosion, hence called ero_flux). The resuspension is
          !  applied to the bottom-most grid box value qc(:,1) so the total mass
          !  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
          !  after both fluxes are applied.
          !
          SED_LOOP: DO ised=1,NST  
             dep_mass(ised)=0.0_sp
             if (SED_MORPH)then
                ! Apply morphology factor.
                sedbed%ero_flux(i,ised)=sedbed%ero_flux(i,ised)*sed(ised)%morph_fac
                sedbed%settling_flux(i,ised)=sedbed%settling_flux(i,ised)*            &
                     &     sed(ised)%morph_fac
             endif

             !  Update bed mass arrays.
             sedbed%bed_mass(i,1,nnew,ised)=MAX(sedbed%bed_mass(i,1,bnew,ised)-     &
                  &     (sedbed%ero_flux(i,ised)-                                       &
                  &     sedbed%settling_flux(i,ised)),                                  &
                  &     0.0_sp)
             DO k=2,Nbed
                sedbed%bed_mass(i,k,nnew,ised)=sedbed%bed_mass(i,k,nstp,ised)
             END DO
          END DO SED_LOOP

          cff3=0.0_sp
          DO ised=1,NST
             cff3=cff3+sedbed%bed_mass(i,1,nnew,ised)
          END DO

          IF (cff3.eq.0.0_sp) THEN 
             cff3=eps 
          END IF
          bed(i,1,ithck)=0.0_sp
          DO ised=1,NST
             sedbed%bed_frac(i,1,ised)=sedbed%bed_mass(i,1,nnew,ised)/cff3
             bed(i,1,ithck)=MAX(bed(i,1,ithck)+                      &
                  &     sedbed%bed_mass(i,1,nnew,ised)/                                &
                  &     (sed(ised)%Srho*                                         &
                  &     (1.0_sp-bed(i,1,iporo))),0.0_sp)
          END DO

       END DO NODE_LOOP
    endif  !/* SUSPLOAD section */

    !
    !-----------------------------------------------------------------------
    !  At this point, all deposition or erosion is complete, and
    !  has been added/subtracted to top layer. Thickness has NOT been corrected.
    !-----------------------------------------------------------------------
    !
    NODE_LOOP_CB : DO i=1,m
       !    Calculate active layer thickness, bottom(i,iactv).
       !    (trunk version allows this to be zero...this has minimum of 6*D50)
       bottom(i,iactv)=MAX(0.0_sp,                                 &
            &                    0.007_sp*                               &
            &                    (tau_w(i)-bottom(i,itauc)))+   &
            &                          6.0_sp*bottom(i,isd50)

       if (COHESIVE_BED)then 
          bottom(i,iactv)=MAX(0.0_sp,                                 &
               &                       0.007_sp*                               &
               &                       (tau_w(i)-bed(i,1,ibtcr)))+    &
               &                       6.0_sp*bottom(i,isd50)
       endif

       if (MIXED_BED)then
          cff1=             MAX(0.0_sp,                                 &
               &                          0.007_sp*                               &
               &                          (tau_w(i)-bed(i,1,ibtcr)))+    &
               &                          6.0_sp*bottom(i,isd50)
          cff2=             MAX(0.0_sp,                                 &
               &                          0.007_sp*                               &
               &                          (tau_w(i)-bottom(i,itauc)))+   &
               &                          6.0_sp*bottom(i,isd50)
          bottom(i,iactv)=MAX(cff1,cff2)
       endif
       if (SED_MORPH)then
          ! Apply morphology factor.
          bottom(i,iactv)=MAX(bottom(i,iactv)*sed(1)%morph_fac,      &
               &          bottom(i,iactv))
       endif
       if (COHESIVE_BED .or. MIXED_BED)then
          !  Find first layer with tc > tb
          !  Remember, the critical stresses apply to the TOP of each layer
          Kbed = 1
          bzactv = 0.0_sp
          frt = 0.0_sp ! CRS
          cff1 = tau_w(i)
          tcb_top = bed(i,1,ibtcr)
          tcb_bot = bed(i,2,ibtcr)
          if (MIXED_BED)then
             !         Calc. tau crit for bottom of next layer
             !         Update cohesive property of seds in top layer
             cff3 = 0.0_sp
             cff4 = 1.0_sp
             DO ised=1,NCS
                cff3=cff3+sedbed%bed_frac(i,1,ised)
                cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,1,ised)
             END DO
             DO ised=NCS+1,NST
                cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,1,ised)
             ENDDO
             pcoh=min(max((cff3-transN)/(transC-transN)        &
                  &      ,0.0_sp),1.0_sp)
             tcb_top = pcoh*bed(i,1,ibtcr)+(pcoh-1.0_sp)*cff4   !BF
          endif

          IF(cff1 .GT. tcb_top)THEN
             !         Calculate tcb_temp for next layer
             tcb_bot = bed(i,Kbed+1,ibtcr)
             if (MIXED_BED)then
                !            Recalculate cohesive fraction and mean grain tau crit
                !            Note that Kbed is used for grain props, and Kbed+1 for
                !            bed tau crit at bottom of layer Kbed
                cff3 = 0.0_sp
                cff4 = 1.0_sp
                DO ised=1,NCS
                   cff3=cff3+sedbed%bed_frac(i,Kbed,ised)
                   cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,Kbed,ised)
                END DO
                DO ised=NCS+1,NST
                   cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,Kbed,ised)
                ENDDO
                !            Calculate cohesive behavior and blended tau crit
                pcoh=min(max((cff3-transN)/(transC-transN),    &
                     &            0.0_sp),1.0_sp)
                tcb_bot =pcoh*bed(i,Kbed+1,ibtcr)+(pcoh-1.0_sp)*cff4  !BF
             endif
             DO WHILE ( (Kbed.LE.(Nbed-1)) .AND.                        &
                  &       (cff1 .GT. tcb_bot))
                !ALA   Instead of adding entire 2nd layer, add just what is needed
                !ALA!         Add entire layer
                !ALA          bzactv = bzactv + bed(i,Kbed,ithck)
                !       Add thickness equivalent to difference between cff1 and tcb_bot
                IF (cff1.GT.bed(i,Kbed+1,ibtcr)) THEN
                   bzactv = bzactv + bed(i,Kbed,ithck)
                   tcb_top = tcb_bot
                   tcb_bot = bed(i,Kbed+1,ibtcr)
                ELSE
                   bzactv = bzactv + MIN(bed(i,Kbed,ithck),           &
                        &                  (MAX(0.0_sp,0.007_sp*(cff1-tcb_bot))))
                   tcb_top = cff1
                   tcb_bot = bed(i,Kbed+1,ibtcr)   
                ENDIF
                if (MIXED_BED)then
                   !               Recalculate cohesive fraction and mean grain tau crit
                   cff3 = 0.0_sp
                   cff4 = 1.0_sp
                   DO ised=1,NCS
                      cff3=cff3+sedbed%bed_frac(i,Kbed,ised)
                      cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,Kbed,ised)
                   END DO
                   DO ised=NCS+1,NST
                      cff4=cff4*sed(ised)%tau_ce**sedbed%bed_frac(i,Kbed,ised)
                   ENDDO
                   !               Calculate cohesive behavior and blended tau crit
                   pcoh=min(max((cff3-transN)/(transC-transN), &
                        &               0.0_sp),1.0_sp)
                   tcb_bot = pcoh*bed(i,Kbed+1,ibtcr)+(pcoh-1.0_sp)*cff4
                endif
                Kbed = Kbed+1
             ENDDO
             frt = MAX(0.0_sp,(cff1-tcb_top) ) /                        &
                  &          MAX( eps,                                             &
                  &       ( tcb_bot-tcb_top ))
          ENDIF
          bzactv = bzactv+frt*bed(i,Kbed,ithck)
          bzactv = MAX( bzactv, 6.0_sp*bottom(i,isd50) ) !CRS
          bottom(i,iactv)=min(bottom(i,iactv),bzactv) !?CRS
          !      bottom(i,iactv)=max(bottom(i,iactv),bzactv) !?CRS
       endif !/* defined COHESIVE_BED || defined MIXED_BED */

    END DO NODE_LOOP_CB

    NODE_LOOP2 : DO i=1,m

       !
       !          Calculate net deposition and erosion
       cff=0.0_sp
       cff2=0.0_sp
       DO ised=1,NST
          cff=cff+sedbed%settling_flux(i,ised)
          cff2=cff2+sedbed%ero_flux(i,ised)
          dep_mass(ised)=0.0_sp
          IF ((sedbed%ero_flux(i,ised)-sedbed%settling_flux(i,ised)).lt.      &
               &       0.0_sp) THEN
             dep_mass(ised)=sedbed%settling_flux(i,ised)-              &
                  &       sedbed%ero_flux(i,ised)
          END IF
       END DO
       bottom(i,idnet)=cff-cff2

       IF ( cff-cff2.GT.0.0_sp) THEN ! NET deposition
          !  Deposition. Determine if we need to create a new bed layer 
          if (COHESIVE_BED .or. MIXED_BED)then
             !
             !  Calculate tau_crit of deposited bed
             !
             !  Calculate new mass in top layer
             bmz(1) = 0.0_sp
             DO ised=1,NST
                bmz(1) = bmz(1)+sedbed%bed_mass(i,1,nnew,ised)
             END DO
             IF (Nbed.GT.1) THEN
                !  Average of cff1 and cff2, where
                !    cff1 = linear extension of previous tcr slope to new surface
                !    cff2 = minimum deposition
                cff1 = bed(i,1,ibtcr) -                               &
                     &             bottom(i,idnet) *                                  &
                     &             (bed(i,2,ibtcr)-bed(i,1,ibtcr)) /                &
                     &             (bmz(1)-bottom(i,idnet))
                cff2 = MAX( tau_w(i) , tcr_min )
                bed(i,1,ibtcr) = MIN(bed(i,1,ibtcr),                &
                     &                         MAX( 0.5_sp*(cff1+cff2), cff2 ))
                !            cff1 = bottom(i,idnet)/MAX(bmz(1),eps)
                !            cff2 = bed(i,1,ibtcr)+bed(i,2,ibtcr)-2*tcr_min
                !            bed(i,1,ibtcr) =bed(i,1,ibtcr) - cff1*cff2 
             ELSE
                !  TODO: Not sure what mud tau_crit should be for single-layer bed.
                !  Try weighted average of dep and old value
                cff1 = (bed(i,1,ibtcr)*(bmz(1)-bottom(i,idnet))     &
                     &              + cff2*bottom(i,idnet)) / bmz(1)
                cff2 = MAX( tau_w(i) , tcr_min )
                bed(i,1,ibtcr) = MAX( cff1, cff2 )
             END IF
             bed(i,1,ibtcr) = 0.1  !neeed be removed!!
          endif

          ! (no test for age here)
          bed(i,1,iaged)=T_model
          IF(bed(i,1,ithck).gt.                                   &
               &       MAX(bottom(i,iactv),newlayer_thick)) THEN
             ! Top layer is too thick
             IF (Nbed.gt.2) THEN
                IF(bed(i,2,ithck).lt.minlayer_thick) THEN
                   ! Layer 2 is smaller than minimum size
                   ! Instead of pushing down all layers, just combine top 2 layers
                   cff=0.0_sp
                   cff1=0.0_sp
                   cff2=0.0_sp
                   DO ised=1,NST
                      cff =cff +dep_mass(ised)
                      cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                      cff2=cff2+sedbed%bed_mass(i,2,nnew,ised)
                   END DO
                   if (COHESIVE_BED .or. MIXED_BED)then
                      !
                      !   Assign new tau_crit at 2nd layer
                      !
                      bed(i,2,ibtcr)= bed(i,2,ibtcr) -                &
                           &       (cff1-cff)/(cff1+cff2)*(bed(i,3,ibtcr)-bed(i,1,ibtcr))
                      !                  bed(i,2,ibtcr) = 0.5_sp*( bed(i,1,ibtcr)+         &
                      !     &                                        bed(i,2,ibtcr) )
                   endif
                   !

                   !  Update bed mass
                   DO ised=1,NST
                      sedbed%bed_mass(i,2,nnew,ised)=                    &
                           &       MAX(sedbed%bed_mass(i,2,nnew,ised)+           &
                           &           sedbed%bed_mass(i,1,nnew,ised)-               &
                           &           dep_mass(ised),0.0_sp)
                      sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                   END DO
                   ! ALA - average time and porosity
                   ! ALA CHECK WITH CRS cff1 or cff1-cff for first layer
                   bed(i,2,iaged)=(bed(i,1,iaged)*cff1+         &
                        &                   bed(i,2,iaged)*cff2)/(cff1+cff2)
                   bed(i,1,iaged)=T_model
                   bed(i,2,iporo)=(bed(i,1,iporo)*cff1+         &
                        &                   bed(i,2,iporo)*cff2)/(cff1+cff2)
                   ! ALA CHECK WITH CRS POROSITY OF 1ST LAYER
                   bed(i,1,iporo)=bed(i,1,iporo)
                ELSE
                   ! Layer 2 is > minlayer thick, need another layer
                   !  Combine bottom layers.
                   cff1=0.0_sp
                   cff2=0.0_sp
                   DO ised=1,NST
                      cff1=cff1+sedbed%bed_mass(i,Nbed-1,nnew,ised)
                      cff2=cff2+sedbed%bed_mass(i,Nbed,nnew,ised)
                   END DO
                   bed(i,Nbed,iporo)=                             &
                        &              (bed(i,Nbed-1,iporo)*cff1+                &
                        &               bed(i,Nbed,iporo)*cff2)/(cff1+cff2)
                   bed(i,Nbed,iaged)=                             &
                        &              (bed(i,Nbed-1,iaged)*cff1+                &
                        &               bed(i,Nbed,iaged)*cff2)/(cff1+cff2)
                   if (COHESIVE_BED.or.MIXED_BED)then
                      !
                      !   Assign tcrit at top of new bottom bed tcrit for Nbed-1
                      bed(i,Nbed,ibtcr)= bed(i,Nbed-1,ibtcr)
                      !              bed(i,Nbed,ibtcr)=bed(i,Nbed-1,ibtcr) +                &
                      !     &            cff2*(bed(i,Nbed,ibtcr)-bed(i,Nbed-1,ibtcr))/cff1
                   endif

                   DO ised=1,NST
                      sedbed%bed_mass(i,Nbed,nnew,ised)=                 &
                           &            sedbed%bed_mass(i,Nbed-1,nnew,ised)+          &
                           &           sedbed%bed_mass(i,Nbed  ,nnew,ised)
                   END DO
                   !
                   !  Push layers down.
                   DO k=Nbed-1,2,-1
                      bed(i,k,iporo)=bed(i,k-1,iporo)
                      bed(i,k,iaged)=bed(i,k-1,iaged)
                      DO ised =1,NST
                         sedbed%bed_mass(i,k,nnew,ised)=                 &
                              &               sedbed%bed_mass(i,k-1,nnew,ised)
                      END DO
                      if (COHESIVE_BED.or.MIXED_BED)then
                         bed(i,k,ibtcr)=bed(i,k-1,ibtcr)
                      endif
                   END DO
                   !  Set new top parameters for top 2 layers
                   if (COHESIVE_BED.or.MIXED_BED)then
                      !  Tau_crit at top already determined
                      !  Set tau_crit at 2nd layer
                      cff=0.0_sp
                      cff1=0.0_sp
                      DO ised=1,NST
                         cff =cff +dep_mass(ised)
                         cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                      END DO

                      cff2 = (bed(i,2,ibtcr)-bed(i,1,ibtcr))/cff1
                      bed(i,2,ibtcr) = bed(i,1,ibtcr)+cff*cff2
                   endif

                   DO ised=1,NST
                      sedbed%bed_mass(i,2,nnew,ised)=                    &
                           &               MAX(sedbed%bed_mass(i,2,nnew,ised)-           &
                           &                   dep_mass(ised),0.0_sp)
                      sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                   END DO
                END IF
             ELSEIF (Nbed.eq.2) THEN 
                !
                ! NBED=2
                !
                if (COHESIVE_BED .and. MIXED_BED)then
                   !  Tau_crit at top already determined 
                   !  Set tau_crit at 2nd layer
                   cff=0.0_sp
                   cff1=0.0_sp
                   DO ised=1,NST
                      cff =cff +dep_mass(ised)
                      cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                   END DO
                   cff2 = (bed(i,2,ibtcr)-bed(i,1,ibtcr))/cff1
                   bed(i,2,ibtcr) = bed(i,1,ibtcr)+cff*cff2
                endif
                cff1=0.0_sp
                cff2=0.0_sp
                DO ised=1,NST
                   cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                   cff2=cff2+sedbed%bed_mass(i,2,nnew,ised)
                END DO
                DO ised=1,NST
                   sedbed%bed_mass(i,2,nnew,ised)=                       &
                        &     MAX(sedbed%bed_mass(i,2,nnew,ised)+              &
                        &         sedbed%bed_mass(i,1,nnew,ised)-                  &
                        &              dep_mass(ised),0.0_sp)
                   sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                END DO
                ! ALA - average time and porosity
                bed(i,2,iaged)=(bed(i,1,iaged)*cff1+            &
                     &           bed(i,2,iaged)*cff2)/(cff1+cff2)
                bed(i,1,iaged)=T_model
                bed(i,2,iporo)=(bed(i,1,iporo)*cff1+            &
                     &           bed(i,2,iporo)*cff2)/(cff1+cff2)
                ! ALA CHECK WITH CRS POROSITY OF 1ST LAYER
                bed(i,1,iporo)=bed(i,1,iporo)                    
             ELSE
                ! NBED=1
             END IF
          ELSE
             ! Net deposition has occurred, but no new bed layer was created 
          END IF


       ELSE  ! Net erosion occurred
          if (COHESIVE_BED .or. MIXED_BED)then
            IF (Nbed.ge.2) THEN 
               bmz(1) = 0.0_sp
               DO ised=1,NST
                  bmz(1) = bmz(1)+sedbed%bed_mass(i,1,nnew,ised)
               END DO
               !  recalc tc for top of bed, based on linear
               !  interpolation and mass removed / orig. mass in top layer
               bed(i,1,ibtcr)=bed(i,1,ibtcr)+                        &
                    &    MIN(1.0_sp,-bottom(i,idnet)/MAX(eps,bmz(1)))*      &
                    &    MAX(0.0_sp,(bed(i,2,ibtcr)-bed(i,1,ibtcr)))     
            endif
          endif
          bed(i,1,iaged)=T_model
          IF (Nbed.eq.2) THEN 
             ! NBED=2
             DO ised=1,NST
                sedbed%bed_mass(i,2,nnew,ised)=                       &
                     &           MAX(sedbed%bed_mass(i,2,nnew,ised)+              &
                     &           sedbed%bed_mass(i,1,nnew,ised)-                  &
                     &           dep_mass(ised),0.0_sp)
                sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
             END DO
          ELSEIF (Nbed.eq.1) THEN
             ! ALF NO NEED TO DO ANYTHING
          ELSE
          END IF

       END IF

       ! Recalculate thickness and fractions for all layers.
       DO k=1,Nbed
          cff3=0.0_sp
          DO ised=1,NST
             cff3=cff3+sedbed%bed_mass(i,k,nnew,ised)
          END DO
          IF (cff3.eq.0.0_sp) THEN 
             cff3=eps 
          END IF
          bed(i,k,ithck)=0.0_sp
          DO ised=1,NST
             sedbed%bed_frac(i,k,ised)=sedbed%bed_mass(i,k,nnew,ised)/cff3
             bed(i,k,ithck)=MAX(bed(i,k,ithck)+                 &
                  &                  sedbed%bed_mass(i,k,nnew,ised)/               &
                  &                     (sed(ised)%Srho*                          &
                  &           (1.0_sp-bed(i,k,iporo))),0.0_sp)
          END DO
       END DO

    END DO NODE_LOOP2

    NODE_LOOP3 : DO i=1,m

       IF (bottom(i,iactv).gt.bed(i,1,ithck)) THEN
          IF (Nbed.eq.1) THEN
             bottom(i,iactv)=bed(i,1,ithck)
          ELSE
             thck_to_add=bottom(i,iactv)-bed(i,1,ithck)
             thck_avail=0.0_sp
             Kbed=1                                        ! initialize
             DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                   thck_avail=thck_avail+bed(i,k,ithck)
                   Kbed=k
                END IF
             END DO
             !
             !  Catch here if there was not enough bed material.
             !
             IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,iactv)=bed(i,1,ithck)+thck_avail
                thck_to_add=thck_avail
             END IF
             !
             !  Update bed mass of top layer and fractional layer.
             !
             cff2=MAX(thck_avail-thck_to_add,0.0_sp)/                  &
                  &       MAX(bed(i,Kbed,ithck),eps)	        
             DO ised=1,NST
                cff1=0.0_sp
                DO k=1,Kbed
                   cff1=cff1+sedbed%bed_mass(i,k,nnew,ised)
                END DO
                cff3=cff2*sedbed%bed_mass(i,Kbed,nnew,ised)
                sedbed%bed_mass(i,1   ,nnew,ised)=cff1-cff3
                sedbed%bed_mass(i,Kbed,nnew,ised)=cff3
             END DO
             !
             !  Update thickness of fractional layer ksource_sed.
             !
             bed(i,Kbed,ithck)=MAX(thck_avail-thck_to_add,0.0_sp)
             if (COHESIVE_BED .or. MIXED_BED)then
                !
                ! Update tau_cr of fractional layer
                if(Kbed<Nbed)bed(i,Kbed,ibtcr) = bed(i,Kbed+1,ibtcr)-              &
                     &        cff2*(bed(i,Kbed+1,ibtcr)-bed(i,Kbed,ibtcr))
             endif
             !
             !  Update bed fraction of top layer.
             !
             cff3=0.0_sp
             DO ised=1,NST
                cff3=cff3+sedbed%bed_mass(i,1,nnew,ised)
             END DO
             IF (cff3.eq.0.0_sp) THEN
                cff3=eps
             END IF
             DO ised=1,NST
                sedbed%bed_frac(i,1,ised)=sedbed%bed_mass(i,1,nnew,ised)/cff3
             END DO
             !
             !  Upate bed thickness of top layer.
             !
             bed(i,1,ithck)=bottom(i,iactv)
             !
             !  Pull all layers closer to the surface.
             !
             DO k=Kbed,Nbed
                ks=Kbed-2
                bed(i,k-ks,ithck)=bed(i,k,ithck)
                bed(i,k-ks,iporo)=bed(i,k,iporo)
                bed(i,k-ks,iaged)=bed(i,k,iaged)

                if(COHESIVE_BED.or.MIXED_BED)then
                   bed(i,k-ks,ibtcr)=bed(i,k,ibtcr)
                endif

                DO ised=1,NST
                   sedbed%bed_frac(i,k-ks,ised)=sedbed%bed_frac(i,k,ised)
                   sedbed%bed_mass(i,k-ks,nnew,ised)=sedbed%bed_mass(i,k,nnew,ised)
                END DO
             END DO

             !
             !  Add new layers onto the bottom. Split what was in the bottom layer to
             !  fill these new empty cells. ("ks" is the number of new layers).
             !
             ks=Kbed-2
             cff=1.0_sp/REAL(ks+1,sp)
             if (COHESIVE_BED.or.MIXED_BED)then

                !alf       cff1 = cff*(tcr_max)-bed(i,Nbed-ks,ibtcr))
                !alf       Use the value from the bottom layer as max.
                cff1 = cff*(bed(i,Nbed,ibtcr)-bed(i,Nbed-ks,ibtcr))


                ! Interpolate bottom layer to tau_crit_max
                ! TODO: should be the reference profile and not tau_crit_max
                DO k=Nbed-1,Nbed-ks,-1

                   bed(i,k,ibtcr)=bed(i,Nbed-ks,ibtcr)+           &
                        &              REAL(k-Nbed+ks,sp) * cff1
                ENDDO

             endif
             ! ALA CHECK WITH CRS about bed_frac 
             nnn=0
             DO ised=1,NST
                nlysm(ised)=newlayer_thick*REAL(ks+1,sp)*          &
                     &            (sed(ised)%Srho*                                 &
                     &            (1.0_sp-bed(i,Nbed-ks,iporo)))*               &
                     &             sedbed%bed_frac(i,Nbed-ks,ised)
                IF (ks.gt.0) THEN
                   IF (sedbed%bed_mass(i,Nbed-ks,nnew,ised).gt.         &
                        &                nlysm(ised)) THEN
                      nnn=nnn+1
                      nlysm(ised)=                                 &
                           &                newlayer_thick*REAL(ks,sp)*        &
                           &                   (sed(ised)%Srho*                         &
                           &                (1.0_sp-bed(i,Nbed-ks,iporo)))*       &
                           &                sedbed%bed_frac(i,Nbed-ks,ised)
                   END IF
                END IF
             END DO
             IF (nnn.eq.NST) THEN
                bed(i,Nbed,ithck)=bed(i,Nbed-ks,ithck)-           &
                     &                   newlayer_thick*REAL(ks,sp)
                DO ised=1,NST
                   sedbed%bed_mass(i,Nbed,nnew,ised)=                      &
                        &            sedbed%bed_mass(i,Nbed-ks,nnew,ised)-nlysm(ised)
                END DO
                DO k=Nbed-1,Nbed-ks,-1
                   bed(i,k,ithck)=newlayer_thick
                   bed(i,k,iaged)=bed(i,Nbed-ks,iaged)
                   DO ised=1,NST
                      sedbed%bed_frac(i,k,ised)=sedbed%bed_frac(i,Nbed-ks,ised)
                      sedbed%bed_mass(i,k,nnew,ised)=                      &
                           &                nlysm(ised)/REAL(ks,sp)
                   END DO
                END DO
             ELSE
                cff=1.0_sp/REAL(ks+1,sp)
                DO k=Nbed,Nbed-ks,-1
                   bed(i,k,ithck)=bed(i,Nbed-ks,ithck)*cff
                   bed(i,k,iaged)=bed(i,Nbed-ks,iaged)
                   DO ised=1,NST
                      sedbed%bed_frac(i,k,ised)=sedbed%bed_frac(i,Nbed-ks,ised)
                      sedbed%bed_mass(i,k,nnew,ised)=                      &
                           &                    sedbed%bed_mass(i,Nbed-ks,nnew,ised)*cff
                   END DO
                END DO
             END IF

          END IF  ! Nbed > 1

       END IF  ! increase top bed layer

       if (MIXED_BED)then
          !
          !     Update cohesive property of seds in top layer
          cff1 = 0.0_sp
          DO ised=1,NCS
             cff1=cff1+sedbed%bed_frac(i,1,ised)
          END DO
          bottom(i,idprp)=min(max((cff1-transN)/                  &
               &     (transC-transN),0.0_sp),1.0_sp)
       endif


       if (SED_FLOCS .and. SED_DEFLOC)then
          alpha = MIN(DTsed/t_dfloc,1.0_sp)
          DO k=2,Nbed
             cff3=0.0_sp
             DO ised=1,NST
                cff3=cff3+sedbed%bed_mass(i,k,nnew,ised)
             END DO
             cff1 = 0.0_sp
             DO ised=1,NCS
                cff1 = cff1+sedbed%bed_mass(i,k,nnew,ised)
             END DO
             DO ised=1,NCS
                masseq(ised)=mud_frac_eq(ised)*cff1
                sedbed%bed_mass(i,k,nnew,ised)=alpha*masseq(ised)+           &
                     &             sedbed%bed_mass(i,k,nnew,ised)*(1.0_sp-alpha)
                sedbed%bed_frac(i,k,ised)=sedbed%bed_mass(i,k,nnew,ised)/cff3
             END DO
          END DO
       endif

       if (COHESIVE_BED.or. MIXED_BED)then
          !
          !  Key to this algorithm: "mass available depth" (mad) [kg/m2]
          !  present tau crit profile (tcp)
          !  representative tau crit profile (tcr)          
          !  Compute depth and mass depth of sediment column
          !  Compute cumulative depth
          !  Compute mass depth
          bmz(1) = 0.0_sp
          DO ised=1,NST
             bmz(1) = bmz(1)+sedbed%bed_mass(i,1,nnew,ised)
          ENDDO
          DO k=2,Nbed
             bmz(k) = bmz(k-1)
             DO ised=1,NST
                bmz(k)=bmz(k)+sedbed%bed_mass(i,k,nnew,ised)
             ENDDO
          ENDDO

          !  Calculate representative critical shear stress profile
          !  Note that the values are for the TOP of the layer...we
          !  assume the bottom of the bottom layer has tcr = tcr_max
          tcr(1) = tcr_min
          DO k=2,Nbed
             tcr(k) = tcr_min
             IF (bmz(k-1).GT.eps) THEN
                !             tcr(k) = exp((log(bmz(k-1))-                            &
                !     &            bottom(i,idoff))/                                &
                !     &            bottom(i,idslp))
                tcr(k) = exp((log(bmz(k-1))-                            &
                     &                 tcr_off)/                                 &
                     &                 tcr_slp)
             ENDIF
             tcr(k) = MIN( MAX( tcr(k), tcr_min), tcr_max )
          ENDDO
          !          ks=Kbed-2
          !          DO k=Nbed,Nbed-ks,-1
          !             bed(i,k,ibtcr)=tcr(k)
          !          ENDDO
          !  Relax tau crit profile bottom(i,k,ibtcr) 
          !  towards representative profile tcr...100 x slower for "swelling"
          !  than consolidation
          !  TODO - make the factor an input parameter
          IF( bed(i,1,ibtcr).LE.tcr(1)) THEN
             !          alpha = MIN(DTsed/bottom(i,idtim),1.0_sp)
             alpha = MIN(DTsed/tcr_tim,1.0_sp)
          ELSE
             !         alpha = MIN(DTsed/(100.0_sp*bottom(i,idtim)),1.0_sp)
             alpha = MIN(DTsed/(100.0_sp*tcr_tim),1.0_sp)
          ENDIF
          DO k=1,Nbed
             bed(i,k,ibtcr)=bed(i,k,ibtcr)+                         &
                  &        alpha*(tcr(k)-bed(i,k,ibtcr))
             bed(i,k,ibtcr)=                                          &
                  &       MIN( MAX( bed(i,k,ibtcr), tcr_min), tcr_max )
          ENDDO
          !ALA    Enforce last layer = tcr_max
          !ALA          bed(i,Nbed,ibtcr)= tcr_max

       endif !/* cohesive bed */

    END DO NODE_LOOP3

    !
    !-----------------------------------------------------------------------
    ! Store old bed thickness.
    !-----------------------------------------------------------------------
    !
    if (SED_MORPH)then
       do i=1,m
          sedbed%bed_thick(i,nnew)=0.0_sp
          DO k=1,Nbed
             sedbed%bed_thick(i,nnew)=sedbed%bed_thick(i,nnew)+                  &
                  &                       bed(i,k,ithck)
          END DO
       end do

# if defined (MULTIPROCESSOR)
       if(par)then
          call aexchange(nc,myid,nprocs,sedbed%bed_thick)
       endif
# endif
       bottom(i,morph) = sedbed%bed_thick(i,nnew)
    endif

# if defined (MULTIPROCESSOR)
    if(par)then
       allocate(temp(0:mt,1:nbed))
       DO ised=1,NST
          temp = 0.0_sp
          temp = sedbed%bed_frac(:,:,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_frac(:,:,ised) = temp
          temp = 0.0_sp
          temp = sedbed%bed_mass(:,:,nnew,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_mass(:,:,nnew,ised) = temp
       end do

       do ibed=1,MBEDP
          temp = 0.0_sp      
          temp = bed(:,:,ibed)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          bed(:,:,ibed)=temp
       end do
       deallocate(temp)
    endif
# endif

    if(dbg_set(dbg_sbr)) write(ipt,*) "End:  calc_sed_bed_cohesive"

    return

  End Subroutine calc_sed_bed_cohesive2


!=======================================================================
!                                                                      !
!  This routine computes sediment bed layer stratigraphy.              !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
  Subroutine calc_sed_bed2

    implicit none

    !
    !  Local variable declarations.
    !
    integer :: Ksed, i, ised,ibed, j, k, ks
    integer :: bnew, nnn, cnt,cc,jj

    real(sp), parameter :: eps = 1.0E-14_sp

    real(sp) :: cff, cff1, cff2, cff3
    real(sp) :: thck_avail, thck_to_add

    real(sp), dimension(NST) :: nlysm

    real(sp), dimension(NST) :: dep_mass
    real(sp), dimension(0:mt) :: tau_w

    real(sp),allocatable,dimension(:,:)::temp

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start:  calc_sed_bed2"

    if(BEDLOAD)then
      bnew=nnew
    else
      bnew=nstp
    endif

    ! KLUDGE alert
    ! minlayer_thick = 0.0005
    minlayer_thick = newlayer_thick
    !
    !-----------------------------------------------------------------------
    ! Compute sediment bed layer stratigraphy.
    !-----------------------------------------------------------------------
    !
    if(BEDLOAD_MPM.or.SUSLOAD)then
       tau_w = taub
    endif
    !
    !-----------------------------------------------------------------------
    !  Update bed properties according to ero_flux and dep_flux.
    !-----------------------------------------------------------------------
    !
    if (SUSLOAD)then
       NODE_LOOP : DO i=1,m
          !
          !  The deposition and resuspension of sediment on the bottom "bed"
          !  is due to precipitation flux FC(:,0), already computed, and the
          !  resuspension (erosion, hence called ero_flux). The resuspension is
          !  applied to the bottom-most grid box value qc(:,1) so the total mass
          !  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
          !  after both fluxes are applied.
          !
          SED_LOOP: DO ised=1,NST  
             dep_mass(ised)=0.0_sp
             if (SED_MORPH)then
                ! Apply morphology factor.
                sedbed%ero_flux(i,ised)=sedbed%ero_flux(i,ised)*sed(ised)%morph_fac
                sedbed%settling_flux(i,ised)=sedbed%settling_flux(i,ised)*            &
                     &      sed(ised)%morph_fac
             endif

             !  Update bed mass arrays.
             sedbed%bed_mass(i,1,nnew,ised)=MAX(sedbed%bed_mass(i,1,bnew,ised)-     &
                  &     (sedbed%ero_flux(i,ised)-                                       &
                  &     sedbed%settling_flux(i,ised)),                                  &
                  &        0.0_sp)
             DO k=2,Nbed
                sedbed%bed_mass(i,k,nnew,ised)=sedbed%bed_mass(i,k,nstp,ised)
             END DO
          END DO SED_LOOP

          cff3=0.0_sp
          DO ised=1,NST
             cff3=cff3+sedbed%bed_mass(i,1,nnew,ised)
          END DO

          IF (cff3.eq.0.0_sp) THEN 
             cff3=eps 
          END IF
          bed(i,1,ithck)=0.0_sp
          DO ised=1,NST
             sedbed%bed_frac(i,1,ised)=sedbed%bed_mass(i,1,nnew,ised)/cff3
             bed(i,1,ithck)=MAX(bed(i,1,ithck)+                      &
                  &     sedbed%bed_mass(i,1,nnew,ised)/                                &
                  &     (sed(ised)%Srho*                                         &
                  &     (1.0_sp-bed(i,1,iporo))),0.0_sp)
          END DO

       END DO NODE_LOOP
    endif  !/* SUSPLOAD section */
    !
    !-----------------------------------------------------------------------
    !  At this point, all deposition or erosion is complete, and
    !  has been added/subtracted to top layer. Thickness has NOT been corrected.
    !-----------------------------------------------------------------------
    !
    NODE_LOOP2 : DO i=1,m
       !   Calculate active layer thickness, bottom(i,iactv).
       !   (trunk version allows this to be zero...this has minimum of 6*D50)
       bottom(i,iactv)=MAX(0.0_sp,                                 &
            &             0.007_sp*                               &
            &             (tau_w(i)-bottom(i,itauc)))+   &
            &             6.0_sp*bottom(i,isd50)
       if(SED_MORPH)then
          ! Apply morphology factor.
          bottom(i,iactv)=MAX(bottom(i,iactv)*sed(1)%morph_fac,      &
               &          bottom(i,iactv))
       endif
       !
       !   Calculate net deposition and erosion
       cff=0.0_sp
       cff2=0.0_sp
       DO ised=1,NST
          cff=cff+sedbed%settling_flux(i,ised)
          cff2=cff2+sedbed%ero_flux(i,ised)
          dep_mass(ised)=0.0_sp
          IF ((sedbed%ero_flux(i,ised)-sedbed%settling_flux(i,ised)).lt.      &
               &     0.0_sp) THEN
             dep_mass(ised)=sedbed%settling_flux(i,ised)-              &
                  &        sedbed%ero_flux(i,ised)
          END IF
       END DO

       IF ( cff-cff2.GT.0.0_sp) THEN ! NET depostion
          !  Deposition. Determine if we need to create a new bed layer 
          ! (no test for age here)
          bed(i,1,iaged)=T_model
          IF(bed(i,1,ithck).gt.                                   &
               &    MAX(bottom(i,iactv),newlayer_thick)) THEN
             ! Top layer is too thick
             IF (Nbed.gt.2) THEN
                IF(bed(i,2,ithck).lt.minlayer_thick) THEN
                   ! Layer 2 is smaller than minimum size
                   ! Instead of pushing down all layers, just combine top 2 layers
                   cff=0.0_sp
                   cff1=0.0_sp
                   cff2=0.0_sp
                   DO ised=1,NST
                      cff =cff +dep_mass(ised)
                      cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                      cff2=cff2+sedbed%bed_mass(i,2,nnew,ised)
                   END DO
                   !  Update bed mass
                   DO ised=1,NST
                      sedbed%bed_mass(i,2,nnew,ised)=                    &
                           &                 MAX(sedbed%bed_mass(i,2,nnew,ised)+           &
                           &                 sedbed%bed_mass(i,1,nnew,ised)-               &
                           &                 dep_mass(ised),0.0_sp)
                      sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                   END DO
                   ! ALA - average T_model and porosity
                   ! ALA CHECK WITH CRS cff1 or cff1-cff for first layer
                   bed(i,2,iaged)=(bed(i,1,iaged)*cff1+         &
                        &               bed(i,2,iaged)*cff2)/(cff1+cff2)
                   bed(i,1,iaged)=T_model
                   bed(i,2,iporo)=(bed(i,1,iporo)*cff1+         &
                        &               bed(i,2,iporo)*cff2)/(cff1+cff2)
                   ! ALA CHECK WITH CRS POROSITY OF 1ST LAYER
                   bed(i,1,iporo)=bed(i,1,iporo)
                ELSE
                   ! Layer 2 is > minlayer thick, need another layer
                   !  Combine bottom layers.
                   cff1=0.0_sp
                   cff2=0.0_sp
                   DO ised=1,NST
                      cff1=cff1+sedbed%bed_mass(i,Nbed-1,nnew,ised)
                      cff2=cff2+sedbed%bed_mass(i,Nbed,nnew,ised)
                   END DO
                   bed(i,Nbed,iporo)=                             &
                        &              (bed(i,Nbed-1,iporo)*cff1+                &
                        &              bed(i,Nbed,iporo)*cff2)/(cff1+cff2)
                   bed(i,Nbed,iaged)=                             &
                        &              (bed(i,Nbed-1,iaged)*cff1+                &
                        &               bed(i,Nbed,iaged)*cff2)/(cff1+cff2)
                   DO ised=1,NST
                      sedbed%bed_mass(i,Nbed,nnew,ised)=                 &
                           &            sedbed%bed_mass(i,Nbed-1,nnew,ised)+          &
                           &            sedbed%bed_mass(i,Nbed  ,nnew,ised)
                   END DO
                   !
                   !  Push layers down.
                   DO k=Nbed-1,2,-1
                      bed(i,k,iporo)=bed(i,k-1,iporo)
                      bed(i,k,iaged)=bed(i,k-1,iaged)
                      DO ised =1,NST
                         sedbed%bed_mass(i,k,nnew,ised)=                 &
                              &               sedbed%bed_mass(i,k-1,nnew,ised)
                      END DO
                   END DO
                   !  Set new top parameters for top 2 layers
                   DO ised=1,NST
                      sedbed%bed_mass(i,2,nnew,ised)=                    &
                           &              MAX(sedbed%bed_mass(i,2,nnew,ised)-           &
                           &              dep_mass(ised),0.0_sp)
                      sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                   END DO
                END IF
             ELSEIF (Nbed.eq.2) THEN 
                ! NBED=2
                cff1=0.0_sp
                cff2=0.0_sp
                DO ised=1,NST
                   cff1=cff1+sedbed%bed_mass(i,1,nnew,ised)
                   cff2=cff2+sedbed%bed_mass(i,2,nnew,ised)
                END DO
                DO ised=1,NST
                   sedbed%bed_mass(i,2,nnew,ised)=                       &
                        &            MAX(sedbed%bed_mass(i,2,nnew,ised)+              &
                        &                sedbed%bed_mass(i,1,nnew,ised)-              &
                        &                dep_mass(ised),0.0_sp)
                   sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                END DO
                ! ALA - average T_model and porosity
                bed(i,2,iaged)=(bed(i,1,iaged)*cff1+            &
                     &              bed(i,2,iaged)*cff2)/(cff1+cff2)
                bed(i,1,iaged)=T_model
                bed(i,2,iporo)=(bed(i,1,iporo)*cff1+            &
                     &              bed(i,2,iporo)*cff2)/(cff1+cff2)
                ! ALA CHECK WITH CRS POROSITY OF 1ST LAYER
                bed(i,1,iporo)=bed(i,1,iporo)                    
             ELSE
                ! NBED=1
             END IF
          ELSE
             ! Net deposition has occured, but no new bed layer was created 
          END IF
       ELSE
          ! Net erosion occurred
          bed(i,1,iaged)=T_model
          IF (Nbed.eq.2) THEN 
             ! NBED=2
             DO ised=1,NST
                sedbed%bed_mass(i,2,nnew,ised)=                       &
                     &           MAX(sedbed%bed_mass(i,2,nnew,ised)+              &
                     &               sedbed%bed_mass(i,1,nnew,ised)-                  &
                     &               dep_mass(ised),0.0_sp)
                sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
             END DO
          ELSEIF (Nbed.eq.1) THEN
             ! ALF NO NEED TO DO ANYTHING
          ELSE
          END IF
       END IF

       ! Recalculate thickness and fractions for all layers.
       DO k=1,Nbed
          cff3=0.0_sp
          DO ised=1,NST
             cff3=cff3+sedbed%bed_mass(i,k,nnew,ised)
          END DO
          IF (cff3.eq.0.0_sp) THEN 
             cff3=eps 
          END IF
          bed(i,k,ithck)=0.0_sp
          DO ised=1,NST
             sedbed%bed_frac(i,k,ised)=sedbed%bed_mass(i,k,nnew,ised)/cff3
             bed(i,k,ithck)=MAX(bed(i,k,ithck)+                 &
                  &          sedbed%bed_mass(i,k,nnew,ised)/               &
                  &             (sed(ised)%Srho*                          &
                  &          (1.0_sp-bed(i,k,iporo))),0.0_sp)
          END DO
       END DO

    END DO NODE_LOOP2


    NODE_LOOP3 : DO i=1,m
       IF (bottom(i,iactv).gt.bed(i,1,ithck)) THEN
          IF (Nbed.eq.1) THEN
             bottom(i,iactv)=bed(i,1,ithck)
          ELSE
             thck_to_add=bottom(i,iactv)-bed(i,1,ithck)
             thck_avail=0.0_sp
             Ksed=1                                        ! initialize
             DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                   thck_avail=thck_avail+bed(i,k,ithck)
                   Ksed=k
                END IF
             END DO
             !
             !  Catch here if there was not enough bed material.
             !
             IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,iactv)=bed(i,1,ithck)+thck_avail
                thck_to_add=thck_avail
             END IF
             !
             !  Update bed mass of top layer and fractional layer.
             !
             cff2=MAX(thck_avail-thck_to_add,0.0_sp)/                  &
                  &             MAX(bed(i,Ksed,ithck),eps)	        
             DO ised=1,NST
                cff1=0.0_sp
                DO k=1,Ksed
                   cff1=cff1+sedbed%bed_mass(i,k,nnew,ised)
                END DO
                cff3=cff2*sedbed%bed_mass(i,Ksed,nnew,ised)
                sedbed%bed_mass(i,1   ,nnew,ised)=cff1-cff3
                sedbed%bed_mass(i,Ksed,nnew,ised)=cff3
             END DO
             !
             !  Update thickness of fractional layer ksource_sed.
             !
             bed(i,Ksed,ithck)=MAX(thck_avail-thck_to_add,0.0_sp)
             !
             !  Update bed fraction of top layer.
             !
             cff3=0.0_sp
             DO ised=1,NST
                cff3=cff3+sedbed%bed_mass(i,1,nnew,ised)
             END DO
             IF (cff3.eq.0.0_sp) THEN
                cff3=eps
             END IF
             DO ised=1,NST
                sedbed%bed_frac(i,1,ised)=sedbed%bed_mass(i,1,nnew,ised)/cff3
             END DO
             !
             !  Upate bed thickness of top layer.
             !
             bed(i,1,ithck)=bottom(i,iactv)
             !
             !  Pull all layers closer to the surface.
             !
             DO k=Ksed,Nbed
                ks=Ksed-2
                bed(i,k-ks,ithck)=bed(i,k,ithck)
                bed(i,k-ks,iporo)=bed(i,k,iporo)
                bed(i,k-ks,iaged)=bed(i,k,iaged)
                DO ised=1,NST
                   sedbed%bed_frac(i,k-ks,ised)=sedbed%bed_frac(i,k,ised)
                   sedbed%bed_mass(i,k-ks,nnew,ised)=sedbed%bed_mass(i,k,nnew,ised)
                END DO
             END DO
             !
             !  Add new layers onto the bottom. Split what was in the bottom layer to
             !  fill these new empty cells. ("ks" is the number of new layers).
             !
             ks=Ksed-2
             ! ALA CHECK WITH CRS about bed_frac 
             nnn=0
             DO ised=1,NST
                nlysm(ised)=newlayer_thick*REAL(ks+1,sp)*          &
                     &  (sed(ised)%Srho*                                 &
                     &  (1.0_sp-bed(i,Nbed-ks,iporo)))*               &
                     &  sedbed%bed_frac(i,Nbed-ks,ised)
                IF (ks.gt.0) THEN
                   IF (sedbed%bed_mass(i,Nbed-ks,nnew,ised).gt.         &
                        & nlysm(ised)) THEN
                      nnn=nnn+1
                      nlysm(ised)=                                 &
                           &  newlayer_thick*REAL(ks,sp)*        &
                           & (sed(ised)%Srho*                         &
                           &  (1.0_sp-bed(i,Nbed-ks,iporo)))*       &
                           &  sedbed%bed_frac(i,Nbed-ks,ised)
                   END IF
                END IF
             END DO
             IF (nnn.eq.NST) THEN
                bed(i,Nbed,ithck)=bed(i,Nbed-ks,ithck)-           &
                     &  newlayer_thick*REAL(ks,sp)
                DO ised=1,NST
                   sedbed%bed_mass(i,Nbed,nnew,ised)=                      &
                        &  sedbed%bed_mass(i,Nbed-ks,nnew,ised)-nlysm(ised)
                END DO
                DO k=Nbed-1,Nbed-ks,-1
                   bed(i,k,ithck)=newlayer_thick
                   bed(i,k,iaged)=bed(i,Nbed-ks,iaged)
                   DO ised=1,NST
                      sedbed%bed_frac(i,k,ised)=sedbed%bed_frac(i,Nbed-ks,ised)
                      sedbed%bed_mass(i,k,nnew,ised)=                      &
                           &  nlysm(ised)/REAL(ks,sp)
                   END DO
                END DO
             ELSE
                cff=1.0_sp/REAL(ks+1,sp)
                DO k=Nbed,Nbed-ks,-1
                   bed(i,k,ithck)=bed(i,Nbed-ks,ithck)*cff
                   bed(i,k,iaged)=bed(i,Nbed-ks,iaged)
                   DO ised=1,NST
                      sedbed%bed_frac(i,k,ised)=sedbed%bed_frac(i,Nbed-ks,ised)
                      sedbed%bed_mass(i,k,nnew,ised)=                      &
                           &  sedbed%bed_mass(i,Nbed-ks,nnew,ised)*cff
                   END DO
                END DO
             END IF
          END IF  ! Nbed > 1
       END IF  ! increase top bed layer
    END DO NODE_LOOP3
    !
    !-----------------------------------------------------------------------
    ! Store old bed thickness.
    !-----------------------------------------------------------------------
    !
    if (SED_MORPH)then
       do i=1,m
          sedbed%bed_thick(i,nnew)=0.0_sp
          DO k=1,Nbed
             sedbed%bed_thick(i,nnew)=sedbed%bed_thick(i,nnew)+                  &
                  &                       bed(i,k,ithck)
          END DO
       end do

# if defined (MULTIPROCESSOR)
       if(par)then
          call aexchange(nc,myid,nprocs,sedbed%bed_thick)
       endif
# endif
       bottom(i,morph) = sedbed%bed_thick(i,nnew)
    endif
    !
    !-----------------------------------------------------------------------
    !  Apply periodic or gradient boundary conditions to property arrays.
    !-----------------------------------------------------------------------
    !

# if defined (MULTIPROCESSOR)
    if(par)then
       allocate(temp(0:mt,1:nbed))
       DO ised=1,NST
          temp = 0.0_sp
          temp = sedbed%bed_frac(:,:,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_frac(:,:,ised) = temp
          temp = 0.0_sp
          temp = sedbed%bed_mass(:,:,nnew,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_mass(:,:,nnew,ised) = temp
       end do

       do ibed=1,MBEDP
          temp = 0.0_sp      
          temp = bed(:,:,ibed)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          bed(:,:,ibed)=temp
       end do
       deallocate(temp)
    endif
# endif


    if(dbg_set(dbg_sbr)) write(ipt,*) "End:  calc_sed_bed2"

    return

  End Subroutine calc_sed_bed2




!=======================================================================
!                                                                      !
!  This routine computes sediment bed layer stratigraphy.              !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
  Subroutine calc_sed_bed

    implicit none
    !
    !  Local variable declarations.
    !
    integer :: Ksed, i, ised,ibed, j, k, ks
    integer :: bnew, nnn, cnt,cc,jj

    real(sp), parameter :: eps = 1.0E-14_sp

    real(sp) :: cff, cff1, cff2, cff3
    real(sp) :: thck_avail, thck_to_add

    real(sp), dimension(NST) :: dep_mass
    real(sp), dimension(0:mt) :: tau_w

    real(sp),allocatable,dimension(:,:)::temp

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start:  calc_sed_bed"


    if(BEDLOAD)then
      bnew=nnew
    else
      bnew=nstp
    endif

    !
    !-----------------------------------------------------------------------
    ! Compute sediment bed layer stratigraphy.
    !-----------------------------------------------------------------------
    !
    if(BEDLOAD_MPM.or.SUSLOAD)then
       tau_w = taub
    endif
    !
    !-----------------------------------------------------------------------
    !  Update bed properties according to ero_flux and dep_flux.
    !-----------------------------------------------------------------------
    !
    if(SUSLOAD)then
       NODE_LOOP : DO i=1,m
          SED_LOOP: DO ised=1,NST
             !
             !  The deposition and resuspension of sediment on the bottom "bed"
             !  is due to precepitation flux FC(:,0), already computed, and the
             !  resuspension (erosion, hence called ero_flux). The resuspension is
             !  applied to the bottom-most grid box value qc(:,1) so the total mass
             !  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
             !  after both fluxes are applied.
             !
             dep_mass(ised)=0.0_sp

             if( SED_MORPH)then
                !
                ! Apply morphology factor.
                !
                sedbed%ero_flux(i,ised)=sedbed%ero_flux(i,ised)*sed(ised)%morph_fac
                sedbed%settling_flux(i,ised)=sedbed%settling_flux(i,ised)*            &
                     &                              sed(ised)%morph_fac
                !
             endif
             IF ((sedbed%ero_flux(i,ised)-sedbed%settling_flux(i,ised)).lt.        &
                  &           0.0_sp) THEN
                !
                !  If first time step of deposit, then store deposit material in
                !  temporary array, dep_mass.
                !
                IF ((T_model.gt.(bed(i,1,iaged)+1.1_sp*DTsed)).and.   &
                     &            (bed(i,1,ithck).gt.newlayer_thick)) THEN
                   dep_mass(ised)=sedbed%settling_flux(i,ised)-               &
                        &                           sedbed%ero_flux(i,ised)
                END IF
                bed(i,1,iaged)=T_model
             END IF
             !
             !  Update bed mass arrays.
             !
             sedbed%bed_mass(i,1,nnew,ised)=MAX(sedbed%bed_mass(i,1,bnew,ised)-    &
                  &  (sedbed%ero_flux(i,ised)-           &
                  &  sedbed%settling_flux(i,ised)),     &
                  &  0.0_sp)
             DO k=2,Nbed
                sedbed%bed_mass(i,k,nnew,ised)=sedbed%bed_mass(i,k,nstp,ised)
             END DO
          END DO SED_LOOP
          !
          !  If first time step of deposit, create new layer and combine bottom
          !  two bed layers.
          !
          cff=0.0_sp
          !
          !  Determine if deposition ocurred here.
          !
          IF (Nbed.gt.1) THEN
             DO ised=1,NST
                cff=cff+dep_mass(ised)
             END DO
             IF (cff.gt.0.0_sp) THEN
                !
                !  Combine bottom layers.
                !
                bed(i,Nbed,iporo)=0.5_sp*(bed(i,Nbed-1,iporo)+        &
                     &                                    bed(i,Nbed,iporo))
                bed(i,Nbed,iaged)=0.5_sp*(bed(i,Nbed-1,iaged)+        &
                     &                                    bed(i,Nbed,iaged))
                DO ised=1,NST
                   sedbed%bed_mass(i,Nbed,nnew,ised)=                           &
                        &  sedbed%bed_mass(i,Nbed-1,nnew,ised)+      &
                        &  sedbed%bed_mass(i,Nbed  ,nnew,ised)
                END DO
                !
                !  Push layers down.
                !
                DO k=Nbed-1,2,-1
                   bed(i,k,iporo)=bed(i,k-1,iporo)
                   bed(i,k,iaged)=bed(i,k-1,iaged)
                   DO ised =1,NST
                      sedbed%bed_mass(i,k,nnew,ised)=sedbed%bed_mass(i,k-1,nnew,ised)
                   END DO
                END DO
                !
                !  Set new top layer parameters.
                !
                DO ised=1,NST
                   sedbed%bed_mass(i,2,nnew,ised)=MAX(sedbed%bed_mass(i,2,nnew,ised)-&
                        &  dep_mass(ised),0.0_sp)
                   sedbed%bed_mass(i,1,nnew,ised)=dep_mass(ised)
                END DO
             END IF
          END IF !NBED=1
          !
          ! Recalculate thickness and fractions for all layers.
          !
          DO k=1,Nbed
             cff3=0.0_sp
             DO ised=1,NST
                cff3=cff3+sedbed%bed_mass(i,k,nnew,ised)
             END DO
             IF (cff3.eq.0.0_sp) THEN
                cff3=eps
             END IF
             bed(i,k,ithck)=0.0_sp
             DO ised=1,NST
                sedbed%bed_frac(i,k,ised)=sedbed%bed_mass(i,k,nnew,ised)/cff3
                bed(i,k,ithck)=MAX(bed(i,k,ithck)+                    &
                     &  sedbed%bed_mass(i,k,nnew,ised)/               &
                     &  (sed(ised)%Srho*                          &
                     &  (1.0_sp-bed(i,k,iporo))),0.0_sp)
             END DO
          END DO
       END DO NODE_LOOP
       !
       !  End of Suspended Sediment only section.
       !
    endif
    !
    !  Ensure top bed layer thickness is greater or equal than active layer
    !  thickness. If need to add sed to top layer, then entrain from lower
    !  levels. Create new layers at bottom to maintain Nbed.
    !
    NODE_LOOP2 : DO i=1,m
       !
       !  Calculate active layer thickness, bottom(i,iactv).
       !
       bottom(i,iactv)=MAX(0.0_sp,                                 &
            &   0.007_sp*                               &
            &   (tau_w(i)-bottom(i,itauc)))+   &
            &    6.0_sp*bottom(i,isd50)
       !
       if (SED_MORPH)then
          !
          ! Apply morphology factor.
          !
          bottom(i,iactv)=MAX(bottom(i,iactv)*sed(1)%morph_fac,      &
               &   bottom(i,iactv))
       endif
       !
       IF (bottom(i,iactv).gt.bed(i,1,ithck)) THEN
          IF (Nbed.eq.1) THEN
             bottom(i,iactv)=bed(i,1,ithck)
          ELSE
             thck_to_add=bottom(i,iactv)-bed(i,1,ithck)
             thck_avail=0.0_sp
             Ksed=1                                        ! initialize
             DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                   thck_avail=thck_avail+bed(i,k,ithck)
                   Ksed=k
                END IF
             END DO
             !
             !  Catch here if there was not enough bed material.
             !
             IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,iactv)=bed(i,1,ithck)+thck_avail
                thck_to_add=thck_avail
             END IF
             !
             !  Update bed mass of top layer and fractional layer.
             !
             cff2=MAX(thck_avail-thck_to_add,0.0_sp)/                  &
                  &   MAX(bed(i,Ksed,ithck),eps)
             DO ised=1,NST
                cff1=0.0_sp
                DO k=1,Ksed
                   cff1=cff1+sedbed%bed_mass(i,k,nnew,ised)
                END DO
                cff3=cff2*sedbed%bed_mass(i,Ksed,nnew,ised)
                sedbed%bed_mass(i,1   ,nnew,ised)=cff1-cff3
                sedbed%bed_mass(i,Ksed,nnew,ised)=cff3
             END DO
             !
             !  Update thickness of fractional layer ksource_sed.
             !
             bed(i,Ksed,ithck)=MAX(thck_avail-thck_to_add,0.0_sp)
             !
             !  Upate bed fraction of top layer.
             !
             cff3=0.0_sp
             DO ised=1,NST
                cff3=cff3+sedbed%bed_mass(i,1,nnew,ised)
             END DO
             IF (cff3.eq.0.0_sp) THEN
                cff3=eps
             END IF
             DO ised=1,NST
                sedbed%bed_frac(i,1,ised)=sedbed%bed_mass(i,1,nnew,ised)/cff3
             END DO
             !
             !  Upate bed thickness of top layer.
             !
             bed(i,1,ithck)=bottom(i,iactv)
             !
             !  Pull all layers closer to the surface.
             !
             DO k=Ksed,Nbed
                ks=Ksed-2
                bed(i,k-ks,ithck)=bed(i,k,ithck)
                bed(i,k-ks,iporo)=bed(i,k,iporo)
                bed(i,k-ks,iaged)=bed(i,k,iaged)
                DO ised=1,NST
                   sedbed%bed_frac(i,k-ks,ised)=sedbed%bed_frac(i,k,ised)
                   sedbed%bed_mass(i,k-ks,nnew,ised)=sedbed%bed_mass(i,k,nnew,ised)
                END DO
             END DO
             !
             !  Add new layers onto the bottom. Split what was in the bottom layer to
             !  fill these new empty cells. ("ks" is the number of new layers).
             !
             ks=Ksed-2
             cff=1.0_sp/REAL(ks+1,sp)
             DO k=Nbed,Nbed-ks,-1
                bed(i,k,ithck)=bed(i,Nbed-ks,ithck)*cff
                bed(i,k,iaged)=bed(i,Nbed-ks,iaged)
                DO ised=1,NST
                   sedbed%bed_frac(i,k,ised)=sedbed%bed_frac(i,Nbed-ks,ised)
                   sedbed%bed_mass(i,k,nnew,ised)=                            &
                        &  sedbed%bed_mass(i,Nbed-ks,nnew,ised)*cff
                END DO
             END DO
          END IF  ! Nbed > 1
       END IF  ! increase top bed layer

    END DO NODE_LOOP2
    !
    !-----------------------------------------------------------------------
    ! Store old bed thickness.
    !-----------------------------------------------------------------------
    !
    if (SED_MORPH)then
       do i=1,m
          sedbed%bed_thick(i,nnew)=0.0_sp
          DO k=1,Nbed
             sedbed%bed_thick(i,nnew)=sedbed%bed_thick(i,nnew)+                  &
                  &                       bed(i,k,ithck)
          END DO
       end do

# if defined (MULTIPROCESSOR)
       if(par)then
          call aexchange(nc,myid,nprocs,sedbed%bed_thick)
       endif
# endif
       bottom(i,morph) = sedbed%bed_thick(i,nnew)
    endif





# if defined (MULTIPROCESSOR)
    if(par)then
       allocate(temp(0:mt,1:nbed))
       DO ised=1,NST
          temp = 0.0_sp
          temp = sedbed%bed_frac(:,:,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_frac(:,:,ised) = temp
          temp = 0.0_sp
          temp = sedbed%bed_mass(:,:,nnew,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_mass(:,:,nnew,ised) = temp
       end do

       do ibed=1,MBEDP
          temp = 0.0_sp      
          temp = bed(:,:,ibed)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          bed(:,:,ibed)=temp
       end do
       deallocate(temp)
    endif
# endif

    if(dbg_set(dbg_sbr)) write(ipt,*) "End:  calc_sed_bed"

    return   

  End Subroutine calc_sed_bed

!=======================================================================
!                                                                      !
!  This routine computes sediment bed layer mixing (biodiffusion)      !
!                                                                      !
!                                                                      !
!=======================================================================
  Subroutine calc_sed_biodiff

    implicit none

 


    !
    !  Local variable declarations.
    !
    integer :: i, j, k, ised,ibed
    real(sp), parameter :: eps = 1.0E-14_sp
    real(sp), parameter :: cmy2ms = 3.1688765E-12_sp ! multiply cm2/yr by this to get m2/s

    real(sp) :: cff, cff1, cff2, cff3

    real(sp), dimension(0:kb) :: dep_mass

    integer :: iu,il,lp,ii
    !      real(sp) :: rtemp, Dbmx, Dbmm, zs, zm, zp
    real(sp) :: rtemp
    real(sp), dimension(Nbed) :: zb
    real(sp), dimension(Nbed) :: zc
    real(sp), dimension(Nbed) :: dzui
    real(sp), dimension(Nbed) :: dzli
    real(sp), dimension(Nbed) :: dzmi
    real(sp), dimension(Nbed) :: Db
    real(sp), dimension(Nbed) :: Dc
    real(sp), dimension(Nbed) :: a
    real(sp), dimension(Nbed) :: d
    real(sp), dimension(Nbed) :: b
    real(sp), dimension(Nbed) :: cc
    real(sp), dimension(Nbed) :: dd
    real(sp),allocatable,dimension(:,:)::temp

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: calc_sed_biodiff "
    !
    !-----------------------------------------------------------------------
    ! Compute sediment bed layer stratigraphy.
    !-----------------------------------------------------------------------
    !
    !      print *, 'Mixing...'
    NODE_LOOP : DO i=1,m
       !
       !  Set mixing coefficient profile
       !  (hardwire uniform mixing)
       !
       !          DO k=1,Nbed
       !             Db(k)=1.0E-8_sp
       !          ENDDO

       IF (Nbed.GT.2) THEN
          !  Compute cumulative depth 
          zb(1)=bed(i,1,ithck)
          DO k=2,Nbed
             zb(k)=zb(k-1)+bed(i,k,ithck)
          END DO
          !  Compute depths of bed centers 
          zc(1)=0.5_sp*(bed(i,1,ithck))
          DO k=2,Nbed
             zc(k)=zb(k-1)+0.5_sp*bed(i,k,ithck)
          END DO
          !
          !  Set biodiffusivity profile
          if (DB_PROFILE)then
             !  Depth-varying biodiffusivity profile
             !   
             !ALF            Dbmx = bottom(i,idbmx)
             !ALF            Dbmm = bottom(i,idbmm)
             !ALF            zs = bottom(i,idbzs)
             !ALF            zm = bottom(i,idbzm)
             !ALF            zp = bottom(i,idbzp)
             DO k=1,Nbed
                Db(k)= Dbmx
                !              IF( zb(k).GT.Dbzp)THEN          ! should be .GE. ?
                IF( zb(k).GE.Dbzp)THEN          
                   Db(k)=0.0_sp
                ELSE
                   IF((zb(k) .GT. Dbzs).AND.                                 &
                        &           (zb(k) .LE. Dbzm)) THEN
                      rtemp= LOG(Dbmm/Dbmx)/(-Dbzm-Dbzs)
                      Db(k)=Dbmx*exp(rtemp*(-zb(k)-Dbzs))
                   ELSEIF((zb(k).GT.Dbzm).AND.                               &
                        &           (zb(k).LT.Dbzp) ) THEN
                      Db(k)=(Dbmm-(Dbmm/(Dbzp-Dbzm)))*            &
                           &                  (zb(k)-Dbzm)
                   ENDIF
                ENDIF
             END DO
          else
             !    Uniform biodiffusivity profile at max value 
             !
             DO k=1,Nbed
                Db(k)=Dbmx
                !              Db(k)=bottom(i,idbmx)
             ENDDO
          endif  !/* defined DB_PROFILE */
          !           write(*,*) 'Db',Db
          !  Calculate finite differences
          dzui(1)=1.0_sp/(zc(2)-zc(1))
          dzli(1)=1.E35_sp       ! should not be needed
          dzmi(1)=1.0_sp/bed(i,1,ithck)
          DO k=2,Nbed-1
             dzui(k)=1.0_sp/(zc(k+1)-zc(k))
             dzli(k)=1.0_sp/(zc(k)-zc(k-1))
             !dzmi(k)=1.0_sp/(zb(k+1)-zb(k))
             ! equivalent:
             dzmi(k)=1.0/bed(i,k,ithck)
             !  Tridiagonal terms
             b(k)= -dtsed*dzmi(k)*Db(k-1)*dzli(k)
             d(k)=1.0_sp+dtsed*dzmi(k)*                               &
                  &          ( Db(k-1)*dzli(k)+Db(k)*dzui(k) )
             a(k)= -dtsed*dzmi(k)*Db(k)*dzui(k)
          ENDDO
          dzui(Nbed)=1.0E35_sp ! should not be needed
          dzli(Nbed)=1.0_sp/(zc(Nbed)-zc(Nbed-1))
          dzmi(Nbed)=1.0_sp/bed(i,Nbed,ithck)
          !  No-flux boundary conditions
          b(1) = 999.9_sp  ! should not be needed
          d(1)=   1.0_sp +dtsed*dzmi(1)*Db(1)*dzui(1)
          a(1)=      -dtsed*dzmi(1)*Db(1)*dzui(1)
          b(Nbed)=   -dtsed*dzmi(Nbed)*Db(Nbed-1)*dzli(Nbed)
          d(Nbed)=1.0_sp+dtsed*dzmi(Nbed)*Db(Nbed-1)*dzli(Nbed)
          a(Nbed)=999.9_sp ! should not be needed
          !
          !   Calculate mixing for each size fraction
          DO ised=1,NST
             !   ...make working copies 
             DO k=1,Nbed
                cc(k) = sedbed%bed_frac(i,k,ised)
                dd(k)= d(k)
             ENDDO
             !   Solve a tridiagonal system of equations using Thomas' algorithm
             !   Anderson, Tannehill, and Pletcher (1984) pp. 549-550
             !   ...establish upper triangular matrix
             il = 1
             iu = Nbed
             lp = il+1
             DO k = lp,iu
                rtemp = b(k)/dd(k-1)
                dd(k)= dd(k)-rtemp*a(k-1);
                cc(k)= cc(k)-rtemp*cc(k-1);
             ENDDO
             !   ...back substitution
             cc(iu) = cc(iu)/dd(iu)
             DO k  = lp,iu
                ii = iu-k+il;
                cc(ii) = (cc(ii)-a(ii)*cc(ii+1))/dd(ii);
             ENDDO
             !   ...solution stored in cc; copy out
             DO k = 1,Nbed
                sedbed%bed_frac(i,k,ised)=cc(k)
             ENDDO
          ENDDO
          ! TODO - Mix porosity or assign it as f(depth)?
          ! TODO - Mix age?

          ! Recompute bed masses
          DO k=1,Nbed
             !             write(*,*) i,k,(bed_frac(i,k,ised),ised=1,NST)
             !            debugging: ensure fracs add up to 1
             cff3 = 0.0_sp
             DO ised=1,NST
                cff3 = cff3+sedbed%bed_frac(i,k,ised)
             ENDDO
             if( abs(1.0_sp-cff3).GT.1e-6 )                            &
                  &         write(*,*) 'error: sum_frac: ',cff3
             cff3=0.0_sp
             DO ised=1,NST
                cff3=cff3+sedbed%bed_mass(i,k,nnew,ised)
             ENDDO
             DO ised=1,NST
                sedbed%bed_mass(i,k,nnew,ised)=sedbed%bed_frac(i,k,ised)*cff3
             ENDDO
          ENDDO


       END IF !NBED.GT.2

    END DO NODE_LOOP

# if defined (MULTIPROCESSOR)
    if(par)then
       allocate(temp(0:mt,1:nbed))
       DO ised=1,NST
          temp = 0.0_sp
          temp = sedbed%bed_frac(:,:,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_frac(:,:,ised) = temp
          temp = 0.0_sp
          temp = sedbed%bed_mass(:,:,nnew,ised)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          sedbed%bed_mass(:,:,nnew,ised) = temp
       end do

       do ibed=1,MBEDP
          temp = 0.0_sp      
          temp = bed(:,:,ibed)
          call node_match(1,nbn,bn_mlt,bn_loc,bnc,mt,nbed,myid,nprocs,temp)
          call aexchange(nc,myid,nprocs,temp)
          bed(:,:,ibed)=temp
       end do
       deallocate(temp)
    endif
# endif

    if(dbg_set(dbg_sbr)) write(ipt,*) "End:  calc_sed_biodiff"

    return

  End Subroutine calc_sed_biodiff


!
!-----------------------------------------------------------------------
!  Horizontal Advection of Sediment Concentration 
!-----------------------------------------------------------------------
!
  Subroutine calc_fluid_mud

    implicit none
    integer :: i

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start: calc_fluid_mud "

# if defined (FLUID_MUD)


    call fluid_mud_source_sink(DTsed,sed(1)%Srho,bottom(:,bflag)&
         &,bed(:,1,ithck),bed(:,1,iporo),sedbed%frac(:,1,1),taub,sed(1)%t_ce,sed(1)%rate)

    call internal_step_fluid_mud

    do ised=1,nst
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

    if(dbg_set(dbg_sbr)) write(ipt,*) "End: calc_fluid_mud "

    Return

  End Subroutine calc_fluid_mud
!
!-----------------------------------------------------------------------
!  Horizontal Advection of Sediment Concentration 
!-----------------------------------------------------------------------
!
  Subroutine Calc_Sed_Advection
    Use Scalar
    implicit none
    integer :: i,k,ised,d_cdis,d_cflx
    if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Sed_Advection "
    if(.not. oned_model)then
       d_cdis = 1
       if(numqbc > 0) d_cdis = numqbc
       d_cflx = 1
       if(Nobc > 0) d_cflx = Nobc

       do i=1,nst

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
       do i=1,nst
          sed(i)%cnew = sed(i)%conc
       end do
    endif

    if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Sed_Advection "

    return

  End Subroutine Calc_Sed_Advection
  !
  !------------------------------------------------------
  !Vertical Diffusion of Sediment Concentrations
  !------------------------------------------------------
  !
  Subroutine Calc_Sed_Diffusion
    Use Scalar
    implicit none
    integer :: i

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Sed_Diffusion "

    !Jianzhong GE 03/05/2013
    if(vert_hindered)call hindered_diffusion
    !Jianzhong GE 03/05/2013

    do i=1,nst
       call vdif_scal(sed(i)%cnew,DTsed) 
    end do

    if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Sed_Diffusion "

    return

  End Subroutine Calc_Sed_Diffusion


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
    do ii=1,nst
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

    return

  End Subroutine Hindered_Diffusion
!Jianzhong Ge 03/05/2013



  !------------------------------------------------------
  ! Set Boundary Conditions
  !------------------------------------------------------
  Subroutine Calc_Sed_Bcond
    Use Scalar
    implicit none
    integer :: i
    real(sp) :: temp

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Sed_Bcond"

    !------------------------------------------------------
    !Exchange Concentration on Processor Boundaries
    !------------------------------------------------------
# if defined (MULTIPROCESSOR)
    if(par)then
       do i=1,nst
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
       do i=1,nst 
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
       do i=1,nst
          call bcond_scal_PTsource(sed(i)%conc, &
               sed(i)%cnew, &
               sed(i)%cdis) 
       end do
    endif

    if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Sed_Bcond"

    return

  End Subroutine Calc_Sed_Bcond




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
        bottom(i,tmass) = bottom(i,tmass) + sedbed%bed_mass(i,k,nnew,ised)
      end do
    end do
  end do
  
  End Subroutine Update_Thickness_Delta




  Subroutine finalize_sed_step
    implicit none
    integer :: i,ised

    if(dbg_set(dbg_sbr)) write(ipt,*) "Start: finalize_sed_step"

    !------------------------------------------------------
    ! Update Concentration Variables to next time level
    !------------------------------------------------------
    do ised=1,nst
       sed(ised)%conc = sed(ised)%cnew
    end do

    !------------------------------------------------------
    ! store the fraction at bed surface
    !------------------------------------------------------
    do ised=1,nst
       sed(ised)%frac(1:m,1:nbed)=sedbed%bed_frac(1:m,1:nbed,ised)
    end do

    !------------------------------------------------------
    ! Limit to Non-Negative (modify to Venkat Limiter) 
    !------------------------------------------------------
    do i=1,nst
       sed(i)%conc = max(sed(i)%cnew,0.0)
    end do

    !------------------------------------------------------
    ! Calculate total sediment concentration (all classes)
    !------------------------------------------------------
    csed = 0.0
    do i=1,nst
       csed = csed + sed(i)%conc
    end do

    !------------------------------------------------------
    !Exchange Concentration on Processor Boundaries
    !------------------------------------------------------
# if defined (MULTIPROCESSOR)
    if(par)then
       do i=1,nst
          call aexchange(nc,myid,nprocs,sed(i)%cnew)
       end do
       call aexchange(nc,myid,nprocs,taub)
    endif
# endif

    !J. Ge for tracer advection
    IF(BACKWARD_ADVECTION==.TRUE.)THEN
       IF(BACKWARD_STEP==1)THEN
          do i=1,nst
             SED0(i)%conc = sed(i)%cnew
          end do
       ELSEIF(BACKWARD_STEP==2)THEN
          do i=1,nst
             SED2(i)%conc  = SED0(i)%conc 
             SED0(i)%conc  = sed(i)%cnew 
          end do
       ENDIF
    ENDIF
    !J. Ge for tracer advection
    if(dbg_set(dbg_sbr)) write(ipt,*) "End: finalize_sed_step"

    Return

  End Subroutine finalize_sed_step


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

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: UPDATE_DEPTH"

  END SUBROUTINE UPDATE_DEPTH




# endif


End Module Mod_Sed_CSTMS
