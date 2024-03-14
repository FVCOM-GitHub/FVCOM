MODULE mod_vegetation 

# if defined VEGETATION
  Use Mod_Par
  Use Mod_Prec 
  Use Mod_Types
  Use Mod_wd
  Use all_vars
  Use Mod_Nctools, Only : NCFILE, NCVAR, NCDIM, FILEHEAD, ADD,  &
                          FIND_FILE, FIND_DIM, FIND_VAR,        &
                          NC_INIT, NC_OPEN, NC_LOAD, NC_CONNECT_AVAR, &
                          NC_READ_VAR, NC_DISCONNECT
# if defined (WAVE_CURRENT_INTERACTION)
  USE MOD_WAVE_CURRENT_INTERACTION
# endif
  !
  implicit none
  !
  !  ************************************************************************
  !  **                                                                    **
  !  ** This module declares vegetation model internal parameters.  Some   **
  !  ** of these parameters are usually read from the appropriate input    **
  !  ** script.                                                            **
  !  **                                                                    **
  !  ** The current design allows the user to have a lot of latitude for   **
  !  ** customizing or expanding the vegetation model.                     **
  !  **                                                                    **
  !  ** The vegetatation model is composed of several files:               **
  !  **                                                                    **
  !  **   *  Vegetation modifies rhs3d terms:                              **
  !  **                                                                    **
  !  **        vegetation_drag.F                                           **
  !  **                                                                    **
  !  **   *  Vegetation modifies turbulence terms:                         **
  !  **                                                                    **
  !  **        vegetation_turb_cal.F                                       **
  !  **                                                                    **
  !  **   *  Vegetation modifies streaming terms:                          **
  !  **                                                                    **
  !  **        vegetation_stream.F                                         **
  !  **                                                                    **
  !  **   *  Wave Thrust on Marsh calculation:                             **
  !  **                                                                    **
  !  **        marsh_wave_thrust.F                                         ** 
  !  **                                                                    **
  !  **   *  Vegetation biomass calculation:                               **
  !  **                                                                    **
  !  **        vegetation_biomass.F                                        **
  !  **                                                                    **
  !  **   *  Internal model parameters declarations:                       **
  !  **                                                                    **
  !  **        vegetation_mod.h                                            **
  !  **                                                                    **
  !  **   *  Model parameters standard input script:                       **
  !  **                                                                    **
  !  **        vegetation.in                                               **
  !  **                                                                    **
  !  **   *  Code to read input model parameters:                          **
  !  **                                                                    **
  !  **        vegetation_inp.h                                            **
  !  **                                                                    **
  !  **   *  Code to assign indices to model variables during the          **
  !  **      reading of metadata information from "varinfo.dat":           **
  !  **                                                                    **
  !  **        vegetation_var.h                                            **
  !  **                                                                    **
  !  **   *  Code to define input model parameters in all output           **
  !  **      NetCDF files:                                                 **
  !  **                                                                    **
  !  **        vegetation_def_his.h                                        **
  !  **                                                                    **
  !  **   *  Code to write out input model parameters in all output        **
  !  **      NetCDF files:                                                 **
  !  **                                                                    **
  !  **        vegetation_wrt_his.h                                        **
  !  **                                                                    **
  !  ** Note that all the files are located in ROMS/Nonlinear/Vegetation   **
  !  ** and the *.h files are included within  <...> to allow the user     **
  !  ** customize any of them in the project directory, while  keeping     **
  !  ** the distributed code intact (check build scripts for details).     **
  !  **                                                                    **
  !  ************************************************************************
  !  */


  !svn $Id: vegetation_mod.h 429 2015-06-10 17:30:26Z arango $
  !================================================== Hernan G. Arango ===
  !  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
  !    Licensed under a MIT/X style license                              !
  !    See License_ROMS.txt                                              !
  !================================================= John C. Warner =====!   
  !================================================= Neil K. Ganju ======!   
  !================================================= Alexis Beudin ======!   
  !================================================= Tarandeep S. Kalra =!
  !=======================================================================
  !                                                                      !
  !  Vegetation Model Kernel Variables:                                  !
  !  =================================                                   !
  !  NVEG          Number of vegetation types                            !
  !  NVEGP         Number of vegetation array properties                 !
  !  CD_VEG        Drag coefficient for each veg type                    ! 
  !  E_VEG         Youngs modulus for each veg type                      !
  !  VEG_MASSDEN   Mass density for each veg type                        !
  !  VEGHMIXCOEF   Viscosity coefficient for vegetation boundary         ! 
  !                                                                      ! 
  !  Plant Property indices:                                             !
  !  ======================                                              !
  !  pdens         Density                                               !
  !  phght         Height                                                !
  !  pdiam         Diameter                                              !
  !  pthck         Thickness                                             !
  !  pabbm         Above ground biomass                                  !
  !  pbgbm         Below ground biomass                                  !
  !                                                                      !
  !  Plant Property indices:                                             !
  !  ======================                                              !
  !  idvprp        Indices for storing plant properties                  ! 
  !                                                                      !
  !  Plant Property Output IDs:                                          !
  !  ==========================                                          !
  !  ipdens         Id to output plant density                           !
  !  iphght         Id to output plant height                            !
  !  ipdiam         Id to output plant diameter                          !
  !  ipthck         Id to output plant thickness                         !
  !  ipupbm         Id to output above ground biomass                    !
  !  ipdwbm         Id to output below ground biomass                    !
  !  idWdvg         Id to output wave dissipation from vegetation        !
  !                                                                      !
  !  Wave Thrust on Marsh Output:                                        !
  !  ==========================                                          !
  !  idTims        Initial masking for the marsh                         ! 
  !  idTmsk        Masking for getting thrust due to waves at rho pts.   ! 
  !  idTmax        Maximum thrust due to waves                           !
  !  idTton        Tonelli masking based thrust due to waves             !
  !=======================================================================
  !

  ! Main controls for vegetation
  LOGICAL :: VEG_BIOMASS
  LOGICAL :: VEG_STREAMING
  LOGICAL :: MARSH_WAVE_THRUST
  LOGICAL :: VEG_FLEX
  LOGICAL :: VEG_HMIXING
  LOGICAL :: VEG_TURB
  LOGICAL :: VEG_DRAG
  LOGICAL :: VEG_SWAN_COUPLING

  integer :: NVEG, NVEGP
  integer :: counter
  integer :: phght, pdens, pdiam, pthck
  integer :: ipdens,iphght,ipdiam,ipthck

  !#ifdef VEG_BIOMASS 
  integer :: pabbm, pbgbm
  integer :: ipabbm, ipbgbm
  !#endif 
  !#ifdef VEG_STREAMING 
  integer :: idWdvg
  !#endif 
  integer, allocatable :: idvprp(:)
  !#ifdef MARSH_WAVE_THRUST 
  integer ::  idTims, idTmsk, idTmax, idTton
  !#endif 
  !
  real(sp), allocatable :: E_VEG(:)
  real(sp), allocatable :: CD_VEG(:)
  real(sp), allocatable :: VEG_MASSDENS(:)
  real(sp), allocatable :: VEGHMIXCOEF(:)
  character(len=20), allocatable ::vegname(:)
  !======================================================================!
  !                                                                      !
  !  Vegetation Model Kernel Variables:                                  !
  !  plant         Vegetation variable properties:                       !
  !                   plant(:,:,:,phght) => height                       !
  !                   plant(:,:,:,pdens) => density                      !
  !                   plant(:,:,:,pthck) => thickness                    !
  !                   plant(:,:,:,pdiam) => diameter                     !
  !                   plant(:,:,:,pabbm) => above ground biomass         !
  !                   plant(:,:,:,pbgbm) => below ground biomass         !
  !  ru_veg         Momentum term for x direction(takes account for all  !
  !                 vegetation types)                                    !
  !  rv_veg         Momentum term for x direction(takes account for all  !
  !                 vegetation types)                                    !
  !  ru_veg_loc     Momentum term for x direction(takes account for only !
  !                 local vegetation type)                               !
  !  rv_veg_loc     Momentum term for x direction(takes account for all  !
  !                 local vegetation types)                              !
  !  step2d_uveg    Momentum term for 2d x direction                     !
  !  step2d_vveg    Momentum term for 2d y direction                     !
  !  bend           Bending for each vegetation                          !
  !  Lveg           Effective blade length                               ! 
  !  tke_veg        Turbulent kinetic energy from vegetation             !
  !  gls_veg        Length scale change from vegetation                  !
  !  dissip_veg     Dissipation from the SWAN model due to vegetation    !
  !  BWDXL_veg      Wave streaming effect due to vegetation              !
  !  BWDYL_veg      Wave streaming effect due to vegetation              !
  !  visc2d_r_veg   Effect of viscosity change at vegetation interface   ! 
  !  visc3d_r_veg   Effect of viscosity change at vegetation interface   ! 
  !  marsh_mask     User input of marsh masking at MSL                   ! 
  !  mask_thrust    Tonellis masking for wave thrust on marshes          !
  !  Thrust_max     Maximum thrust from wave to marshes                  !
  !  Thrust_tonelli Reduced thrust from tonelli's masking                !
  !                                                                      !
  !======================================================================!
  TYPE T_VEG
     !
     !  Nonlinear model state.
     !
     real(sp), allocatable :: plant(:,:,:)
     character(len=20),pointer :: vegname(:)
     character(len=20),pointer :: paraname(:)
     character(len=20),pointer :: vegunits(:)

     !  Momentum terms go back to act as sink in rhs
     real(sp), allocatable :: ru_veg(:,:)
     real(sp), allocatable :: rv_veg(:,:)

     !  Momentum terms feed to the turbulence model 
     real(sp), allocatable :: ru_loc_veg(:,:,:)
     real(sp), allocatable :: rv_loc_veg(:,:,:)
     real(sp), allocatable :: step2d_uveg(:)
     real(sp), allocatable :: step2d_vveg(:)
     real(sp), allocatable :: Lveg(:,:)
     !# ifdef VEG_FLEX 
     real(sp), allocatable :: bend(:,:)
     !# endif         
     !# ifdef VEG_TURB
     real(sp), allocatable :: tke_veg(:,:)
     real(sp), allocatable :: gls_veg(:,:)
     !# endif 
     !# ifdef VEG_HMIXING
     real(sp), allocatable :: visc2d_r_veg(:)
     real(sp), allocatable :: visc3d_r_veg(:,:)
     !# endif 
     !# if defined VEG_SWAN_COUPLING && defined VEG_STREAMING
     real(sp), allocatable :: dissip_veg(:)
     real(sp), allocatable :: BWDXL_veg(:,:)
     real(sp), allocatable :: BWDYL_veg(:,:)
     !# endif 
     !# ifdef MARSH_WAVE_THRUST
     real(sp), allocatable :: marsh_mask(:)
     real(sp), allocatable :: mask_thrust(:)
     real(sp), allocatable :: Thrust_max(:)
     real(sp), allocatable :: Thrust_tonelli(:)
     !# endif

  END TYPE T_VEG

  TYPE (T_VEG), target :: VEG
  !

  character(len=500) :: vegfile ! input file for vegetation type and parameters

  real(sp),parameter :: rhow = 1025.0 ! constant sea water density

  !-----------------------------------------------------------------------
  !  Generic Length Scale parameters.
  !-----------------------------------------------------------------------
  !
  !    gls_Gh0
  !    gls_Ghcri
  !    gls_Ghmin
  !    gls_Kmin      Minimum value of specific turbulent kinetic energy.
  !    gls_Pmin      Minimum Value of dissipation.
  !    gls_cmu0      Stability coefficient (non-dimensional).
  !    gls_c1        Shear production coefficient (non-dimensional).
  !    gls_c2        Dissipation coefficient (non-dimensional).
  !    gls_c3m       Buoyancy production coefficient (minus).
  !    gls_c3p       Buoyancy production coefficient (plus).
  !    gls_E2
  !    gls_m         Turbulent kinetic energy exponent (non-dimensional).
  !    gls_n         Turbulent length scale exponent (non-dimensional).
  !    gls_p         Stability exponent (non-dimensional).
  !    gls_sigk      Constant Schmidt number (non-dimensional) for
  !                    turbulent kinetic energy diffusivity.
  !    gls_sigp      Constant Schmidt number (non-dimensional) for
  !                    turbulent generic statistical field, "psi".
  real(sp) :: gls_m
  real(sp) :: gls_n
  real(sp) :: gls_p
  real(sp) :: gls_sigk
  real(sp) :: gls_sigp
  real(sp) :: gls_cmu0
  real(sp) :: gls_cmupr
  real(sp) :: gls_c1
  real(sp) :: gls_c2
  real(sp) :: gls_c3m
  real(sp) :: gls_c3p
  real(sp) :: gls_Kmin
  real(sp) :: gls_Pmin

  !  Momentum terms go back to act as sink 
  real(sp), allocatable :: ru_veg(:,:)
  real(sp), allocatable :: rv_veg(:,:)
  real(sp), allocatable :: step2d_uveg(:)
  real(sp), allocatable :: step2d_vveg(:)
 !  Node-based velocity
  real(sp), allocatable :: uu_veg(:,:)
  real(sp), allocatable :: vv_veg(:,:)


  LOGICAL VEGETATION_MODEL
  CHARACTER(LEN=80) VEGETATION_MODEL_FILE
  CHARACTER(LEN=80) VEGETATION_FIELDS_TYPE
  CHARACTER(LEN=80) VEGETATION_FIELDS_FILE
  INTEGER,PARAMETER :: VegStrLen = 20

  NAMELIST /NML_VEGETATION/               &
         & VEGETATION_MODEL,              &
         & VEGETATION_MODEL_FILE,         &
         & VEGETATION_FIELDS_TYPE,        &
         & VEGETATION_FIELDS_FILE
  
CONTAINS

  SUBROUTINE NAME_LIST_INITIALIZE_VEGETATION
    
    IMPLICIT NONE

    !--Parameters in Namelist NML_VEGETATION
    VEGETATION_MODEL = .FALSE.
    VEGETATION_MODEL_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
    VEGETATION_FIELDS_TYPE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
    VEGETATION_FIELDS_FILE = "DO NOT ADD UNTILL FVCOM IS RUNNING BY ITS SELF FIRST"
 
  END SUBROUTINE NAME_LIST_INITIALIZE_VEGETATION

  SUBROUTINE NAME_LIST_PRINT_VEGETATION

    IMPLICIT NONE

    write(UNIT=IPT,NML=NML_VEGETATION)

    RETURN
  END SUBROUTINE NAME_LIST_PRINT_VEGETATION

  SUBROUTINE NAME_LIST_READ_VEGETATION

    IMPLICIT NONE
    integer :: ios, i
    Character(Len=120):: FNAME
    character(len=160) :: pathnfile
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: Read_Name_List_Vegetation;"

    ios = 0

    FNAME = "./"//trim(casename)//"_run.nml"

    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read Additional Models Settings
    READ(UNIT=NMLUNIT, NML=NML_VEGETATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_VEGETATION)
       Call Fatal_Error("Can Not Read NameList NML_VEGETATION from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_VEGETATION)

    CLOSE(NMLUNIT)

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Ends: Read_Name_List_Vegetation;"
  END SUBROUTINE NAME_LIST_READ_VEGETATION

  SUBROUTINE OPEN_FORCING_VEGETATION

    IMPLICIT NONE

    TYPE(NCFILE), POINTER :: NCF
    INTEGER :: charnum
    LOGICAL :: back
    CHARACTER(LEN=160) :: pathnfile

    back = .true.

        ! LOAD VEGETATION FILEDS
    IF (VEGETATION_FIELDS_TYPE /= UNIFORM) THEN

       ! TEST FILE NAME
       charnum = index (VEGETATION_FIELDS_FILE,".nc",back)
       if (charnum /= len_trim(VEGETATION_FIELDS_FILE)-2)&
            & CALL WARNING("Sediment parameter file does not end in .nc", &
            & trim(VEGETATION_FIELDS_FILE))

       ! INITIALIZE TYPE TO HOLD FILE METADATA
       pathnfile= trim(INPUT_DIR)//trim(VEGETATION_FIELDS_FILE)
       CALL  NC_INIT(NCF,pathnfile)

       ! OPEN THE FILE AND LOAD METADATA
       if(.not. NCF%OPEN) then
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)
          FILEHEAD => ADD(FILEHEAD,NCF)
       end if

    END IF

  END SUBROUTINE OPEN_FORCING_VEGETATION
  ! 
  !==========================================================================
  ! Allocate Data, Initialize Vegetation Fields and Parameters       
  !==========================================================================
  Subroutine Setup_Vegetation
    character(len=80)    :: fname
    logical             :: fexist

    fname = VEGETATION_MODEL_FILE
    !ensure vegetation parameter file exists
    if(fname == "none" .or. fname == "NONE" .or. fname == "None")then
       vegfile = "./"//trim(input_dir)//"/"//trim(casename)//'_vegetation.inp'
    else
       vegfile = "./"//trim(input_dir)//"/"//trim(fname)
    endif

    inquire(file=trim(vegfile),exist=fexist)
    if(.not.fexist)then
       write(*,*)'vegetation parameter file: ',trim(vegfile),' does not exist'
       write(*,*)'stopping'
       call pstop
    end if

    Call read_VegPar

    Call initialize_vegetation

    Call allocate_vegarr

    Call read_VegFields



  End Subroutine Setup_Vegetation


  Subroutine Update_Vegetation
    !
    implicit none 
    !

    CALL E2N3D(U,UU_VEG)
    CALL E2N3D(V,VV_VEG)

    if(VEG_DRAG)then
       Call Calc_vegetation_drag
    end if

    if(VEG_TURB)then
       Call Calc_vegetation_turb
    end if

    if(VEG_HMIXING)then
       Call Calc_vegetation_hmixing
    end if
# if defined (WAVE_CURRENT_INTERACTION)
    if(VEG_SWAN_COUPLING .and. VEG_STREAMING)then
       Call Calc_vegetation_stream
    end if
# endif

  End Subroutine Update_Vegetation



  SUBROUTINE initialize_vegetation
    !

    !
    implicit none 
    !
    !     Setup property indices 
    ! 
    counter = 1
    pdens   = counter 
    counter = counter+1 
    phght   = counter
    counter = counter+1 
    pdiam   = counter
    counter = counter+1 
    pthck   = counter
    if(VEG_BIOMASS)then 
       counter = counter+1 
       pabbm   = counter
       counter = counter+1 
       pbgbm   = counter 
    endif
    NVEGP = counter

    RETURN
  END SUBROUTINE initialize_vegetation

  SUBROUTINE allocate_vegarr
    !
    !=======================================================================
    !                                                                      !
    !  This routine allocates all variables in the module for all nested   !
    !  grids.                                                              !
    !                                                                      !
    !=======================================================================
    !


    implicit none 
    !                       
    !  Imported variable declarations.
    !

    !
    !-----------------------------------------------------------------------
    !  Allocate structure variables.
    !-----------------------------------------------------------------------
    !

    allocate ( VEG % plant(0:mt,NVEG,NVEGP) )    ; VEG % plant = 0.0_sp
    allocate ( VEG % vegname(NVEG))              ; VEG % vegname = ""
    allocate ( VEG % paraname(NVEGP))            ; VEG % paraname = ""
    allocate ( VEG % vegunits(NVEGP))            ; VEG % vegunits = ""

    VEG % vegunits(pdens) = "plant-meter2"
    VEG % vegunits(phght) = "meter"
    VEG % vegunits(pdiam) = "meter"
    VEG % vegunits(pthck) = "meter"
    if(VEG_BIOMASS)then
       VEG % vegunits(pabbm) = "kg/m^2"
       VEG % vegunits(pbgbm) = "kg/m^2"
    endif

    VEG % vegname = vegname

    VEG % paraname(pdens) = "plant_density"
    VEG % paraname(phght) = "plant_height"
    VEG % paraname(pdiam) = "plant_diameter"
    VEG % paraname(pthck) = "plant_thickness"
    if(VEG_BIOMASS)then
       VEG % paraname(pabbm) = "plant_biomass_above_ground"
       VEG % paraname(pbgbm) = "plant_biomass_below_ground"
    endif


    allocate ( VEG % ru_veg(0:mt,kb) )           ; VEG % ru_veg = 0.0_sp
    allocate ( VEG % rv_veg(0:mt,kb) )           ; VEG % rv_veg = 0.0_sp
    allocate ( VEG % ru_loc_veg(0:mt,kb,NVEG) )  ; VEG % ru_loc_veg = 0.0_sp
    allocate ( VEG % rv_loc_veg(0:mt,kb,NVEG) )  ; VEG % rv_loc_veg = 0.0_sp
    allocate ( VEG % step2d_uveg(0:mt) )         ; VEG % step2d_uveg = 0.0_sp
    allocate ( VEG % step2d_vveg(0:mt) )         ; VEG % step2d_vveg = 0.0_sp
    allocate ( VEG % Lveg(0:mt,kb) )             ; VEG % Lveg = 0.0_sp
    if(VEG_FLEX)then
       allocate ( VEG % bend(0:mt,NVEG) )           ; VEG % bend = 0.0_sp
    endif
    if(VEG_HMIXING)then
       allocate ( VEG % visc2d_r_veg(0:mt) )        ; VEG % visc2d_r_veg = 0.0_sp
       allocate ( VEG % visc3d_r_veg(0:mt,kb) )     ; VEG % visc3d_r_veg = 0.0_sp
    endif
    if(VEG_TURB)then
       allocate ( VEG % tke_veg(0:mt,kb) )          ; VEG % tke_veg = 0.0_sp
       allocate ( VEG % gls_veg(0:mt,kb) )          ; VEG % gls_veg = 0.0_sp
    endif
    if(VEG_SWAN_COUPLING .and. VEG_STREAMING)then
       allocate ( VEG % dissip_veg(0:mt) )          ; VEG % dissip_veg = 0.0_sp
       allocate ( VEG % BWDXL_veg(0:mt,kb) )        ; VEG % BWDXL_veg = 0.0_sp
       allocate ( VEG % BWDYL_veg(0:mt,kb) )        ; VEG % BWDYL_veg = 0.0_sp
    endif
    if(MARSH_WAVE_THRUST)then
       allocate ( VEG % marsh_mask(0:mt) )          ; VEG % marsh_mask = 0.0_sp
       allocate ( VEG % mask_thrust(0:mt) )         ; VEG % mask_thrust = 0.0_sp
       allocate ( VEG % Thrust_max(0:mt) )          ; VEG % Thrust_max = 0.0_sp
       allocate ( VEG % Thrust_tonelli(0:mt) )      ; VEG % Thrust_tonelli = 0.0_sp
    endif

    !
    !-----------------------------------------------------------------------
    !  Allocate various cell-based variables for momentum calculation
    !-----------------------------------------------------------------------
    !
    allocate ( ru_veg(0:nt,kb) )           ; ru_veg = 0.0_sp
    allocate ( rv_veg(0:nt,kb) )           ; rv_veg = 0.0_sp
    allocate ( step2d_uveg(0:nt) )         ; step2d_uveg = 0.0_sp
    allocate ( step2d_vveg(0:nt) )         ; step2d_vveg = 0.0_sp
    !
    !-----------------------------------------------------------------------
    !  Allocate various node-based variables for momentum calculation
    !-----------------------------------------------------------------------
    !
    allocate ( uu_veg(0:mt,kb) )              ; uu_veg = 0.0_sp
    allocate ( vv_veg(0:mt,kb) )              ; vv_veg = 0.0_sp
    !
    !-----------------------------------------------------------------------
    !  Allocate various input variables for vegetation module.
    !-----------------------------------------------------------------------
    !

    RETURN
  END SUBROUTINE allocate_vegarr


  SUBROUTINE read_VegPar
    !
    !=======================================================================
    !                                                                      !
    !  This routine reads in vegetation model parameters.                  !
    !  Equivalent of read_phypar.F so it gets read in that                 !
    !  This routine also outputs vegetation model parameters.              !
    !=======================================================================
    !
    Use Input_Util
    !
    implicit none
    !
    !  Imported variable declarations
    !
    !
    !  Local variable declarations.
    !
    integer linenum,i,k1,iscan
    real(sp)           :: ftemp
    character(len=120) :: stemp
    real(sp),allocatable:: vegtmp(:)

    integer :: Npts, Nval
    integer :: iveg, ng, status
    !
    real(sp), dimension(200) :: Rval
    real(sp), allocatable :: Rveg(:,:)
    !
    character (len=40 ) :: KeyWord
    character (len=256) :: line
    character (len=256), dimension(200) :: Cval
    character(len=80),allocatable:: strtmp(:)
    !
    !-----------------------------------------------------------------------
    !  Read input parameters from vegetation.in
    !-----------------------------------------------------------------------
    !

    !read in   NVEG          Number of vegetation types     
    Call Get_Val(NVEG,vegfile,'NVEG',line=linenum,echo=.true.) 

    IF (.not.allocated(CD_VEG))       allocate(CD_VEG(NVEG))
    IF (.not.allocated(E_VEG))        allocate(E_VEG(NVEG))
    IF (.not.allocated(VEG_MASSDENS)) allocate(VEG_MASSDENS(NVEG))
    IF (.not.allocated(VEGHMIXCOEF))  allocate(VEGHMIXCOEF(NVEG))
    IF (.not.allocated(vegname))  allocate(vegname(NVEG))

    allocate(vegtmp(nveg))

    !read in   CD_VEG        Drag coefficient for each veg type 
    if(nveg>1)then
       Call Get_Val_Array(vegtmp,vegfile,'CD_VEG',nveg,echo=.true.) 
       CD_VEG(1:NVEG)=vegtmp
    else
       Call Get_Val(CD_VEG(1),vegfile,'CD_VEG',line=linenum,echo=.true.)
    end if

    !read in vegname         Name for eache veg type
    if(nveg>1)then
       allocate(strtmp(nveg))
       Call Get_Val_Array(strtmp,vegfile,'VEG_NAME',nveg,echo=.true.)
       vegname(1:nveg)=strtmp
    else
       Call Get_Val(vegname(1),vegfile,'VEG_NAME',line=linenum,echo=.true.)
    end if

    !read in   E_VEG         Youngs modulus for each veg type 
    if(nveg>1)then
       Call Get_Val_Array(vegtmp,vegfile,'E_VEG',nveg,echo=.true.) 
       E_VEG(1:NVEG)=vegtmp
    else
       Call Get_Val(E_VEG(1),vegfile,'E_VEG',line=linenum,echo=.true.)
    end if

    !read in   VEG_MASSDEN   Mass density for each veg type                        !
    if(nveg>1)then
       Call Get_Val_Array(vegtmp,vegfile,'VEG_MASSDENS',nveg,echo=.true.) 
       VEG_MASSDENS(1:NVEG)=vegtmp
    else
       Call Get_Val(VEG_MASSDENS(1),vegfile,'VEG_MASSDENS',line=linenum,echo=.true.)
    end if

    !read in   VEGHMIXCOEF   Viscosity coefficient for vegetation boundary         ! 
    if(nveg>1)then
       Call Get_Val_Array(vegtmp,vegfile,'VEGHMIXCOEF',nveg,echo=.true.) 
       VEGHMIXCOEF(1:NVEG)=vegtmp
    else
       Call Get_Val(VEGHMIXCOEF(1),vegfile,'VEGHMIXCOEF',line=linenum,echo=.true.)
    end if

    ! read in VEG_FLEX switch 
    Call Get_Val(VEG_FLEX,vegfile,'VEG_FLEX',line=linenum,echo=.true.)

    ! read in VEG_TURB switch 
    Call Get_Val(VEG_TURB,vegfile,'VEG_TURB',line=linenum,echo=.true.)

    ! read in VEG_SWAN_COUPLING switch 
    Call Get_Val(VEG_SWAN_COUPLING,vegfile,'VEG_SWAN_COUPLING',line=linenum,echo=.true.)

    ! read in VEG_HMIXING  switch 
    Call Get_Val(VEG_HMIXING,vegfile,'VEG_HMIXING',line=linenum,echo=.true.)

    ! read in VEG_DRAG switch 
    Call Get_Val(VEG_DRAG,vegfile,'VEG_DRAG',line=linenum,echo=.true.)

    ! read in VEG_BIOMASS switch 
    Call Get_Val(VEG_BIOMASS,vegfile,'VEG_BIOMASS',line=linenum,echo=.true.)

    ! read in VEG_STREAMING switch 
    Call Get_Val(VEG_STREAMING,vegfile,'VEG_STREAMING',line=linenum,echo=.true.)

    ! read in MARSH_WAVE_THRUST switch 
    Call Get_Val(MARSH_WAVE_THRUST,vegfile,'MARSH_WAVE_THRUST',line=linenum,echo=.true.)

    !read in gls_p
    Call Get_Val(gls_p,vegfile,'GLS_P',line=linenum,echo=.true.)

    !read in gls_m
    Call Get_Val(gls_m,vegfile,'GLS_M',line=linenum,echo=.true.)

    !read in gls_n
    Call Get_Val(gls_n,vegfile,'GLS_N',line=linenum,echo=.true.)

    !read in gls_cmu0
    Call Get_Val(gls_cmu0,vegfile,'GLS_CMU0',line=linenum,echo=.true.)

    !read in gls_Kmin
    Call Get_Val(gls_Kmin,vegfile,'GLS_Kmin',line=linenum,echo=.true.)

    !read in gls_Pmin
    Call Get_Val(gls_Pmin,vegfile,'GLS_Pmin',line=linenum,echo=.true.)

    !read in gls_c1
    Call Get_Val(gls_c1,vegfile,'GLS_C1',line=linenum,echo=.true.)

    !read in 
    Call Get_Val(gls_c2,vegfile,'GLS_C2',line=linenum,echo=.true.)

    RETURN
  END SUBROUTINE read_VegPar

  Subroutine read_VegFields
    USE MOD_UTILS
    USE MOD_NCTOOLS
    USE MOD_INPUT
    USE ALL_VARS
    Use Input_Util

    implicit none

    integer :: iveg
    integer linenum,i,k1,iscan
    real(sp)           :: ftemp
    real(sp),allocatable:: vegtmp(:)

    character (len=256), dimension(200) :: Cval


    if(VEGETATION_FIELDS_TYPE/=UNIFORM)then
       call LOAD_VEGETATION_FIELDS(veg%plant,NVEG,NVEGP)
    else
       !  plant         Vegetation variable properties:                       !
       !                   plant(:,:,phght) => height                       !
       !                   plant(:,:,pdens) => density                      !
       !                   plant(:,:,pthck) => thickness                    !
       !                   plant(:,:,pdiam) => diameter                     !
       !                   plant(:,:,pabbm) => above ground biomass         !
       !                   plant(:,:,pbgbm) => below ground biomass         !
       !    counter = 1
       !    pdens   = counter 
       !    counter = counter+1 
       !    phght   = counter
       !    counter = counter+1 
       !    pdiam   = counter
       !    counter = counter+1 
       !    pthck   = counter
       !    if(VEG_BIOMASS)then 
       !       counter = counter+1 
       !       pabbm   = counter
       !       counter = counter+1 
       !       pbgbm   = counter 
       !    endif
       !    NVEGP = counter
       allocate(vegtmp(NVEGP))
       do iveg=1,nveg
          Call Get_Val_Array(vegtmp,vegfile,trim(vegname(iveg)),NVEGP,echo=.true.) 
          veg%plant(:,iveg,pdens)=vegtmp(1)
          veg%plant(:,iveg,phght)=vegtmp(2)
          veg%plant(:,iveg,pdiam)=vegtmp(3)
          veg%plant(:,iveg,pthck)=vegtmp(4)
          if(VEG_BIOMASS)then 
             veg%plant(:,iveg,pabbm)=vegtmp(5)
             veg%plant(:,iveg,pbgbm)=vegtmp(6)
          end if
       end do
    end if

  End Subroutine read_VegFields

  !=======================================================================
  !                                                                      !
  !  This routine computes the vegetation (posture-dependent) drag       !
  !  for rhs3d.F                                                         !
  !                                                                      !  
  !  References:                                                         !
  !                                                                      !
  !  Luhar M., and H. M. Nepf (2011), Flow-induced reconfiguration of    !
  !  buoyant and flexible aquatic vegetation, Limnology and Oceanography,!
  !   56(6): 2003-2017.                                                  !
  !                                                                      !
  !=======================================================================

  SUBROUTINE Calc_vegetation_drag 

    !
    !  Imported variable declarations.
    !
    !

    !
    !  Local variable declarations.
    !
    integer :: i, j, k, iveg
    ! 
    real(sp), parameter :: one_third  = 1.0_sp/3.0_sp
    real(sp), parameter :: one_twelfth  = 1.0_sp/12.0_sp
    real(dp), parameter :: rad2deg = 180.0_dp / pi
    real(sp), parameter :: Inival  = 0.0_sp
    real(sp) :: cff, cff1, cff2, cff3, cff4, Hz_inverse
    real(sp) :: sma, buoy, Umag, Ca, cflex
    real(sp) :: Lveg_loc, plant_height_eff
    real(sp), dimension(0:mt,kb) :: dab
    real(sp), dimension(0:mt) :: wrk

    !
    !-----------------------------------------------------------------------
    !  Resistance imposed on the flow by vegetation.
    !-----------------------------------------------------------------------
    !
    dab=Inival
    VEG%ru_veg=Inival
    VEG%rv_veg=Inival
    VEG%Lveg=Inival
    !
# if defined (WET_DRY)
    !
    !  Set limiting factor for drag force. The drag force is adjusted
    !  to not change the direction of momentum.  It only should slow down
    !  to zero. The value of 0.75 is arbitrary limitation assigment
    !  (same as for bottom stress).
    !
    if(WETTING_DRYING_ON)then
       cff=0.75_sp/dti
    end if
# endif

    VEG_LOOP: DO iveg=1,NVEG
       K_LOOP: DO k=kbm1,1,-1
          DO i=1,m
             ! cycle the rest calculation if no vegetation covered
             if(VEG%plant(i,iveg,phght)==0)cycle

             if(VEG_FLEX)then
                !
                ! Flexible vegetation
                !
                ! Second moment of area
                ! 
                sma=(VEG%plant(i,iveg,pdiam)*                               &
                     &             VEG%plant(i,iveg,pthck)**3.0_sp)*(one_twelfth)
                !
                ! Buoyancy parameter
                !
                buoy=(rhow-veg_massdens(iveg))*grav*VEG%plant(i,iveg,pdiam)*&
                     &              VEG%plant(i,iveg,pthck)*                              &
                     &              VEG%plant(i,iveg,phght)**3.0_sp/(E_veg(iveg)*sma)
                !
                ! Current speed at rho points
                !
                cff2=uu_veg(i,k)
                cff3=vv_veg(i,k)
                Umag=SQRT(cff2*cff2+cff3*cff3)
                !
                ! Cauchy number
                !
                Ca=0.5_sp*rhow*Cd_veg(iveg)*VEG%plant(i,iveg,pdiam)*     &
                     &                   Umag**2.0_sp*VEG%plant(i,iveg,phght)**3.0_sp/    &
                     &                       (E_veg(iveg)*sma)
                !
                cflex=1.0_sp-((1.0_sp-0.9_sp*Ca**(-one_third))/           &
                     &             (1.0_sp+(Ca**(-1.5_sp)*(8.0_sp+buoy**(1.5_sp)))))
                !
                ! To avoid NaN value when Ca is zero
                !
                cflex=MIN(cflex,1.0_sp)
                !
                ! Effective blade length
                !
                plant_height_eff=cflex*VEG%plant(i,iveg,phght)
                !
                ! Blade bending angle
                !
                VEG%bend(i,iveg)=ACOS(cflex**one_third)*rad2deg
             else
                !
                ! For stiff vegetation
                !
                plant_height_eff=VEG%plant(i,iveg,phght)
             endif
             !
             ! Select the grid cell (full or part) within the canopy layer
             !
             dab(i,k)=dab(i,k+1)+dt(i)*dz(i,k)
             Hz_inverse=1.0_sp/(dt(i)*dz(i,k))
             cff1=MIN((dab(i,k)-plant_height_eff)*Hz_inverse,1.0_sp)
             Lveg_loc=MIN(1.0_sp-cff1,1.0_sp)
             !
             ! Prepare drag term (at rho points)
             !       
             wrk(i)=0.5_sp*cd_veg(iveg)*VEG%plant(i,iveg,pdiam)*    &
                  &                 VEG%plant(i,iveg,pdens)*(dt(i)*dz(i,k))*Lveg_loc
             !
             ! Store Lveg_loc for all vegetation types
             !
             VEG%Lveg(i,k)=Lveg_loc+VEG%Lveg(i,k)
          END DO
          !
          ! Compute friction force (at cell faces)
          !
          DO i=1,m
             ! cycle the rest calculation if no vegetation covered
             if(VEG%plant(i,iveg,phght)==0)cycle

             cff2=SQRT(uu_veg(i,k)*uu_veg(i,k)+vv_veg(i,k)*vv_veg(i,k))
             cff3=uu_veg(i,k)*cff2
             VEG%ru_loc_veg(i,k,iveg)=wrk(i)*cff3
             !
             !  Add the ru_iveg from this veg type to another veg type
             !  which can be there at the same point (i,k)
             !  Alexis's comment: not confident in what is happening when
             !                     multiple vegetation types are concomitant
             !
             VEG%ru_veg(i,k)=VEG%ru_loc_veg(i,k,iveg)+VEG%ru_veg(i,k)

# ifdef WET_DRY
             if(WETTING_DRYING_ON)then
                cff4=cff*dt(i)*dz(i,k)
                VEG%ru_veg(i,k)=SIGN(1.0_sp, VEG%ru_veg(i,k))*                &
                     &               MIN(ABS(VEG%ru_veg(i,k)),                            &
                     &                   ABS(uu_veg(i,k))*cff4)
             end if
# endif


          END DO

          DO i=1,m
             ! cycle the rest calculation if no vegetation covered
             if(VEG%plant(i,iveg,phght)==0)cycle

             cff2=SQRT(uu_veg(i,k)*uu_veg(i,k)+vv_veg(i,k)*vv_veg(i,k))
             cff3=vv_veg(i,k)*cff2
             VEG%rv_loc_veg(i,k,iveg)=wrk(i)*cff3
             !
             !   Add the rv_iveg from this veg type to another veg type
             !   which can be there at the same point (i,k)
             !
             VEG%rv_veg(i,k)=VEG%rv_loc_veg(i,k,iveg)+VEG%rv_veg(i,k)

# ifdef WET_DRY
             if(WETTING_DRYING_ON)then
                cff4=cff*dt(i)*dz(i,k)
                VEG%rv_veg(i,k)=SIGN(1.0_sp, VEG%rv_veg(i,k))*                &
                     &               MIN(ABS(VEG%rv_veg(i,k)),                            &
                     &                   ABS(vv_veg(i,k))*cff4)
             endif
# endif


          END DO
       END DO K_LOOP
    END DO VEG_LOOP
    !
    !-----------------------------------------------------------------------
    !  Add in resistance imposed on the flow by the vegetation (3D->2D).
    !  Changes feedback in Nonlinear/step2d_LF_AM3.F
    !-----------------------------------------------------------------------
    !
    DO i=1,m
       cff=dt(i)*dz(i,kbm1)
       cff2=cff*VEG%ru_veg(i,kbm1)
       DO k=kbm1-1,1,-1
          cff=dt(i)*dz(i,k)
          cff2=cff2+cff*VEG%ru_veg(i,k)
       END DO
       VEG%step2d_uveg(i)=cff2
    END DO
    !
    DO i=1,m
       cff=dt(i)*dz(i,kbm1)
       cff2=cff*VEG%rv_veg(i,kbm1)
       DO k=kbm1-1,1,-1
          cff=dt(i)*dz(i,k)
          cff2=cff2+cff*VEG%rv_veg(i,k)
       END DO
       VEG%step2d_vveg(i)=cff2
    END DO

    CALL N2E2D(VEG%step2d_uveg,step2d_uveg)
    CALL N2E2D(VEG%step2d_vveg,step2d_vveg)
    CALL N2E3D(VEG%ru_veg,ru_veg)
    CALL N2E3D(VEG%rv_veg,rv_veg)
    !
# if defined (MULTIPROCESSOR)
    if(par)then
       call aexchange(nc,myid,nprocs,VEG%Lveg)
       call aexchange(nc,myid,nprocs,VEG%step2d_uveg,VEG%step2d_vveg )
       call aexchange(nc,myid,nprocs,VEG%ru_veg,VEG%rv_veg)
       call aexchange(ec,myid,nprocs,step2d_uveg,step2d_vveg)
       call aexchange(ec,myid,nprocs,ru_veg,rv_veg)
    endif
# endif
    RETURN
  END SUBROUTINE Calc_vegetation_drag


  !=======================================================================
  !                                                                      !
  !  This routine computes the turbulent kinetic energy and length scale !
  !  modifications due to vegetation for gls_corstep.F                   !
  !                                                                      !
  !  References:                                                         !
  !                                                                      !
  !   Uittenbogaard R. (2003): Modelling turbulence in vegetated aquatic !
  !   flows. International workshop on RIParian FORest vegetated         !
  !   channels: hydraulic, morphological and ecological aspects,         !
  !   20-22 February 2003, Trento, Italy.                                !
  !                                                                      !
  !   Warner J.C., C.R. Sherwood, H.G. Arango, and R.P. Signell (2005):  !
  !   Performance of four turbulence closure models implemented using a  !
  !   generic length scale method, Ocean Modelling 8: 81-113.            !
  !                                                                      !
  !=======================================================================
  !***********************************************************************
  SUBROUTINE Calc_vegetation_turb
    !***********************************************************************
    !

    !
    !
    !  Imported variable declarations.
    !

    !
    !  Local variable declarations.
    !
    integer :: i, j, k, iveg
    !
    real(sp), parameter :: one_half=1.0_sp/2.0_sp 
    real(sp), parameter :: one_third=1.0_sp/3.0_sp 
    real(sp), parameter :: Inival=0.0_sp
    real(sp), parameter :: cl_veg=1.0_sp, ck=0.09_sp
    real(sp), parameter :: max_L=10.0e10_sp 
    real(sp), parameter :: eps=1.0e-12_sp 
    real(sp) :: wrku1, wrku2, wrku3, wrku4, wrku
    real(sp) :: wrkv1, wrkv2, wrkv3, wrkv4, wrkv
    real(sp) :: wrk, cff1, cff2, cff3, dissip, inverse_dissip
    real(sp) :: solid, L2, eqvegT
    real(sp) :: taufree, tauveg, taueff
    real(sp) :: tke_loc_veg, gls_loc_veg
    real(sp), dimension(0:mt) :: vegu
    real(sp), dimension(0:mt) :: vegv 

    real(sp), dimension(0:mt,kb) :: gls
# if !defined (GOTM)
    real(sp), dimension(0:mt,kb) :: tke
# endif

    VEG % tke_veg=Inival
    VEG % gls_veg=Inival

    cff1=3.0_sp+gls_p/gls_n
    cff2=1.5_sp+gls_m/gls_n
    cff3=-1.0_sp/gls_n

    VEG_LOOP: DO iveg=1,NVEG
       DO k=KBM1,1,-1
          DO i=1,m
             ! cycle the rest calculation if no vegetation covered
             if(VEG%plant(i,iveg,phght)==0)cycle

             !
             !-----------------------------------------------------------------------
             ! Additional turbulence generated by the vegetation = 
             ! work spent by the fluid against the plants (in m3/s3)
             !-----------------------------------------------------------------------
             !
             wrku=VEG%ru_loc_veg(i,k,iveg)*uu_veg(i,k)
             wrkv=VEG%rv_loc_veg(i,k,iveg)*vv_veg(i,k)
             tke_loc_veg=MAX(sqrt(wrku*wrku+wrkv*wrkv),eps)
             !
             !-----------------------------------------------------------------------
             ! Dissipation due to vegetation
             !-----------------------------------------------------------------------
             ! Dissipation in GLS (Eq. 12 in Warner et al., 2005)
             !
# if defined (GOTM)
             wrk=MAX(tke(i,k),gls_Kmin)
             !  gls= Turbulent energy squared times turbulent length scale       !

             gls(i,k)=tke(i,k)**2*l(i,k)
             dissip=(gls_cmu0**cff1)*(wrk**cff2)*(gls(i,k)**cff3)
# else
             tke(i,k)=q2(i,k)
             gls(i,k)=Q2L(i,k)
             wrk=MAX(tke(i,k),gls_Kmin)
             dissip=(gls_cmu0**cff1)*(wrk**cff2)*(gls(i,k)**cff3)
# endif

             inverse_dissip=1.0_sp/MAX(dissip,eps)
             !
             ! Dissipation time-scale for free turbulence
             !
             taufree=wrk*inverse_dissip
             !
             !# ifdef VEG_FLEX
             !
             ! Equivalent thickness: horizontal projection of the bending plant
             ! 
             !              eqvegT=VEG%plant(i,iveg,pthck)+sin(VEG%bend(i,iveg))*         &
             !     &                                       VEG%plant(i,iveg,phght)
             !# else
             eqvegT=VEG%plant(i,iveg,pthck)
             !# endif                 
             !
             !
             ! Solidity:cross-sectional area of a plant the number of plants per m2
             !
             !
             solid=VEG%plant(i,iveg,pdiam)*eqvegT*VEG%plant(i,iveg,pdens)
             !
             ! Eddies typical size constrained by distance in between the plants
             !
             L2=cl_veg*((1.0_sp-MIN(solid,1.0_sp))/VEG%plant(i,iveg,pdens))**one_half
             L2=MIN(L2,max_L)
             !
             ! Dissipation time-scale of eddies in between the plants
             !
             tauveg=(L2**2.0_sp/(ck**2.0_sp*tke_loc_veg))**one_third
             !
             ! Effective dissipation time-scale
             ! 
             taueff=MIN(taufree,tauveg)
             gls_loc_veg=gls_c2*tke_loc_veg/taueff
             if(i==500)print*,tauveg,L2,ck,tke_loc_veg,one_third
             !
             !-----------------------------------------------------------------------
             ! Add the tke and gls changes from all vegetation types
             !-----------------------------------------------------------------------
             ! 
             VEG % tke_veg(i,k)=tke_loc_veg + VEG % tke_veg(i,k)
             VEG % gls_veg(i,k)=gls_loc_veg + VEG % gls_veg(i,k)
          END DO
       END DO
    END DO VEG_LOOP
    !
# if defined (MULTIPROCESSOR)
    if(par)then
       call aexchange(nc,myid,nprocs,VEG % tke_veg,VEG % gls_veg)
    endif
# endif
    RETURN
  END SUBROUTINE Calc_vegetation_turb


# if defined (WAVE_CURRENT_INTERACTION)
  !***********************************************************************
  !
  ! Since FVCOM does not contain VORTEX force, this subroutine hence
  !  is not included in the iteration.
  !
  SUBROUTINE Calc_vegetation_stream
    !***********************************************************************
    !

    !
    !  Imported variable declarations.
    !

    !  Local variable declarations.
    !
    integer :: i, j, k, iveg
    real(sp) :: cff1, cff2
    real(sp) :: EWD_veg
    real(sp), parameter :: Lwave_min = 1.0_sp

    real(sp) :: Dstp
    real(sp) :: waven, wavenx, waveny
    real(sp) :: sigma, osigma
    !----------------------------------------------------------------------
    !----------Executing the code------------------------------------------
    !----------------------------------------------------------------------
    !     
    DO k=KBM1,1,-1
       DO i=1,m
          Dstp= dt(i)*(zz(i,kbm1)-zz(i,1)) !z_w(i,N)-z_w(i,0)
          !----------------------------------------------------------------------
          !  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
          !----------------------------------------------------------------------
          waven=2.0_sp*pi/MAX(wlen(i),Lwave_min)
          cff1=1.5_sp*pi-Dwave(i)
          wavenx=waven*COS(cff1)
          waveny=waven*SIN(cff1)
          sigma=MIN(SQRT(grav*waven*TANH(waven*Dstp)),2.0_sp)
          osigma=1.0_sp/sigma
          ! 
          !----------------------------------------------------------------------
          !   Note: Alexis - check if we need a local VEG%dissip_veg here 
          !   Also Lveg is for 1 veg type only 
          !----------------------------------------------------------------------
          !
          EWD_veg=VEG % dissip_veg(i)
          cff2=EWD_veg*osigma*VEG % Lveg(i,k)
          VEG % BWDXL_veg(i,k)=cff2*wavenx
          VEG % BWDYL_veg(i,k)=cff2*waveny
          !           
       END DO
    END DO
    !
  END SUBROUTINE Calc_vegetation_stream

  !=======================================================================
  !                                                                      ! 
  !  This routine computes the wave thrust on marshes. Marsh thrust      !
  !  values are computed with correction from the wave angle. For each   !
  !  cell if one side is sheltered from other cells, that side is not    !
  !  exposed to waves. Each cell has four cell normals directed towards  !
  !  the center of the cell. The angle of the normals is with respect to !
  !  the North and clockwise direction. For a submerged marsh,           !
  !  "Tonelli mask" is used to reduce the value of the wave thrust.      !
  !                                                                      !
  !  References:                                                         !   
  !                                                                      !
  !=======================================================================
  !                                                                      !
  !  Tonelli, M., Fagherazzi, Sergio., and Petti., M., 2010: Modeling    !
  !  wave impact on salt marsh boundaries, Journal of Geophysical        !
  !  Research, 115, 0148-0227.                                           !   
  !                                                                      !
  !  Dean, R.G. and Dalrymple, R.A., 1991: Water Wave Mechanics for      !
  !  Engineers and Scientists, World Scientific Publications             !
  !                                                                      !
  !=======================================================================

  SUBROUTINE Calc_marsh_wave_thrust
    !***********************************************************************
    !

    !  Local variable declarations.
    !
    integer :: i,j,cnt,jj,sum,cc
    integer, parameter :: zero = 0
    real(dp), parameter :: rad2deg = 180.0_dp / pi
    integer :: bound_wd1,bound_wd2,cffi
    integer :: sum_l,sum_r,sum_u,sum_d
    integer :: Isolated_node

    real(sp) :: Kw,Integral_Kp,Fw1,Fw2,Fw,cff
    real(sp) :: depth_all,mask_local_tonelli

    real(sp) :: angler_deg,eft_angle
    real(sp) :: dir_l,dir_r,dir_u,dir_d
    real(sp) :: exp
    real(sp) :: thrust_w
    !
    !  Mean density (Kg/m3) used when the Boussinesq approximation is
    !  inferred.
    !
    real(dp) :: rho0 = 1025.0_dp

    integer, dimension(0:mt) :: mask_wd
    integer, dimension(0:mt) :: bound_wdv
    integer, dimension(0:mt) :: bound_wdh
    integer, dimension(0:mt) :: bound_wd
    integer, dimension(0:mt) :: bound_wd_nan

    integer, dimension(0:mt) :: bound_wd_l
    integer, dimension(0:mt) :: bound_wd_r
    integer, dimension(0:mt) :: bound_wd_u
    integer, dimension(0:mt) :: bound_wd_d

    !----------------------------------------------------------------------
    !----------Executing the code---------------------------------------
    !----------------------------------------------------------------------
    DO i=1,m
       mask_wd(i)     = VEG % marsh_mask(i) 
       bound_wdv(i)   = zero          
       bound_wdh(i)   = zero
       bound_wd(i)    = zero   
       bound_wd_nan(i)= zero  
    END DO
    !----------------------------------------------------------------------
    !             FIND THE LAST WET AND DRY BOUNDARY POINT
    !----------------------------------------------------------------------
    !----Point=1 is last wet point; point=-1 are last dry point------------
    !----------------------------------------------------------------------
    ! 

    !----------------------------------------------------------------------
    !----Do the sweep around a node first ------------------------
    !----------------------------------------------------------------------
    do i=1,m
       cnt = 0.
       do jj=1,NTSN(i)
          cc = NBSN(i,jj)
          cnt = cnt + mask_wd(cc)
       end do
       ! if cnt = 0, means fully dry, isolated
       ! if cnt = ntve(i), means fully wet
       if(cnt == NTSN(i))cnt = 0
       bound_wd(i)=cnt
       !----------------------------------------------------------------------
       !-------These are the points used to compute thrust -------------------
       !----------------------------------------------------------------------
       IF(bound_wd(i).ne.0) THEN 
          bound_wd(i)=1 
       ENDIF
       bound_wd_nan(i) = bound_wd(i) 

       !----------------------------------------------------------------------
       !------"10" is corresponding to a flag to not include these pts.-------
       !------ it can be changed to any number greater than 4 ----------------
       !----------------------------------------------------------------------
       IF(bound_wd(i).eq.0) THEN     
          bound_wd_nan(i)=10          
       ENDIF
    end do


    !----------------------------------------------------------------------
    !                              FIND NEIGHBOURS
    !----------------------------------------------------------------------

    DO i=1,m
       cnt = 0
       Isolated_node = 0
       do jj=1,ntsn(i)
          cc = nbsn(i,jj)
          cnt = cnt + mask_wd(cc)
       end do
       sum = cnt + mask_wd(i)
       !----------------------------------------------------------------------
       !---- Have different conditions for isolated cell masking 
       !----------------------------------------------------------------------
       !----------------------------------------------------------------------
       !          if it is not equal to ntsn(i)+1, then it is an isolated node
       !----------------------------------------------------------------------
       if (sum.ne.ntsn(i)+1)then
          isolated_node = 1
          sum = 10
       else
          !----------------------------------------------------------------------
          ! if the cell has a wet neighbor make sum = 1 
          !----------------------------------------------------------------------
          sum = 1
       end if

       Isolated_node = Isolated_node + mask_wd(i)
       !----------------------------------------------------------------------
       !------if isolated cell=2 this means node is isolated + wet 
       !----------------------------------------------------------------------
       IF (Isolated_node.ne.2) THEN 
          Isolated_node = 0
       ELSE 
          Isolated_node = 1 
       ENDIF

       IF (Isolated_node.eq.1) THEN 
          bound_wd(i) = 1 
       ENDIF

       !*********** This subroutine was not finalized from here.******
       !***
       IF (Isolated_node.eq.1) THEN 
          bound_wd(i) = 1 
          sum  = 1
       ENDIF
       !---------------------------------------------------------------------
       ! angler is in radians between xi and actual East ! convert to degrees
       !---------------------------------------------------------------------
       ! no angler in FVCOM
       !angler_deg = angler(i)*180.0_sp/pi

       !---------------------------------------------------------------------
       !                     CALCULATE THE MARSH THRUST
       !---------------------------------------------------------------------
       kw = 2.0_sp*pi/Wlen(i)              
       Integral_kp = sinh(kw*h(i))/(kw*cosh(h(i)*kw))   
       Fw1 = rho0*grav*HSC1(i)*Integral_kp*0.001_sp   
       Fw2 = (rho0*grav*HSC1(i))*HSC1(i)*0.5_sp*0.001_sp  
       !---------------------------------------------------------------------
       !         Total wave thrust at mean sea level 
       !---------------------------------------------------------------------
       Fw  = Fw1 + Fw2 

       IF (sum.ne.10) THEN            
          exp=0.0_sp !angler_deg*sum*bound_wd(i)
          eft_angle = exp - Dwave(i)*rad2deg
          thrust_w=abs(Fw*cos(eft_angle*deg2rad))
          !---------------------------------------------------------------------
          !       If 90>angle>270 then waves arrive from opposite direction with 
          !       respect to the cell side so impact other side of cell 
          !---------------------------------------------------------------------
          !IF (abs(eft_angle).ge.90.0_sp.and.                          &
          !     &          abs(eft_angle).le.270.0_sp) THEN
          !   thrust_w = 0.0_sp 
          !ENDIF
       ELSE  
          thrust_w = 0.0_sp
       ENDIF

       !---------------------------------------------------------------------
       !   if marsh is submerged in water depending on depth, wave thrust to 
       !   be reduced 
       !---------------------------------------------------------------------
       depth_all=d(i)
       IF (depth_all.lt.0.2_sp) THEN
          cff=1.0_sp-0.45_sp*depth_all*5.0_sp
       ELSEIF (0.2_sp.lt.depth_all.and.depth_all.lt.0.4_sp) THEN
          cff=0.55_sp*(1.0_sp-2.5_sp*(depth_all-0.2_sp))
       ELSE 
          cff=0.275_sp 
       ENDIF
       VEG % mask_thrust(i)=cff
       !
       IF (bound_wd_nan(i).ne.10) THEN   
          VEG % Thrust_tonelli(i) = bound_wd_nan(i)*thrust_w*           &
               &                            VEG % mask_thrust(i)  

          VEG % Thrust_max(i) = bound_wd_nan(i)*thrust_w 
       ELSE
          VEG % Thrust_tonelli(i) = 0.0_sp  
          VEG % Thrust_max(i)     = 0.0_sp 
       ENDIF


    END DO


# if defined (MULTIPROCESSOR)
    if(par)then
       call aexchange(nc,myid,nprocs,VEG % mask_thrust,VEG % Thrust_max,VEG % Thrust_tonelli)
    endif
# endif


  END SUBROUTINE Calc_marsh_wave_thrust
# endif
  !                                                                      !
  !  Calculate viscosity change at vegetation iterface and add in        !
  !  hmixing.F                                                           ! 
  !                                                                      !
  !  References:                                                         !   
  !                                                                      !
  !=======================================================================
  !                                                                      !
  !                                                                      !
  !=======================================================================

  !***********************************************************************
  SUBROUTINE Calc_vegetation_hmixing
    !***********************************************************************
    !


    !  Local variable declarations.
    !
    integer :: min_loc(1),max_loc(1)
    integer :: i, j, k, iveg
    real(sp) :: visc3d_r_veg_loc, cff 
    real(sp) :: maxgradvegdens, gradu, gradv, wrku, wrkv
    real(sp), parameter :: Inival=0.0_sp
    !

    !----------------------------------------------------------------------
    !----------Executing the code------------------------------------------
    !----------------------------------------------------------------------
    ! 
    !   This is currently assuming that vegetation types do not overlap
    ! 
    VEG % visc2d_r_veg=Inival
    VEG % visc3d_r_veg=Inival
    visc3d_r_veg_loc=Inival

    DO iveg=1,NVEG
       DO k=kbm1,1,-1  
          DO i=1,m
             ! cycle the rest calculation if no vegetation covered
             if(VEG%plant(i,iveg,phght)==0)cycle

             maxgradvegdens=VEG % plant(i,iveg,pdens)
             maxgradvegdens=MAX(maxgradvegdens,VEG % plant(i,iveg,pdens))
          END DO
       ENDDO
    ENDDO
    !
    DO iveg=1,NVEG
       DO k=kbm1,1,-1 
          DO i=1,m
             ! cycle the rest calculation if no vegetation covered
             if(VEG%plant(i,iveg,phght)==0)cycle

             cff=VegHMixCoef(iveg)*VEG % Lveg(i,k)/maxgradvegdens
             max_loc = maxloc(vx(nbsn(i,1:ntsn(i))))
             min_loc = minloc(vx(nbsn(i,1:ntsn(i))))
             gradu=ABS(VEG % plant(nbsn(i,max_loc(1)),iveg,pdens)-VEG % plant(nbsn(i,min_loc(1)),iveg,pdens))
             wrku=cff*gradu

             max_loc = maxloc(vy(nbsn(i,1:ntve(i))))
             min_loc = minloc(vy(nbsn(i,1:ntve(i))))
             gradv=ABS(VEG % plant(nbsn(i,max_loc(1)),iveg,pdens)-VEG % plant(nbsn(i,min_loc(1)),iveg,pdens))
             wrkv=cff*gradv
             !
             visc3d_r_veg_loc=SQRT(wrku**2+wrkv**2)
             !
             !   Adding from multiple veg types
             ! 
             VEG % visc3d_r_veg(i,k)=visc3d_r_veg_loc+VEG % visc3d_r_veg(i,k)
             !
             !   For overlapping veg types, maximum limit 
             !
             VEG % visc3d_r_veg(i,k)=MIN(VEG % visc3d_r_veg(i,k),              &
                  cff*maxgradvegdens)
          END DO
       END DO
    END DO
# if defined (MULTIPROCESSOR)
    if(par)then
       call aexchange(nc,myid,nprocs,VEG % visc3d_r_veg)
    endif
# endif
  END SUBROUTINE Calc_vegetation_hmixing

  SUBROUTINE LOAD_VEGETATION_FIELDS(PLANT,NVEG2,NVEGP2)
    USE CONTROL
    IMPLICIT NONE
    REAL(SP),ALLOCATABLE  :: PLANT(:,:,:)
    INTEGER,INTENT(IN)    :: NVEG2
    INTEGER,INTENT(IN)    :: NVEGP2
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    TYPE(NCDIM),  POINTER :: DIM3
    integer status,I,IERR

    LOGICAL FOUND

    ! FIND THE VEGETATION FIELDS FILE OBJECT
    NCF => FIND_FILE(FILEHEAD,trim(VEGETATION_FIELDS_FILE),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND VEGETATION FIELDS FILE OBJECT",&
         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))

    DIM1 => FIND_DIM(NCF,'node',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND VEGETATION_FIELDS_FILE DIMENSION 'node' in:",&
         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))
    IF (DIM1%DIM /= MGL)CALL FATAL_ERROR &
         & ("Dimension 'node' in the VEGETATION_FIELDS_FILE does not match MGL for this model?",&
         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))

    DIM2 => FIND_DIM(NCF,'nveg',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND VEGETATION_FIELDS_FILE DIMENSION 'nveg' in:",&
         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))
    IF (DIM2%DIM /= nveg2)CALL FATAL_ERROR &
         & ("Dimension 'nveg' in the VEGETATION_FIELDS_FILE does not match MGL for this model?",&
         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))

    DIM3 => FIND_DIM(NCF,'nvegp',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND VEGETATION_FIELDS_FILE DIMENSION 'nveg' in:",&
         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))
    IF (DIM3%DIM /= nvegp2)CALL FATAL_ERROR &
         & ("Dimension 'nvegp' in the VEGETATION_FIELDS_FILE does not match MGL for this model?",&
         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))

    ! FIND THE 'plant' variable
    VAR => FIND_VAR(NCF,'plant',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND VEGETATION_FIELDS_FILE VARIABLE 'plant' in:",&
         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))

    CALL NC_CONNECT_AVAR(VAR,plant)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)



!!$    ! FIND THE 'density' variable
!!$    VAR => FIND_VAR(NCF,'plant_density',FOUND)
!!$    IF(.not. FOUND) CALL FATAL_ERROR &
!!$         & ("COULD NOT FIND VEGETATION_FIELDS_FILE VARIABLE 'plant_density' in:",&
!!$         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))
!!$
!!$    CALL NC_CONNECT_AVAR(VAR,plant(:,:,1))
!!$    CALL NC_READ_VAR(VAR)
!!$    CALL NC_DISCONNECT(VAR)
!!$
!!$    ! FIND THE 'height' variable
!!$    VAR => FIND_VAR(NCF,'plant_height',FOUND)
!!$    IF(.not. FOUND) CALL FATAL_ERROR &
!!$         & ("COULD NOT FIND VEGETATION_FIELDS_FILE VARIABLE 'plant_height' in:",&
!!$         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))
!!$
!!$    CALL NC_CONNECT_AVAR(VAR,plant(:,:,2))
!!$    CALL NC_READ_VAR(VAR)
!!$    CALL NC_DISCONNECT(VAR)
!!$
!!$    ! FIND THE 'diameter' variable
!!$    VAR => FIND_VAR(NCF,'plant_diameter',FOUND)
!!$    IF(.not. FOUND) CALL FATAL_ERROR &
!!$         & ("COULD NOT FIND VEGETATION_FIELDS_FILE VARIABLE 'plant_diameter' in:",&
!!$         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))
!!$
!!$    CALL NC_CONNECT_AVAR(VAR,plant(:,:,3))
!!$    CALL NC_READ_VAR(VAR)
!!$    CALL NC_DISCONNECT(VAR)
!!$
!!$    ! FIND THE 'thickness' variable
!!$    VAR => FIND_VAR(NCF,'plant_thickness',FOUND)
!!$    IF(.not. FOUND) CALL FATAL_ERROR &
!!$         & ("COULD NOT FIND VEGETATION_FIELDS_FILE VARIABLE 'plant_thickness' in:",&
!!$         & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))
!!$
!!$    CALL NC_CONNECT_AVAR(VAR,plant(:,:,4))
!!$    CALL NC_READ_VAR(VAR)
!!$    CALL NC_DISCONNECT(VAR)
!!$
!!$    IF(VEG_BIOMASS)THEN
!!$       ! FIND THE 'thickness' variable
!!$       VAR => FIND_VAR(NCF,'plant_biomass_above_ground',FOUND)
!!$       IF(.not. FOUND) CALL FATAL_ERROR &
!!$            & ("COULD NOT FIND VEGETATION_FIELDS_FILE VARIABLE 'plant_biomass_above_ground' in:",&
!!$            & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))
!!$
!!$       CALL NC_CONNECT_AVAR(VAR,plant(:,:,5))
!!$       CALL NC_READ_VAR(VAR)
!!$       CALL NC_DISCONNECT(VAR)
!!$
!!$       ! FIND THE 'thickness' variable
!!$       VAR => FIND_VAR(NCF,'plant_biomass_below_ground',FOUND)
!!$       IF(.not. FOUND) CALL FATAL_ERROR &
!!$            & ("COULD NOT FIND VEGETATION_FIELDS_FILE VARIABLE 'plant_biomass_below_ground' in:",&
!!$            & "FILE NAME: "//TRIM(VEGETATION_FIELDS_FILE))
!!$
!!$       CALL NC_CONNECT_AVAR(VAR,plant(:,:,6))
!!$       CALL NC_READ_VAR(VAR)
!!$       CALL NC_DISCONNECT(VAR)
!!$
!!$    END IF

  END SUBROUTINE LOAD_VEGETATION_FIELDS

#endif

END MODULE mod_vegetation
