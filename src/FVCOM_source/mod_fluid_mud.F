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
! FVCOM Fluid Mud Module 
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
! Contact:      J. Ge 
!               East China Normal University, Shanghai China
!
! Based on mathematical model for the fluid mud layer given by Wang
!  and Winterwerp (1992).
!
! Comments:     Fluid-Mud Layer Dynamics Module 
!
! Current FVCOM dependency
!

!
! History

!=======================================================================
Module Mod_Fluid_Mud  
#  if defined (FLUID_MUD)
Use Mod_Time
Use Mod_Par
Use Mod_Prec 
Use Mod_Types
Use Mod_wd
Use Lims
Use Control
Use all_vars
implicit none 

save

 TYPE(TIME) :: ExtTime_fml
 TYPE(TIME) :: IMDTE_fml 
 Type(TIME) :: RKTime_FML
 REAL(SP)   :: DTE_fml        !!EXTERNAL TIME STEP (Seconds)
 REAL(DP)  EXTSTEP_SECONDS_fml

 real(sp),public :: cbed    !=> sediment concentration of the bed [g/L]
 real(sp),public :: cmud    !=> sediment concentration of the fluid mud layer [g/L]
 real(sp),public :: fmud    !=> friction coefficient between consolidated bed and fluid mud layer [-]
 real(sp),public :: fwat    !=> friction coefficient between suspension layer and fluid mud layer [-]
 real(sp),public :: mers    !=> bulk erison coefficient [kg/m^2/s]
 real(sp),public :: rhosus  !=> density of the suspension layer [kg/m^3]
 real(sp),public :: rhomud  !=> density of the mud layer [kg/m^3]
 real(sp),public :: taubng  !=> Bingham yield stress [pa]
 real(sp),public :: vdew    !=> dewatering velocity [m/s]
 integer,public :: dte_ratio !=> ration between ocean and fluid mud models

 character(len=120) :: mudfile
 real(sp),parameter,public :: min_depth_fml = 0.05_sp

 PRIVATE EROSION
 PUBLIC

 REAL(SP), PUBLIC, ALLOCATABLE :: Settling(:)    !=> suspension settling at suspension-mud interface
 REAL(SP), PUBLIC, ALLOCATABLE :: Entrainment(:) !=> sediment entrainment at suspension-mud interface
 REAL(SP), PUBLIC, ALLOCATABLE :: BedErosion(:)  !=> sediment erosion at mud-bed interface
 REAL(SP), PUBLIC, ALLOCATABLE :: Dewatering(:)  !=> sediment dewatering and consolidation at mud-bed interface

 REAL(SP), PUBLIC, ALLOCATABLE :: tempx(:)  !=> sediment dewatering and consolidation at mud-bed interface


!---------------2-d flow variable arrays at elements-------------------------------!
!
!
!      FML  => Fluid Mud Layer
!
!----------------------------------------------------------------------------------!
   REAL(SP), PUBLIC, ALLOCATABLE :: UA_FML(:)            !!VERTICALLY AVERAGED X-VELOC
   REAL(SP), PUBLIC, ALLOCATABLE :: VA_FML(:)            !!VERTICALLY AVERAGED Y-VELOC
   REAL(SP), PUBLIC, ALLOCATABLE :: UAF_FML(:)           !!UA FROM PREVIOUS RK STAGE 
   REAL(SP), PUBLIC, ALLOCATABLE :: VAF_FML(:)           !!VA FROM PREVIOUS RK STAGE 
   REAL(SP), PUBLIC, ALLOCATABLE :: UARK_FML(:)          !!UA FROM PREVIOUS TIMESTEP 
   REAL(SP), PUBLIC, ALLOCATABLE :: VARK_FML(:)          !!VA FROM PREVIOUS TIMESTEP 
   REAL(SP), PUBLIC, ALLOCATABLE :: UARD_FML(:)          !!UA AVERAGED OVER EXTERNAL INT
   REAL(SP), PUBLIC, ALLOCATABLE :: VARD_FML(:)          !!VA AVERAGED OVER EXTERNAL INT
 
   REAL(SP), PUBLIC, ALLOCATABLE :: H1_FML(:)            !!BATHYMETRIC DEPTH   
   REAL(SP), PUBLIC, ALLOCATABLE :: D1_FML(:)            !!CURRENT DEPTH
   REAL(SP), PUBLIC, ALLOCATABLE :: DT1_FML(:)           !!DEPTH AT PREVIOUS TIME STEP
   REAL(SP), PUBLIC, ALLOCATABLE :: EL1_FML(:)           !!CURRENT SURFACE ELEVATION
   REAL(SP), PUBLIC, ALLOCATABLE :: ET1_FML(:)           !!SURFACE ELEVATION AT PREVIOUS TIME STEP
   REAL(SP), PUBLIC, ALLOCATABLE :: ELRK1_FML(:)         !!SURFACE ELEVATION AT BEGINNING OF RK INT 
   REAL(SP), PUBLIC, ALLOCATABLE :: DRK1_FML(:)          !!SURFACE ELEVATION AT BEGINNING OF RK INT 
   REAL(SP), PUBLIC, ALLOCATABLE :: ELF1_FML(:)          !!SURFACE ELEVATION STORAGE FOR RK INT
   REAL(SP), PUBLIC, ALLOCATABLE :: DTFA_FML(:)          !!ADJUSTED DEPTH FOR MASS CONSERVATION
!---------------2-d flow variable arrays at nodes----------------------------------!

   REAL(SP), PUBLIC, ALLOCATABLE :: SOURCE_FLUID_MUD(:) !! source/sink term of fluid mud layer 

   REAL(SP), PUBLIC, ALLOCATABLE :: H_FML(:)     !!BATHYMETRIC DEPTH   
   REAL(SP), PUBLIC, ALLOCATABLE :: D_FML(:)             !!CURRENT DEPTH (m)
   REAL(SP), PUBLIC, ALLOCATABLE :: D_FMLcm(:)           !!CURRENT DEPTH (cm)   
   REAL(SP), PUBLIC, ALLOCATABLE :: DT_FML(:)            !!DEPTH AT PREVIOUS TIME STEP
   REAL(SP), PUBLIC, ALLOCATABLE :: EL_FML(:)            !!CURRENT SURFACE ELEVATION
   REAL(SP), PUBLIC, ALLOCATABLE :: ET_FML(:)            !!SURFACE ELEVATION AT PREVIOUS TIME STEP
   REAL(SP), PUBLIC, ALLOCATABLE :: EGF_FML(:)           !!AVERAGE SURFACE ELEVATION OVER EXTERNAL INT
   REAL(SP), PUBLIC, ALLOCATABLE :: ELF_FML(:)           !!SURFACE ELEVATION STORAGE FOR RK INT
   REAL(SP), PUBLIC, ALLOCATABLE :: ELRK_FML(:)          !!SURFACE ELEVATION AT BEGINNING OF RK INT
   REAL(SP), PUBLIC, ALLOCATABLE :: DRK_FML(:)           !!SURFACE ELEVATION AT BEGINNING OF RK INT
   REAL(SP), PUBLIC, ALLOCATABLE :: WUBOT_FML(:)         !!BOTTOM FRICTION
   REAL(SP), PUBLIC, ALLOCATABLE :: WVBOT_FML(:)         !!BOTTOM FRICTION
   REAL(SP), PUBLIC, ALLOCATABLE :: TAUBM_FML(:)         !!BOTTOM FRICTION MAGNITUDE(Caution is Tau' [no Rho])
   REAL(SP), PUBLIC, ALLOCATABLE :: WUBOT_N_FML(:)       !!BOTTOM FRICTION ON NODES (Caution is Tau' [no Rho])
   REAL(SP), PUBLIC, ALLOCATABLE :: WVBOT_N_FML(:)       !!BOTTOM FRICTION ON NODES (Caution is Tau' [no Rho])
   REAL(SP), PUBLIC, ALLOCATABLE :: TAUBM_N_FML(:)       !!BOTTOM FRICTION MAGNITUDE ON NODES (Caution is Tau' [no Rho])
   REAL(SP), PUBLIC, ALLOCATABLE :: WUSURF_FML(:)        !!SURFACE FRICTION FOR INT
   REAL(SP), PUBLIC, ALLOCATABLE :: WVSURF_FML(:)        !!SURFACE FRICTION FOR INT
   REAL(SP), PUBLIC, ALLOCATABLE :: WUSURF2_FML(:)       !!SURFACE FRICTION FOR EXT
   REAL(SP), PUBLIC, ALLOCATABLE :: WVSURF2_FML(:)       !!SURFACE FRICTION FOR EXT
!-----------------------2-d flow fluxes--------------------------------------------!

   REAL(SP), PUBLIC, ALLOCATABLE :: PSTX_FML(:)           !!EXT MODE BAROTROPIC TERMS
   REAL(SP), PUBLIC, ALLOCATABLE :: PSTY_FML(:)           !!EXT MODE BAROTROPIC TERMS
   REAL(SP), PUBLIC, ALLOCATABLE :: ADVUA_FML(:) 
   REAL(SP), PUBLIC, ALLOCATABLE :: ADVVA_FML(:) 
   REAL(SP), PUBLIC, ALLOCATABLE :: ADX2D_FML(:) 
   REAL(SP), PUBLIC, ALLOCATABLE :: ADY2D_FML(:) 
   REAL(SP), PUBLIC, ALLOCATABLE :: DRX2D_FML(:) 
   REAL(SP), PUBLIC, ALLOCATABLE :: DRY2D_FML(:) 
   REAL(SP), PUBLIC, ALLOCATABLE :: ADVX_FML(:,:)      
   REAL(SP), PUBLIC, ALLOCATABLE :: ADVY_FML(:,:) 


   REAL(SP), PUBLIC, ALLOCATABLE :: CBC_FML(:)           !!BOTTOM FRICTION   
   REAL(SP), PUBLIC, ALLOCATABLE :: TM(:)          
   REAL(SP), PUBLIC, ALLOCATABLE :: DELTA_U(:)   
   REAL(SP), PUBLIC, ALLOCATABLE :: DELTA_V(:) 

!
!--Parameters for Wet/Dry Treatment                 
!

!-----variables controlling porosities through wet/dry determination----------------!
   
   INTEGER ,PUBLIC, ALLOCATABLE :: ISWETN_FML(:)  !!NODE POROSITY AT NODES FOR TIME N
   INTEGER ,PUBLIC, ALLOCATABLE :: ISWETC_FML(:)  !!CELL POROSITY AT CELLS FOR TIME N
   INTEGER ,PUBLIC, ALLOCATABLE :: ISWETNT_FML(:) !!NODE POROSITY AT NODES FOR TIME N-1 INTERNAL
   INTEGER ,PUBLIC, ALLOCATABLE :: ISWETCT_FML(:) !!CELL POROSITY AT CELLS FOR TIME N-1 INTERNAL
   INTEGER ,PUBLIC, ALLOCATABLE :: ISWETCE_FML(:) !!CELL POROSITY AT CELLS FOR TIME N-1 EXTERNAL

   REAL(SP) ,PUBLIC, ALLOCATABLE :: WETNODES_FML(:)  !!NODE POROSITY AT NODES FOR TIME N
   REAL(SP) ,PUBLIC, ALLOCATABLE :: WETCELLS_FML(:)  !!CELL POROSITY AT CELLS FOR TIME N
   REAL(SP) ,PUBLIC, ALLOCATABLE :: WETNODES_CUR(:)  !!NODE POROSITY AT NODES FOR TIME N
   REAL(SP) ,PUBLIC, ALLOCATABLE :: WETCELLS_CUR(:)  !!CELL POROSITY AT CELLS FOR TIME N
contains
!==============================================================================|
!    Allocate and Initialize Most Arrays                                       !
!==============================================================================|
SUBROUTINE INITIALIZE_FLUID_MUD
use control, only : msr,casename,input_dir,SEDIMENT_MODEL_FILE
Use Lims
IMPLICIT NONE
INTEGER :: IERR
REAL(SP) :: max_depth
REAL(SP) :: SBUF,RBUF1,RBUF2,RBUF3
INTEGER :: DEST, SOURCE
# if defined(MULTIPROCESSOR)
   INTEGER :: STAT(MPI_STATUS_SIZE)
# endif
!-----------------------determine the maximum of the depth-------------------------!

   SBUF = 0.0_DP
   SBUF  = MAXVAL(H(1:M))
   RBUF2 = SBUF
# if defined (MULTIPROCESSOR)
    IF(PAR)CALL MPI_ALLREDUCE(SBUF,RBUF2,1,MPI_F,MPI_MAX,MPI_FVCOM_GROUP,IERR)
# endif
  max_depth=RBUF2

  !print*,h(1),real(h(1),dp),max_depth,real(max_depth,dp)

!---------------2-d flow variable arrays at elements-------------------------------!

   ALLOCATE(UA_FML(0:NT))            ;UA_FML        = ZERO  !!VERTICALLY AVERAGED X-VELOC
   ALLOCATE(VA_FML(0:NT))            ;VA_FML        = ZERO  !!VERTICALLY AVERAGED Y-VELOC
   ALLOCATE(UAF_FML(0:NT))           ;UAF_FML       = ZERO  !!UA FROM PREVIOUS RK STAGE 
   ALLOCATE(VAF_FML(0:NT))           ;VAF_FML       = ZERO  !!VA FROM PREVIOUS RK STAGE 

   ALLOCATE(UARK_FML(0:NT))          ;UARK_FML      = ZERO  !!UA FROM PREVIOUS TIMESTEP 
   ALLOCATE(VARK_FML(0:NT))          ;VARK_FML      = ZERO  !!VA FROM PREVIOUS TIMESTEP 
   ALLOCATE(UARD_FML(0:NT))          ;UARD_FML      = ZERO  !!UA AVERAGED OVER EXTERNAL INT
   ALLOCATE(VARD_FML(0:NT))          ;VARD_FML      = ZERO  !!VA AVERAGED OVER EXTERNAL INT
   ALLOCATE(H1_FML(0:NT))            ;H1_FML        = H1-max_depth !!BATHYMETRIC DEPTH   
   ALLOCATE(D1_FML(0:NT))            ;D1_FML        = ZERO !!DEPTH
   ALLOCATE(DT1_FML(0:NT))           ;DT1_FML       = ZERO  !!DEPTH  
   ALLOCATE(EL1_FML(0:NT))           ;EL1_FML       = ZERO  !!SURFACE ELEVATION
   ALLOCATE(ELF1_FML(0:NT))          ;ELF1_FML      = ZERO  !!SURFACE ELEVATION
   ALLOCATE(DTFA_FML(0:MT))          ;DTFA_FML      = ZERO  !!ADJUSTED DEPTH FOR MASS CONSERVATION
   ALLOCATE(ET1_FML(0:NT))           ;ET1_FML       = ZERO  !!SURFACE ELEVATION
   ALLOCATE(ELRK1_FML(0:NT))         ;ELRK1_FML     = ZERO  !!SURFACE ELEVATION
   ALLOCATE(DRK1_FML(0:NT))          ;DRK1_FML     = ZERO  !!SURFACE ELEVATION
!---------------2-d sediment variable arrays at nodes----------------------------------!
   ALLOCATE(Settling(0:MT))          ;Settling    = ZERO
   ALLOCATE(Entrainment(0:MT))       ;Entrainment = ZERO
   ALLOCATE(BedErosion(0:MT))        ;BedErosion  = ZERO
   ALLOCATE(Dewatering(0:MT))        ;Dewatering  = ZERO
   ALLOCATE(SOURCE_FLUID_MUD(0:MT))  ;SOURCE_FLUID_MUD  = ZERO !! SOURCE/SINK TERM

   ALLOCATE(tempx(0:MT))       ;tempx = ZERO

!---------------2-d flow variable arrays at nodes----------------------------------!
   ALLOCATE(H_FML(0:MT))             ;H_FML    = H-max_depth !!BATHYMETRIC DEPTH   
   ALLOCATE(D_FML(0:MT))             ;D_FML    = ZERO        !!DEPTH (m)  
   ALLOCATE(D_FMLcm(0:MT))           ;D_FMLcm  = ZERO        !!DEPTH (cm)
   ALLOCATE(DT_FML(0:MT))            ;DT_FML   = ZERO        !!DEPTH   
   ALLOCATE(ET_FML(0:MT))            ;ET_FML   = ZERO        !!SURFACE ELEVATION
   ALLOCATE(EL_FML(0:MT))            ;EL_FML   = ZERO        !!SURFACE ELEVATION
   ALLOCATE(ELF_FML(0:MT))           ;ELF_FML  = ZERO       !!SURFACE ELEVATION
   ALLOCATE(ELRK_FML(0:MT))          ;ELRK_FML = ZERO        !!SURFACE ELEVATION
   ALLOCATE(DRK_FML(0:MT))           ;DRK_FML  = ZERO        !!SURFACE ELEVATION
   ALLOCATE(WUSURF2_FML(0:NT))       ;WUSURF2_FML= ZERO     !!SURFACE FRICTION FOR EXT
   ALLOCATE(WVSURF2_FML(0:NT))       ;WVSURF2_FML= ZERO     !!SURFACE FRICTION FOR EXT

   ALLOCATE(WUBOT_FML(0:NT))         ;WUBOT_FML  = ZERO     !!BOTTOM FRICTION
   ALLOCATE(WVBOT_FML(0:NT))         ;WVBOT_FML  = ZERO     !!BOTTOM FRICTION
   ALLOCATE(TAUBM_FML(0:NT))         ;TAUBM_FML  = ZERO     !!BOTTOM FRICTION
   ALLOCATE(WUBOT_N_FML(0:MT))       ;WUBOT_N_FML  = ZERO   !!U-Component bottom shear stress on nodes  
   ALLOCATE(WVBOT_N_FML(0:MT))       ;WVBOT_N_FML  = ZERO   !!V-Component bottom shear stress on nodes
   ALLOCATE(TAUBM_N_FML(0:MT))       ;TAUBM_N_FML  = ZERO   !!Magnitude bottom shear stress on nodes
   ALLOCATE(WUSURF_FML(0:NT))        ;WUSURF_FML = ZERO     !!SURFACE FRICTION FOR INT
   ALLOCATE(WVSURF_FML(0:NT))        ;WVSURF_FML = ZERO     !!SURFACE FRICTION FOR INT

!-----------------------2-d flow fluxes---------------------------------------------!

   ALLOCATE(CBC_FML(0:NT))           ;CBC_FML   = ZERO 
   ALLOCATE(TM(0:NT))                ;TM        = ZERO 
   ALLOCATE(DELTA_U(0:NT))           ;DELTA_U   = ZERO 
   ALLOCATE(DELTA_V(0:NT))           ;DELTA_V   = ZERO 

   ALLOCATE(PSTX_FML(0:NT))          ;PSTX_FML  = ZERO       !!EXT MODE BAROTROPIC TERMS
   ALLOCATE(PSTY_FML(0:NT))          ;PSTY_FML  = ZERO       !!EXT MODE BAROTROPIC TERMS
   ALLOCATE(ADVUA_FML(0:NT))         ;ADVUA_FML = ZERO 
   ALLOCATE(ADVVA_FML(0:NT))         ;ADVVA_FML = ZERO
   ALLOCATE(ADX2D_FML(0:NT))         ;ADX2D_FML = ZERO
   ALLOCATE(ADY2D_FML(0:NT))         ;ADY2D_FML = ZERO
   ALLOCATE(DRX2D_FML(0:NT))         ;DRX2D_FML = ZERO
   ALLOCATE(DRY2D_FML(0:NT))         ;DRY2D_FML = ZERO
   ALLOCATE(ADVX_FML(0:NT,KB))       ;ADVX_FML  = ZERO 
   ALLOCATE(ADVY_FML(0:NT,KB))       ;ADVY_FML  = ZERO 



   CALL SETUP_FLUID_MUD(SEDIMENT_MODEL_FILE)

   ExtTime_fml = StartTime
   EXTSTEP_SECONDS_fml = EXTSTEP_SECONDS/dte_ratio
   IMDTE_fml = SECONDS2TIME(EXTSTEP_SECONDS_fml)
   DTE_fml=seconds(IMDTE_fml)

END SUBROUTINE INITIALIZE_FLUID_MUD
!==============================================================================!


!==============================================================================!
 Subroutine Setup_Fluid_Mud(fname)
 use control, only : msr,casename,input_dir
 use lims,    only : m
 Use Mod_Utils, only: pstop
 Use Input_Util
 implicit none
 character(len=80)    :: fname
 logical             :: fexist
 integer linenum,i,k1,iscan

 if(dbg_set(dbg_sbr)) write(ipt,*) "Start: setup_fluid_mud" 

  if(fname == "none" .or. fname == "NONE" .or. fname == "None")then
    mudfile = "./"//trim(input_dir)//"/"//trim(casename)//'_sediment.inp'
  else
    mudfile = "./"//trim(input_dir)//"/"//trim(fname)
  endif
    
  inquire(file=trim(mudfile),exist=fexist)
  if(.not.fexist)then
    write(*,*)'fluid mud parameter file: ',trim(mudfile),' does not exist'
    write(*,*)'stopping'
    call pstop
  end if

 !read in sediment concentration of the bed 
  Call Get_Val(cbed,mudfile,'CBED',line=linenum,echo=.false.)  

 !read in sediment concentration of the fluid mud layer 
  Call Get_Val(cmud,mudfile,'CMUD',line=linenum,echo=.false.)  

 !read in friction coefficient between consolidated bed and fluid mud layer 
  Call Get_Val(fmud,mudfile,'FMUD',line=linenum,echo=.false.)  

 !read in friction coefficient between suspension layer and fluid mud layer 
  Call Get_Val(fwat,mudfile,'FWAT',line=linenum,echo=.false.)  

 !read in bulk erosion coefficient
  Call Get_Val(mers,mudfile,'MERS',line=linenum,echo=.false.) 

 !read in density of the suspension layers
  Call Get_Val(rhosus,mudfile,'RHOSUS',line=linenum,echo=.false.) 

 !read in density of the fluid mud layers
  Call Get_Val(rhomud,mudfile,'RHOMUD',line=linenum,echo=.false.)

 !read in Bingham yield stress
  Call Get_Val(taubng,mudfile,'TAUBNG',line=linenum,echo=.false.)

 !read in dewatering velocity
  Call Get_Val(vdew,mudfile,'VDEW',line=linenum,echo=.false.)

 !read in ration between ocean and fluid mud models
  Call Get_Val(dte_ratio,mudfile,'DTE_RATIO',line=linenum,echo=.false.)

!------  ALLOCATE AND INITIALIZE WET/DRY TREATMENT ARRAYS----------------------|
  Call ALLOC_WD_DATA_FML

!------  INITIALIZE ARRAYS USED FOR WET/DRY TREATMENT--------------------------| 
  CALL SET_WATER_DEPTH_FML

  Call SET_WD_DATA_FML

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: setup_fluid_mud" 
 
  Return
 
End Subroutine Setup_Fluid_Mud

!==============================================================================!


SUBROUTINE INTERNAL_STEP_FLUID_MUD
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_OBCS
  USE MOD_TIME
  USE EQS_OF_STATE
  USE MOD_WD
  USE MOD_PAR
  USE MOD_NORTHPOLE


  IMPLICIT NONE
  INTEGER :: I,J,K
  REAL(SP) :: UTMP,VTMP

# if defined (SEDIMENT)
  REAL(DP) :: DPDAYS
  REAL(SP) :: SPDAYS
# endif
  integer :: i1,i2,IRA

!------------------------------------------------------------------------------|

  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "Start: internal_step_fluid_mud"

!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
  RAMP = 0.0_SP
  IF(IRAMP /= 0) THEN
     RAMP=TANH(real(IINT,sp)/real(IRAMP,sp))
  ELSE
     RAMP = 1.0_SP
  END IF

  DO IRA=1,DTE_RATIO
  
 !----SPECIFY THE BOTTOM ROUGHNESS AND CALCULATE THE BOTTOM STRESSES------------!
    CALL FRICTION_FLUID_MUD
 ! CALL BOTTOM_ROUGHNESS_FML

  ! New Open Boundary Condition ----4

  !==============================================================================!
  !  LOOP OVER EXTERNAL TIME STEPS                                               !
  !==============================================================================!
    DO IEXT=1,ISPLIT

       IF (DBG_SET(DBG_SBRIO)) WRITE(IPT,*) "/// EXT SETP: ",IEXT

       ExtTime_fml = ExtTime_fml + IMDTE_fml

       CALL EXTERNAL_STEP_FLUID_MUD
     
     
    END DO

  !
  !----SHIFT SEA SURFACE ELEVATION AND DEPTH TO CURRENT TIME LEVEL---------------!
  !
    ET_FML  = EL_FML  
    DT_FML  = D_FML 
    ET1_FML = EL1_FML
    DT1_FML = D1_FML

    CALL WD_UPDATE_FML(3)
  END DO

  D_FMLcm=D_FML*100.0 ! convert fml depth to centmeter.

  WETCELLS_FML = ISWETC_FML*1.0
  WETNODES_FML = ISWETN_FML*1.0
  WETCELLS_CUR = ISWETC*1.0
  WETNODES_CUR = ISWETN*1.0

!  if(nlid(287)>0)then
!    print*,EL_FML(nlid(287)),ELF_FML(nlid(287)),H_FML(nlid(287))
!    print*,D_FML(nlid(287)),WETNODES_CUR(nlid(287)),WETNODES_FML(nlid(287)),Entrainment(nlid(287))
!    print*,Settling(nlid(287)),BedErosion(nlid(287)),dewatering(nlid(287)),Source_Fluid_Mud(nlid(287))
!  endif
  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "End: internal_step_fluid_mud"

END SUBROUTINE INTERNAL_STEP_FLUID_MUD




SUBROUTINE EXTERNAL_STEP_FLUID_MUD

  USE MOD_UTILS
  USE ALL_VARS
  USE MOD_TIME
  USE MOD_OBCS
  USE MOD_WD
# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif

  IMPLICIT NONE
  REAL(SP) :: TMP
  INTEGER :: K, I, J, JN, J1,i1,i2

!------------------------------------------------------------------------------|

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: external_step_fluid_mud"
 

!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
  IF(IRAMP /= 0) THEN
     TMP = real(IINT-1,sp)+real(IEXT,sp)/real(ISPLIT,sp)
     RAMP=TANH(TMP/real(IRAMP,sp))
  ELSE
     RAMP = 1.0_SP
  END IF


!
!------SAVE VALUES FROM CURRENT TIME STEP--------------------------------------!
!
  ELRK1_FML = EL1_FML
  ELRK_FML  = EL_FML
  UARK_FML  = UA_FML
  VARK_FML  = VA_FML

!
!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
!

  DO K=1,4
   
     RKTIME_FML = ExtTime_fml + IMDTE_fml * (ALPHA_RK(K) - 1.0_DP)

!     CALL PRINT_REAL_TIME(RKTIME,IPT,"RUNGE-KUTTA")

     
     
!FREE SURFACE AMPLITUDE UPDATE  --> ELF

     CALL EXTEL_EDGE_FLUID_MUD(K)

# if defined (MULTIPROCESSOR)
     IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ELF_FML)
# endif
     
     DO I=1,IBCN(1)
        JN = OBC_LST(1,I)
        J=I_OBC_N(JN)
        ELF_FML(J)=ELRK_FML(J)+ALPHA_RK(K)*(ELF_FML(J)-ELRK_FML(J))
     END DO
     
     ! DAVID ADDED THIS EXCHANGE CALL:
     ! IT SEEMS LIKELY THAT THE HALO VALUES OF ELF WILL BE USED
     ! BEFORE THEY ARE SET CORRECTLY OTHERWISE
# if defined (MULTIPROCESSOR)
     IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ELF_FML)
# endif

     CALL N2E2D(ELF_FML,ELF1_FML)

     CALL WET_JUDGE_FML
          
     CALL FLUX_OBN_FML(K)

    !CALCULATE ADVECTIVE, DIFFUSIVE, AND BAROCLINIC MODES --> UAF ,VAF
     CALL ADVAVE_EDGE_GCN_FLUID_MUD(ADVUA_FML,ADVVA_FML)           !Compute Ext Mode Adv/Diff
     CALL EXTUV_EDGE_FLUID_MUD(K)
  
     CALL BCOND_GCN_2_FML

# if defined (MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,ELF_FML)
# endif
  
     EL_FML  = ELF_FML
     EL1_FML = ELF1_FML

     !!INTERPOLATE DEPTH FROM NODE-BASED TO ELEMENT-BASED VALUES
     CALL N2E2D(EL_FML,EL1_FML)
     
     !UPDATE DEPTH AND VERTICALLY AVERAGED VELOCITY FIELD
     D_FML   = H_FML + EL_FML
     D1_FML  = H1_FML + EL1_FML
     UA_FML  = UAF_FML
     VA_FML  = VAF_FML
     DTFA_FML = D_FML

     !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
# if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UA_FML,VA_FML,D1_FML)
# endif

     !UPDATE WET/DRY FACTORS
     CALL WD_UPDATE_FML(1)

  END DO     !! END RUNGE-KUTTA LOOP
  

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: external_step_fluid_mud"

END SUBROUTINE EXTERNAL_STEP_FLUID_MUD



!==============================================================================|
!  CALCULATE FLUXES OF FREE SURFACE ELEVATION (CONTINUITY) EQUATION            |
!==============================================================================|
   SUBROUTINE EXTEL_EDGE_FLUID_MUD(K)       
!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS

#  if defined (SPHERICAL)
   USE MOD_NORTHPOLE
#  endif


   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP) :: XFLUX(0:MT)
   REAL(SP) :: DIJ,UIJ,VIJ,DTK,UN,EXFLUX
   INTEGER  :: I,J,I1,IA,IB,JJ,J1,J2

!==============================================================================|
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: extel_edge_fluid_mud",K
!----------INITIALIZE FLUX ARRAY ----------------------------------------------!

   XFLUX = 0.0_SP

!---------ACCUMULATE FLUX BY LOOPING OVER CONTROL VOLUME HALF EDGES------------!

   DO I=1,NCV
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     DIJ = D1_FML(I1)

     UIJ = UA_FML(I1)
     VIJ = VA_FML(I1)
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + VIJ*DLTXE(I))  

     XFLUX(IA) = XFLUX(IA)-EXFLUX
     XFLUX(IB) = XFLUX(IB)+EXFLUX
   END DO

#  if defined (SPHERICAL)
   CALL EXTEL_EDGE_XY_FML(K,XFLUX)
#  endif

   
!--ADD SOURCE/SINK TERM OF FLUID MUD-------------------------------------------------------!
!   SOURCE_FLUID_MUD = 0.0
!   do i=1,M
!     if(ngid(i)==3206)SOURCE_FLUID_MUD(i) = 0.00001
!   end do

   XFLUX = XFLUX - SOURCE_FLUID_MUD*ART1


!--SAVE ACCUMULATED FLUX ON OPEN BOUNDARY NODES AND ZERO OUT OPEN BOUNDARY FLUX!

   IF(IOBCN > 0) THEN  
     DO I=1,IOBCN
       XFLUX_OBCN(I)=XFLUX(I_OBC_N(I))
       XFLUX(I_OBC_N(I)) = 0.0_SP
     END DO
   END IF


!----------PERFORM UPDATE ON ELF-----------------------------------------------!

   DTK = ALPHA_RK(K)*DTE_fml
   ELF_FML = ELRK_FML - DTK*XFLUX/ART1
 
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: extel_edge_fluid_mud"

   RETURN
   END SUBROUTINE EXTEL_EDGE_FLUID_MUD
!==============================================================================|






!==============================================================================|
!     ACCUMLATE FLUXES FOR EXTERNAL MODE                                       |
!==============================================================================|

   SUBROUTINE EXTUV_EDGE_FLUID_MUD(K)       
!==============================================================================|
   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_WD

   USE MOD_NORTHPOLE


   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:NT) :: RESX,RESY,TMP
   REAL(SP) :: UAFT,VAFT
   INTEGER  :: I

!==============================================================================|

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: extuv_edge_fluid_mud"

!
!--ACCUMULATE RESIDUALS FOR EXTERNAL MODE EQUATIONS----------------------------|
!
   UAFT = UAF_FML(0)
   VAFT = VAF_FML(0)

   ! THIS APPEARS TO BE TO PREVENT DIVISION BY ZERO, BUT IT IS A
   ! STRANGE WAY TO DO IT!
   H1_FML(0)= H1_FML(1)

   DO I=1,NT
#    if defined (WET_DRY)
     IF( ISWETCE_FML(I)*ISWETC_FML(I) == 1.AND.((WETTING_DRYING_ON.AND.ISWETCE(I)*ISWETC(I) == 1).OR.WETTING_DRYING_ON==.FALSE.))THEN
#    else
     IF( ISWETCE_FML(I)*ISWETC_FML(I) == 1)THEN
#    endif

       RESX(I) = ADX2D_FML(I)+ADVUA_FML(I)+DRX2D_FML(I)+PSTX_FML(I)&
            &-COR(I)*VA_FML(I)*D1_FML(I)*ART(I)-(WUSURF2_FML(I)&
            &+WUBOT_FML(I))*ART(I)

       RESY(I) = ADY2D_FML(I)+ADVVA_FML(I)+DRY2D_FML(I)+PSTY_FML(I)&
            &+COR(I)*UA_FML(I)*D1_FML(I)*ART(I)-(WVSURF2_FML(I)&
            &+WVBOT_FML(I))*ART(I)

#  if defined (SPHERICAL)
       RESX(I) = RESX(I)     &
                 -UA_FML(I)*VA_FML(I)/REARTH*TAN(DEG2RAD*YC(I))*D1_FML(I)*ART(I)
       RESY(I) = RESY(I)     &
                 +UA_FML(I)*UA_FML(I)/REARTH*TAN(DEG2RAD*YC(I))*D1_FML(I)*ART(I)
#  endif

!
!--UPDATE----------------------------------------------------------------------|
!

       UAF_FML(I)  =  (UARK_FML(I)*(H1_FML(I)+ELRK1_FML(I))&
            &-ALPHA_RK(K)*DTE_fml*RESX(I)/ART(I))/(H1_FML(I)+ELF1_FML(I))
       VAF_FML(I)  =  (VARK_FML(I)*(H1_FML(I)+ELRK1_FML(I))&
            &-ALPHA_RK(K)*DTE_fml*RESY(I)/ART(I))/(H1_FML(I)+ELF1_FML(I))


     ELSE
       UAF_FML(I) = 0.0_SP
       VAF_FML(I) = 0.0_SP
     END IF

   END DO

#  if defined (SPHERICAL)
   CALL EXTUV_EDGE_XY_FML(K)
#  endif
   
   VAF_FML(0) = VAFT
   UAF_FML(0) = UAFT

!
!--ADJUST EXTERNAL VELOCITY IN SPONGE REGION-----------------------------------|
!
   UAF_FML = UAF_FML-CC_SPONGE*UAF_FML
   VAF_FML = VAF_FML-CC_SPONGE*VAF_FML

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: extuv_edge_fluid_mud"

   END SUBROUTINE EXTUV_EDGE_FLUID_MUD
!==============================================================================|




!==============================================================================|
!   CALCULATE CONVECTION AND DIFFUSION FLUXES FOR EXTERNAL MODE                !
!==============================================================================|
   SUBROUTINE ADVAVE_EDGE_GCN_FLUID_MUD(XFLUX,YFLUX)
!==============================================================================|

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE BCS
   USE MOD_OBCS
   USE MOD_WD

   IMPLICIT NONE
   INTEGER  :: I,J,K,IA,IB,J1,J2,K1,K2,K3,I1,I2,II
   REAL(SP) :: DIJ,ELIJ,XIJ,YIJ,UIJ,VIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN_TMP
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP
   REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)
   REAL(SP) :: FACT,FM1,ISWETTMP
   INTEGER  :: STAT_MF    
   REAL(SP) :: TPA,TPB

#  if defined (SPHERICAL)
   REAL(DP) :: XTMP,XTMP1
#  endif      

#  if defined (LIMITED_NO)
   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
#  else
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UALFA,VALFA
   REAL(SP) :: UALFA_TMP,VALFA_TMP
   INTEGER :: ERROR
   REAL(SP) :: EPS
#  endif

!#  if defined (THIN_DAM)
   REAL(SP) :: A1UIA1,A1UIA2,A1UIA3,A1UIA4,A2UIA1,A2UIA2,A2UIA3,A2UIA4
   REAL(SP) :: A1UIB1,A1UIB2,A1UIB3,A1UIB4,A2UIB1,A2UIB2,A2UIB3,A2UIB4
   INTEGER  :: J11,J12,J21,J22,E1,E2,ISBCE1,ISBC_TMP,IB_TMP
   LOGICAL  :: ISMATCH
!#  endif

   REAL(SP) :: BTPS
   REAL(SP) :: U_TMP,V_TMP,UAC_TMP,VAC_TMP,WUSURF_TMP,WVSURF_TMP,WUBOT_TMP,WVBOT_TMP,UAF_TMP,VAF_TMP 

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: advave_edge_gcn_fluid_mud"

!------------------------------------------------------------------------------!

   SELECT CASE(HORIZONTAL_MIXING_TYPE)
   CASE ('closure')
      FACT = 1.0_SP
      FM1  = 0.0_SP
   CASE('constant')
      FACT = 0.0_SP
      FM1  = 1.0_SP
   CASE DEFAULT
      CALL FATAL_ERROR("UNKNOW HORIZONTAL MIXING TYPE:",&
           & TRIM(HORIZONTAL_MIXING_TYPE) )
   END SELECT

!
!-------------------------INITIALIZE FLUXES------------------------------------!
!
   XFLUX = 0.0_SP
   YFLUX = 0.0_SP
   PSTX_FML  = 0.0_SP
   PSTY_FML  = 0.0_SP

!
!-------------------------ACCUMULATE FLUX OVER ELEMENT EDGES-------------------!
!
#  if !defined (LIMITED_NO)
   ALLOCATE(UIJ1(NE),VIJ1(NE),UIJ2(NE),VIJ2(NE),STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("The arrays UIJ1,VIJ1,UIJ2 and VIJ2 can not be allocated.")
   UIJ1=0.0_SP;VIJ1=0.0_SP;UIJ2=0.0_SP;VIJ2=0.0_SP
   
   ALLOCATE(UALFA(0:NT),VALFA(0:NT),STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("The arrays UALFA,VALFA can not be allocated.")
   UALFA=1.0_SP;VALFA=1.0_SP
   
   ALLOCATE(FXX(NE),FYY(NE),STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("The arrays FXX,FYY can not be allocated.")
   FXX=0.0_SP;FYY=0.0_SP

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

     DIJ=0.5_SP*(D_FML(J1)+D_FML(J2))

#    if defined (WET_DRY)
     IF((ISWETCE_FML(IA)*ISWETC_FML(IA) == 1 .OR. ISWETCE_FML(IB)*ISWETC_FML(IB) == 1).AND.((WETTING_DRYING_ON.AND.(ISWETCE(IA)*ISWETC(IA) == 1 .OR. ISWETCE(IB)*ISWETC(IB) == 1)).OR.WETTING_DRYING_ON==.FALSE.))THEN
#    else
     IF(ISWETCE_FML(IA)*ISWETC_FML(IA) == 1 .OR. ISWETCE_FML(IB)*ISWETC_FML(IB) == 1)THEN
#    endif
!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)

     IB_TMP = IB

     A1UIA1 = A1U(IA,1)
     A1UIA2 = A1U(IA,2)
     A1UIA3 = A1U(IA,3)
     A1UIA4 = A1U(IA,4)
     A2UIA1 = A2U(IA,1)
     A2UIA2 = A2U(IA,2)
     A2UIA3 = A2U(IA,3)
     A2UIA4 = A2U(IA,4)
!---------------------------------------------------------------
     COFA1=A1UIA1*UA_FML(IA)+A1UIA2*UA_FML(K1)+A1UIA3*UA_FML(K2)+A1UIA4*UA_FML(K3)
     COFA2=A2UIA1*UA_FML(IA)+A2UIA2*UA_FML(K1)+A2UIA3*UA_FML(K2)+A2UIA4*UA_FML(K3)
     COFA5=A1UIA1*VA_FML(IA)+A1UIA2*VA_FML(K1)+A1UIA3*VA_FML(K2)+A1UIA4*VA_FML(K3)
     COFA6=A2UIA1*VA_FML(IA)+A2UIA2*VA_FML(K1)+A2UIA3*VA_FML(K2)+A2UIA4*VA_FML(K3)

#    if defined (SPHERICAL)
     UIJ1(I)=COFA1*DLTXNE(I,1)+COFA2*DLTYNE(I,1)
     VIJ1(I)=COFA5*DLTXNE(I,1)+COFA6*DLTYNE(I,1)
#    else
     XIJ=XIJC(I)-XC(IA)
     YIJ=YIJC(I)-YC(IA)
     UIJ1(I)=COFA1*XIJ+COFA2*YIJ
     VIJ1(I)=COFA5*XIJ+COFA6*YIJ
#    endif
     UALFA_TMP=ABS(UA_FML(IA)-UA_FML(IB_TMP))/ABS(UIJ1(I)+EPSILON(EPS))
     VALFA_TMP=ABS(VA_FML(IA)-VA_FML(IB_TMP))/ABS(VIJ1(I)+EPSILON(EPS))
     IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
     IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
     UALFA(IA)=MIN(UALFA(IA),UALFA_TMP)
     VALFA(IA)=MIN(VALFA(IA),VALFA_TMP)

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)

     A1UIB1 = A1U(IB_TMP,1)
     A1UIB2 = A1U(IB_TMP,2)
     A1UIB3 = A1U(IB_TMP,3)
     A1UIB4 = A1U(IB_TMP,4)
     A2UIB1 = A2U(IB_TMP,1)
     A2UIB2 = A2U(IB_TMP,2)
     A2UIB3 = A2U(IB_TMP,3)
     A2UIB4 = A2U(IB_TMP,4)
        
     COFA3=A1UIB1*UA_FML(IB_TMP)+A1UIB2*UA_FML(K1)+A1UIB3*UA_FML(K2)+A1UIB4*UA_FML(K3)
     COFA4=A2UIB1*UA_FML(IB_TMP)+A2UIB2*UA_FML(K1)+A2UIB3*UA_FML(K2)+A2UIB4*UA_FML(K3)
     COFA7=A1UIB1*VA_FML(IB_TMP)+A1UIB2*VA_FML(K1)+A1UIB3*VA_FML(K2)+A1UIB4*VA_FML(K3)
     COFA8=A2UIB1*VA_FML(IB_TMP)+A2UIB2*VA_FML(K1)+A2UIB3*VA_FML(K2)+A2UIB4*VA_FML(K3)

!     COFA3=A1U(IB,1)*UA(IB)+A1U(IB,2)*UA(K1)+A1U(IB,3)*UA(K2)+A1U(IB,4)*UA(K3)
!     COFA4=A2U(IB,1)*UA(IB)+A2U(IB,2)*UA(K1)+A2U(IB,3)*UA(K2)+A2U(IB,4)*UA(K3)
!     COFA7=A1U(IB,1)*VA(IB)+A1U(IB,2)*VA(K1)+A1U(IB,3)*VA(K2)+A1U(IB,4)*VA(K3)
!     COFA8=A2U(IB,1)*VA(IB)+A2U(IB,2)*VA(K1)+A2U(IB,3)*VA(K2)+A2U(IB,4)*VA(K3)
     
#    if defined (SPHERICAL)
     UIJ2(I)=COFA3*DLTXNE(I,2)+COFA4*DLTYNE(I,2)
     VIJ2(I)=COFA7*DLTXNE(I,2)+COFA8*DLTYNE(I,2)
#    else
     XIJ=XIJC(I)-XC(IB_TMP)
     YIJ=YIJC(I)-YC(IB_TMP)
     UIJ2(I)=COFA3*XIJ+COFA4*YIJ
     VIJ2(I)=COFA7*XIJ+COFA8*YIJ
#    endif
     UALFA_TMP=ABS(UA_FML(IA)-UA_FML(IB_TMP))/ABS(UIJ2(I)+EPSILON(EPS))
     VALFA_TMP=ABS(VA_FML(IA)-VA_FML(IB_TMP))/ABS(VIJ2(I)+EPSILON(EPS))
     IF(UALFA_TMP > 1)UALFA_TMP = 1.0_SP
     IF(VALFA_TMP > 1)VALFA_TMP = 1.0_SP
     UALFA(IB_TMP)=MIN(UALFA(IB_TMP),UALFA_TMP)
     VALFA(IB_TMP)=MIN(VALFA(IB_TMP),VALFA_TMP)

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB_TMP)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)/HPRNU
     ! David moved HPRNU and added HVC
     VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB_TMP)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB_TMP)))/HPRNU

!    SHEAR STRESSES
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX(I)=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
     FYY(I)=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

     ENDIF

   END DO

   DO I=1, NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     IB_TMP = IB

     ELIJ=0.5_SP*(EL_FML(J1)+EL_FML(J2))*(rhomud-rhosus)/rhomud
     DIJ=0.5_SP*(D_FML(J1)+D_FML(J2))

! J. Ge for extra surface water level
     ELIJ =ELIJ-0.5_SP*(EL(J1)+EL(J2))*rhosus/rhomud
! J. Ge for extra surface water level

     UIJ1(I)=UA_FML(IA)+UALFA(IA)*UIJ1(I)
     VIJ1(I)=VA_FML(IA)+VALFA(IA)*VIJ1(I)
     UIJ2(I)=UA_FML(IB_TMP)+UALFA(IB_TMP)*UIJ2(I)
     VIJ2(I)=VA_FML(IB_TMP)+VALFA(IB_TMP)*VIJ2(I)

#    if defined (LIMITED_1)
     IF(UIJ1(I)> MAX(UA_FML(IA),UA_FML(IB_TMP)) .OR. UIJ1(I) < MIN(UA_FML(IA),UA_FML(IB_TMP)) .OR.  &
        UIJ2(I)> MAX(UA_FML(IA),UA_FML(IB_TMP)) .OR. UIJ2(I) < MIN(UA_FML(IA),UA_FML(IB_TMP)))THEN
       UIJ1(I)=UA_FML(IA)
       UIJ2(I)=UA_FML(IB_TMP)
     END IF

     IF(VIJ1(I)> MAX(VA_FML(IA),VA_FML(IB_TMP)) .OR. VIJ1(I) < MIN(VA_FML(IA),VA_FML(IB_TMP)) .OR. &
        VIJ2(I)> MAX(VA_FML(IA),VA_FML(IB_TMP)) .OR. VIJ2(I) < MIN(VA_FML(IA),VA_FML(IB_TMP)))THEN
       VIJ1(I)=VA_FML(IA)
       VIJ2(I)=VA_FML(IB_TMP)
     END IF
#    endif


!    NORMAL VELOCITY
     UIJ=0.5_SP*(UIJ1(I)+UIJ2(I))
     VIJ=0.5_SP*(VIJ1(I)+VIJ2(I))
     UN_TMP=-UIJ*DLTYC(I) + VIJ*DLTXC(I)
          
#    if defined (WET_DRY)
     IF((ISWETCE_FML(IA)*ISWETC_FML(IA) == 1 .OR. ISWETCE_FML(IB)*ISWETC_FML(IB) == 1).AND.((WETTING_DRYING_ON.AND.(ISWETCE(IA)*ISWETC(IA) == 1 .OR. ISWETCE(IB)*ISWETC(IB) == 1)).OR.WETTING_DRYING_ON==.FALSE.))THEN
#    else
     IF(ISWETCE_FML(IA)*ISWETC_FML(IA) == 1 .OR. ISWETCE_FML(IB)*ISWETC_FML(IB) == 1)THEN
#    endif

!    ADD CONVECTIVE AND VISCOUS FLUXES
     XADV=DIJ*UN_TMP*&
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1(I))*0.5_SP
     YADV=DIJ*UN_TMP* &
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2(I)+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1(I))*0.5_SP

!    ACCUMULATE FLUX
#  if !defined (MEAN_FLOW)
     ISBC_TMP = ISBC(I)
     XFLUX(IA)=XFLUX(IA)+(XADV+FXX(I)*EPOR(IA))*(1.0_SP-ISBC_TMP)*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY(I)*EPOR(IA))*(1.0_SP-ISBC_TMP)*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX(I)*EPOR(IB))*(1.0_SP-ISBC_TMP)*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY(I)*EPOR(IB))*(1.0_SP-ISBC_TMP)*IUCP(IB)
#    endif


     END IF


!    ACCUMULATE BAROTROPIC FLUX
!for spherical coordinator and domain across 360^o latitude         
#    if defined (SPHERICAL)
     XTMP  = VX(J2)*TPI-VX(J1)*TPI
     XTMP1 = VX(J2)-VX(J1)
     IF(XTMP1 >  180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP =  360.0_SP*TPI+XTMP
     END IF  

     PSTX_FML(IA)=PSTX_FML(IA)-F_ALFA(IA)*GRAV_E(IA)*D1_FML(IA)*ELIJ*DLTYC(I)
     PSTY_FML(IA)=PSTY_FML(IA)+F_ALFA(IA)*GRAV_E(IA)*D1_FML(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
     PSTX_FML(IB)=PSTX_FML(IB)+F_ALFA(IB)*GRAV_E(IB)*D1_FML(IB)*ELIJ*DLTYC(I)
     PSTY_FML(IB)=PSTY_FML(IB)-F_ALFA(IB)*GRAV_E(IB)*D1_FML(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))

#    else
     PSTX_FML(IA)=PSTX_FML(IA)-GRAV_E(IA)*D1_FML(IA)*ELIJ*DLTYC(I)
     PSTY_FML(IA)=PSTY_FML(IA)+GRAV_E(IA)*D1_FML(IA)*ELIJ*DLTXC(I)
     PSTX_FML(IB)=PSTX_FML(IB)+GRAV_E(IB)*D1_FML(IB)*ELIJ*DLTYC(I)
     PSTY_FML(IB)=PSTY_FML(IB)-GRAV_E(IB)*D1_FML(IB)*ELIJ*DLTXC(I)
#    endif  

   END DO

#  else
   
   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)

     J1=IENODE(I,1)
     J2=IENODE(I,2)

     DIJ=0.5_SP*(D_FML(J1)+D_FML(J2))
     ELIJ=0.5_SP*(EL_FML(J1)+EL_FML(J2))*(rhomud-rhosus)/rhomud

! J. Ge for extra surface water level
     ELIJ =ELIJ-0.5_SP*(EL(J1)+EL(J2))*rhosus/rhomud
! J. Ge for extra surface water level

#    if defined (WET_DRY)
     IF((ISWETCE_FML(IA)*ISWETC_FML(IA) == 1 .OR. ISWETCE_FML(IB)*ISWETC_FML(IB) == 1).AND.((WETTING_DRYING_ON.AND.(ISWETCE(IA)*ISWETC(IA) == 1 .OR. ISWETCE(IB)*ISWETC(IB) == 1)).OR.WETTING_DRYING_ON==.FALSE.))THEN
#    else
     IF(ISWETCE_FML(IA)*ISWETC_FML(IA) == 1 .OR. ISWETCE_FML(IB)*ISWETC_FML(IB) == 1)THEN
#    endif

!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)

     A1UIA1 = A1U(IA,1)
     A1UIA2 = A1U(IA,2)
     A1UIA3 = A1U(IA,3)
     A1UIA4 = A1U(IA,4)
     A2UIA1 = A2U(IA,1)
     A2UIA2 = A2U(IA,2)
     A2UIA3 = A2U(IA,3)
     A2UIA4 = A2U(IA,4)
!---------------------------------------------------------------
     COFA1=A1UIA1*UA_FML(IA)+A1UIA2*UA_FML(K1)+A1UIA3*UA_FML(K2)+A1UIA4*UA_FML(K3)
     COFA2=A2UIA1*UA_FML(IA)+A2UIA2*UA_FML(K1)+A2UIA3*UA_FML(K2)+A2UIA4*UA_FML(K3)
     COFA5=A1UIA1*VA_FML(IA)+A1UIA2*VA_FML(K1)+A1UIA3*VA_FML(K2)+A1UIA4*VA_FML(K3)
     COFA6=A2UIA1*VA_FML(IA)+A2UIA2*VA_FML(K1)+A2UIA3*VA_FML(K2)+A2UIA4*VA_FML(K3)
    
#    if defined (SPHERICAL)
     UIJ1=UA_FML(IA)+COFA1*DLTXNE(I,1)+COFA2*DLTYNE(I,1)
     VIJ1=VA_FML(IA)+COFA5*DLTXNE(I,1)+COFA6*DLTYNE(I,1)
#    else
     XIJ=XIJC(I)-XC(IA)
     YIJ=YIJC(I)-YC(IA)
     UIJ1=UA_FML(IA)+COFA1*XIJ+COFA2*YIJ
     VIJ1=VA_FML(IA)+COFA5*XIJ+COFA6*YIJ
#    endif

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)
     IB_TMP = IB

     A1UIB1 = A1U(IB_TMP,1)
     A1UIB2 = A1U(IB_TMP,2)
     A1UIB3 = A1U(IB_TMP,3)
     A1UIB4 = A1U(IB_TMP,4)
     A2UIB1 = A2U(IB_TMP,1)
     A2UIB2 = A2U(IB_TMP,2)
     A2UIB3 = A2U(IB_TMP,3)
     A2UIB4 = A2U(IB_TMP,4)

     COFA3=A1UIB1*UA_FML(IB_TMP)+A1UIB2*UA_FML(K1)+A1UIB3*UA_FML(K2)+A1UIB4*UA_FML(K3)
     COFA4=A2UIB1*UA_FML(IB_TMP)+A2UIB2*UA_FML(K1)+A2UIB3*UA_FML(K2)+A2UIB4*UA_FML(K3)
     COFA7=A1UIB1*VA_FML(IB_TMP)+A1UIB2*VA_FML(K1)+A1UIB3*VA_FML(K2)+A1UIB4*VA_FML(K3)
     COFA8=A2UIB1*VA_FML(IB_TMP)+A2UIB2*VA_FML(K1)+A2UIB3*VA_FML(K2)+A2UIB4*VA_FML(K3)
     
#    if defined (SPHERICAL)
     UIJ2=UA_FML(IB_TMP)+COFA3*DLTXNE(I,2)+COFA4*DLTYNE(I,2)
     VIJ2=VA_FML(IB_TMP)+COFA7*DLTXNE(I,2)+COFA8*DLTYNE(I,2)
#    else
     XIJ=XIJC(I)-XC(IB_TMP)
     YIJ=YIJC(I)-YC(IB_TMP)
     UIJ2=UA_FML(IB_TMP)+COFA3*XIJ+COFA4*YIJ
     VIJ2=VA_FML(IB_TMP)+COFA7*XIJ+COFA8*YIJ
#    endif


!    NORMAL VELOCITY
     UIJ=0.5_SP*(UIJ1+UIJ2)
     VIJ=0.5_SP*(VIJ1+VIJ2)
     UN_TMP=-UIJ*DLTYC(I) + VIJ*DLTXC(I)

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB_TMP)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)
     ! David moved HPRNU and added HVC
     VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB_TMP)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB_TMP)))/HPRNU

!    SHEAR STRESSES
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
     FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

!    ADD CONVECTIVE AND VISCOUS FLUXES
     XADV=DIJ*UN_TMP*&
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1)*0.5_SP
     YADV=DIJ*UN_TMP* &
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1)*0.5_SP

!    ACCUMULATE FLUX
#  if !defined (MEAN_FLOW)
     ISBC_TMP = ISBC(I)
     XFLUX(IA)=XFLUX(IA)+(XADV+FXX*EPOR(IA))*(1.0_SP-ISBC_TMP)*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY*EPOR(IA))*(1.0_SP-ISBC_TMP)*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX*EPOR(IB))*(1.0_SP-ISBC_TMP)*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY*EPOR(IB))*(1.0_SP-ISBC_TMP)*IUCP(IB)
#   endif


     END IF


!    ACCUMULATE BAROTROPIC FLUX
!for spherical coordinator and domain across 360^o latitude         
#    if defined (SPHERICAL)
     XTMP  = VX(J2)*TPI-VX(J1)*TPI
     XTMP1 = VX(J2)-VX(J1)
     IF(XTMP1 >  180.0_SP)THEN
       XTMP = -360.0_SP*TPI+XTMP
     ELSE IF(XTMP1 < -180.0_SP)THEN
       XTMP =  360.0_SP*TPI+XTMP
     END IF  

     PSTX_FML(IA)=PSTX_FML(IA)-F_ALFA(IA)*GRAV_E(IA)*D1_FML(IA)*ELIJ*DLTYC(I)
     PSTY_FML(IA)=PSTY_FML(IA)+F_ALFA(IA)*GRAV_E(IA)*D1_FML(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))
     PSTX_FML(IB)=PSTX_FML(IB)+F_ALFA(IB)*GRAV_E(IB)*D1_FML(IB)*ELIJ*DLTYC(I)
     PSTY_FML(IB)=PSTY_FML(IB)-F_ALFA(IB)*GRAV_E(IB)*D1_FML(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))
#    else
     PSTX_FML(IA)=PSTX_FML(IA)-GRAV_E(IA)*D1_FML(IA)*ELIJ*DLTYC(I)
     PSTY_FML(IA)=PSTY_FML(IA)+GRAV_E(IA)*D1_FML(IA)*ELIJ*DLTXC(I)
     PSTX_FML(IB)=PSTX_FML(IB)+GRAV_E(IB)*D1_FML(IB)*ELIJ*DLTYC(I)
     PSTY_FML(IB)=PSTY_FML(IB)-GRAV_E(IB)*D1_FML(IB)*ELIJ*DLTXC(I)
#    endif     

   END DO
#  endif

#  if defined (SPHERICAL)
   CALL ADVAVE_EDGE_XY_FML(XFLUX,YFLUX,0.0_SP)
#  endif


   DO I = 1,N
     ISWETTMP = ISWETCE_FML(I)*ISWETC_FML(I)
     XFLUX(I) = XFLUX(I)*ISWETTMP
     YFLUX(I) = YFLUX(I)*ISWETTMP
   END DO
 


!
!-------------------------SET BOUNDARY VALUES----------------------------------!
!

!  MODIFY BOUNDARY FLUX
      DO I=1,N
        IF(ISBCE(I) == 2) THEN
          XFLUX(I)=(XFLUX(I)+Fluxobn(I)*UA_FML(I))*IUCP(I)
          YFLUX(I)=(YFLUX(I)+Fluxobn(I)*VA_FML(I))*IUCP(I)
        ENDIF
      END DO

#  if !defined (LIMITED_NO) 
   DEALLOCATE(UIJ1,VIJ1,UIJ2,VIJ2,STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("Unexpected deallocation error for UIJ1,VIJ1,UIJ2 and VIJ2.")
   DEALLOCATE(UALFA,VALFA,STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("Unexpected deallocation error for UALFA,VALFA.")
   DEALLOCATE(FXX,FYY,STAT=ERROR)
   IF(ERROR /= 0) &
   & CALL FATAL_ERROR("Unexpected deallocation error for FXX,FYY.")
#  endif

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: advave_edge_gcn_fluid_mud"

   END SUBROUTINE ADVAVE_EDGE_GCN_FLUID_MUD
!==============================================================================|




!==============================================================================|
!  Calculate Bottom Drag Coefficient based on Bottom Roughness                 !
!   note:                                                                      !
!   when the log function derived from the constant stress log-viscous         !
!   layer is applied to an estuary, if the value of z0 is close to             !
!   (zz(kbm1)-z(kb)*dt1, drag coefficient "cbc" could become a huge            !
!   number due to near-zero value of alog function. In our application         !
!   we simply cutoff at cbc=0.005. One could adjust this cutoff value          !
!   based on observations or his or her experiences.                           !   
!   CALCULATES:   WUBOT(N), WVBOT(N) : BOTTOM SHEAR STRESSES                   !
!==============================================================================|
SUBROUTINE BOTTOM_ROUGHNESS_FML

  !==============================================================================!
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_WD
  USE MOD_PAR


  IMPLICIT NONE
  INTEGER :: I,II
  REAL(SP), PARAMETER  :: KAPPA = .40_SP   !!VON KARMAN LENGTH SCALE
  REAL(SP), PARAMETER  :: VK2   = .160_SP  !!KAPPA SQUARED
  REAL(SP)             :: ZTEMP,BTPS,RR,U_TAUB,Z0B_GOTM

! USED IN 2D MODEL ONLY
!   REAL(SP), PARAMETER  :: CONST_CD=.0015_SP      !! CD SET CONSTANT TO THIS VALUE
   REAL(SP), PARAMETER  :: ALFA =  .166667_SP, &   !! POWER OF WATER DEPTH
                           NN   = 0.02_SP          !! FACTOR TO DIVIDE
!   REAL(SP), PARAMETER  :: CFMIN   = .0025_SP, &  !! DEEP WATER VALUE
!                           H_BREAK = 1.0_SP,    & !! 
!                           THETA   = 10._SP,   &  !!
!                           LAMB    = 0.3333333333_SP
  !==============================================================================!

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: BOTTOM_ROUGHNESS"

!   1.CONSTANT CD
!    CBC = CONST_CD
!   2.  formula 1
    WHERE (DT1_FML < 3.0_SP)
      CBC_FML = 0.0027_SP
    ELSEWHERE
      CBC_FML = GRAV*(DT1_FML**ALFA/NN)**(-2)
    END WHERE
!   3. formula 2
!    CBC = CFMIN*(1+(H_BREAK/DT1)**THETA)**(LAMB/THETA)

  !==============================================================================|
  !  CALCULATE SHEAR STRESS ON BOTTOM  --> WUBOT/WVBOT                           |
  !==============================================================================|
  DO  I = 1, N
     BTPS= CBC_FML(I)*SQRT(UA_FML(I)**2+VA_FML(I)**2)
     WUBOT_FML(I) = -BTPS * UA_FML(I)
     WVBOT_FML(I) = -BTPS * VA_FML(I)   
  END DO

  !==============================================================================|
  !  Calculate shear stress on nodes (x-component, y-component, magnitude) 
  !==============================================================================|
# if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUBOT_FML,WVBOT_FML)
# endif
  TAUBM_FML = SQRT(WUBOT_FML**2 + WVBOT_FML**2) 

     CALL E2N2D(WUBOT_FML,WUBOT_N_FML)
     CALL E2N2D(WVBOT_FML,WVBOT_N_FML)

  TAUBM_N_FML = SQRT(WUBOT_N_FML**2 + WVBOT_N_FML**2) 
# if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,TAUBM_N_FML)
# endif

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: BOTTOM_ROUGHNESS"

  RETURN
END SUBROUTINE BOTTOM_ROUGHNESS_FML
!==============================================================================|



!==============================================================================|
!  Calculate Bottom and Surface Drag Friction                                  !
                          !   
!   CALCULATES:   WUBOT(N), WVBOT(N) : BOTTOM SHEAR STRESSES                   !
!==============================================================================|

SUBROUTINE FRICTION_FLUID_MUD

  !==============================================================================!
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_WD
  USE MOD_PAR

  IMPLICIT NONE
  INTEGER :: I,II
  REAL(SP), PARAMETER  :: KAPPA = .40_SP   !!VON KARMAN LENGTH SCALE
  REAL(SP), PARAMETER  :: VK2   = .160_SP  !!KAPPA SQUARED
  REAL(dP)             :: ZTEMP,BTPS,RR,U_TAUB,Z0B_GOTM
  real(dp), parameter :: eps2=0.000001
  real(dp) :: eps2t,alfa2,umod,fmrat,tau
! USED IN 2D MODEL ONLY
!   REAL(SP), PARAMETER  :: CONST_CD=.0015_SP      !! CD SET CONSTANT TO THIS VALUE
   REAL(SP), PARAMETER  :: ALFA =  .166667_SP, &   !! POWER OF WATER DEPTH
                           NN   = 0.02_SP          !! FACTOR TO DIVIDE
!   REAL(SP), PARAMETER  :: CFMIN   = .0025_SP, &  !! DEEP WATER VALUE
!                           H_BREAK = 1.0_SP,    & !! 
!                           THETA   = 10._SP,   &  !!
!                           LAMB    = 0.3333333333_SP
  !==============================================================================!

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: FRICTION_FLUID_MUD"

  eps2t = 10.*eps2

!   1.CONSTANT CD
  CBC_FML = fmud
!   2.  formula 1
!    WHERE (DT1 < 3.0_SP)
!      CBC = 0.0027_SP
!    ELSEWHERE
!      CBC = GRAV*(DT1**ALFA/NN)**(-2)
!    END WHERE
!   3. formula 2
!    CBC = CFMIN*(1+(H_BREAK/DT1)**THETA)**(LAMB/THETA)

  !==============================================================================|
  !  CALCULATE SHEAR STRESS ON BOTTOM  --> WUBOT/WVBOT                           |
  !==============================================================================|
  DO  I = 1, N
!     Tm(I)=taubng+fmud*rhomud/8.0*(UA_FML(I)**2+VA_FML(I)**2)
!     IF(UA_FML(I)**2+VA_FML(I)**2<eps2t)THEN
!       alfa2 = taubng/eps2t + fmud*rhomud/8.0
!       BTPS=alfa2*SQRT(UA_FML(I)**2+VA_FML(I)**2)/rhomud
!     ELSE
!       BTPS= Tm(I)/SQRT(UA_FML(I)**2+VA_FML(I)**2)/rhomud
!     END IF

     umod  = sqrt(UA_FML(I)**2+VA_FML(I)**2)
     IF(umod**2<eps2t)THEN
       alfa2 = taubng/eps2t + fmud*rhomud/8.0
       BTPS=alfa2*umod/rhomud
     ELSE
       Tau=taubng+fmud*rhomud/8.0*umod**2
       BTPS= tau/umod/rhomud
     END IF

     WUBOT_FML(I) = -BTPS * UA_FML(I)
     WVBOT_FML(I) = -BTPS * VA_FML(I)

     DELTA_U(I)=U(I,KBM1)-UA_FML(I)
     DELTA_V(I)=V(I,KBM1)-VA_FML(I)

     WUSURF2_FML(I)=DELTA_U(I)*fwat*sqrt(DELTA_U(I)**2+DELTA_V(I)**2)/8.0
     WVSURF2_FML(I)=DELTA_V(I)*fwat*sqrt(DELTA_U(I)**2+DELTA_V(I)**2)/8.0

  END DO

  !==============================================================================|
  !  Calculate shear stress on nodes (x-component, y-component, magnitude) 
  !==============================================================================|
# if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUBOT_FML,WVBOT_FML)
    IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUSURF2_FML,WVSURF2_FML)
# endif
    TAUBM_FML = SQRT(WUBOT_FML**2 + WVBOT_FML**2) 

    CALL E2N2D(WUBOT_FML,WUBOT_N_FML)
    CALL E2N2D(WVBOT_FML,WVBOT_N_FML)

    TAUBM_N_FML = SQRT(WUBOT_N_FML**2 + WVBOT_N_FML**2) 

# if defined (MULTIPROCESSOR)
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,TAUBM_N_FML,WUBOT_N_FML,WVBOT_N_FML)
# endif

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: FRICTION_FLUID_MUD"

  RETURN
END SUBROUTINE FRICTION_FLUID_MUD
!==============================================================================|


!==============================================================================|
! source/sink term for mass balance equation of fluid mud layer
!==============================================================================|
subroutine fluid_mud_source_sink(DTsed,Srho,bed_flag,bed_thck,bed_por,bed_frac,taub,tau_ce,erate)
implicit none
real(sp), intent(in ) :: DTsed,Srho
real(sp), dimension(0:MT),intent(in ) :: bed_flag,bed_thck,bed_por,bed_frac,taub,tau_ce,erate

real(sp) :: erat,poro,frac,flag,dep,dz_bed,eflux_avail

integer  :: i,cnt,jj,cc,k
real(SP) :: fac, fac1, fac2,zd1, zd2,cff1,cff2,cff3
real(SP) :: tmp
real,parameter :: c1=0.25,c2=0.42

real(SP), dimension(0:MT) :: Ub
real(SP), dimension(0:MT) :: Vb
real(SP), dimension(0:MT) :: Ubm
real(SP), dimension(0:MT) :: Vbm
real(SP), dimension(0:MT) :: Zr
real(sp), dimension(1:kbm1) :: Un,Vn



real(sp)       :: alfa
real(sp)       :: buoyan
real(sp)       :: cu
real(sp)       :: cu2
real(sp)       :: cv
real(sp)       :: cv2
real(sp)       :: dewat
real(sp)       :: dmud
real(sp)       :: dmud2
real(sp)       :: dtau
real(sp)       :: eps2t
real(sp)       :: erosi
real(sp)       :: fact
real(sp)       :: h0sus
real(sp)       :: re1
real(sp)       :: re2
real(sp)       :: renold
real(sp)       :: tau
real(sp)       :: taum
real(sp)       :: taux
real(sp)       :: tauy
real(sp)       :: uabs
real(sp)       :: um
real(sp)       :: us
real(sp)       :: ustbe2
real(sp)       :: ustin2
real(sp)       :: usttot
real(sp)       :: usty2
real(sp)       :: uv2
real(sp)       :: uvmag2
real(sp)       :: vm
real(sp)       :: vs
real(sp)       :: tauers
real(sp)       :: excbed
real(sp)       :: entr
real(sp)       :: sett
real(sp)       :: soumud
real(sp), parameter :: eps2=0.000001
real(sp), parameter :: vismud = 5.E-6, rencri = 600., tauly = 0.5

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: fluid_mud_source_sink"

!      fact   = (1./cmud - 1./cbed)
!
!  The definition of fact above takes care of a rising bed level when
!  consolidation occurs. This approach requires updating the bed level
!  later on, which has not been implemented yet, so the thickness of
!  mud layer is calculated wrongly. Since changes in the bed level are
!  insignificant anyhow it is easiest to ignore these. This is done by
!  the folowing definition of fact:(L. Merckelbach)
!

fact  = 1./cmud
eps2t = 10.*eps2

tempx=0.0

do i=1,m
# if defined (WET_DRY)
   if(iswetn(i)==1)then
# endif

     Ub(i) = sum(u(nbve(i,1:ntve(i)),kbm1),dim=1)/float(ntve(i))
     Vb(i) = sum(v(nbve(i,1:ntve(i)),kbm1),dim=1)/float(ntve(i))   

     Ubm(i) = sum(ua_fml(nbve(i,1:ntve(i))),dim=1)/float(ntve(i))
     Vbm(i) = sum(va_fml(nbve(i,1:ntve(i))),dim=1)/float(ntve(i))

     dmud = d_fml(i)
     entr  = 0.0_sp
     dewat = 0.0_sp
     sett  = 0.0_sp
     erosi = 0.0_sp
     soumud= 0.0_sp
     taum  = 0.0_sp
     tau   = 0.0_sp
     excbed= 0.0_sp
     ustbe2= 0.0_sp
     ustin2= 0.0_sp
     usttot= 0.0_sp
     usty2 = 0.0_sp
     uv2   = 0.0_sp
     uvmag2= 0.0_sp
     vm    = 0.0_sp
     vs    = 0.0_sp
     tauers= 0.0_sp
     h0sus = 0.0_sp
     dtau  = 0.0_sp
     buoyan= 0.0_sp

!     Source_Fluid_Mud(i) = 0.0_sp
!     if(ngid(i)==3206)Source_Fluid_Mud(i) = 0.0001_sp
!     cycle

     if(iswetn_fml(i)==1)then
       !
       !           *********************************************
       !           * a. mud layer present at actual grid point *
       !           *********************************************
       !
       !           *********************************************
       !           * a.1 computation of source term for effect *
       !           *     of erosion/dewater                    *
       !           *     (exchange mud layer and bed)          *
       !           *********************************************
       !
       !           source term is defined in zeta-point (+);
       !           so first u- and v-velocities will be transformed to zeta-point
       !           by averaging. (kmax = bottom layer)
       um = Ubm(i)
       vm = Vbm(i)
       !           a.1.1 computation of source term for dewater
       !
       uvmag2 = um**2 + vm**2
       !
       ! Erosie only occurs when there is turbulent flow in the mud-layer (Wang)
       ! Dewater occurs when effective Reynolds number is smaller than 400.
       !

       if (uvmag2<eps2) then
          renold = 10.
       else
          re1    = vismud/sqrt(uvmag2)/max(dmud-min_depth_fml, 0.00001_sp)
          re2    = tauly/2/rhomud/uvmag2
          renold = 1./(re1 + re2)
       endif
       !

       dmud2 = dmud
       if (renold<400.0) then
          dewat = -vdew*cmud*DTsed
          if(abs(dewat)>(dmud-min_depth_fml)*cmud)then
             dewat = -(dmud-min_depth_fml)*cmud
             dmud2  = min_depth_fml
          endif
       else
          dewat = 0.0
       endif

       !
       !           a.1.2. computation of shear stress -taum-
       !
       if (uvmag2<=eps2t) then
       !
       !              for velocities < eps, shear stress -taum- is proportional
       !              to magnitude of velocity.
       !
          alfa = taubng/eps2t + fmud*rhomud/8.0
          taum = alfa*uvmag2
       else
       !
       !              for velocities > eps
       !
          taum = taubng + (fmud*rhomud*uvmag2)/8.0
       endif

       !
       !           a.1.3. compute the source term for erosion
       !
       !           dtau = relative shear stress
       !
       tauers = tau_ce(i)
       erat = erate(i)
       poro = bed_por(i)
       frac = bed_frac(i)
       flag = bed_flag(i)
       dz_bed =bed_thck(i)
       dep = Settling(i)
       dtau = taum - tauers
       !
       if (dtau> 0.0 .and. renold>rencri) then
       !
       !              if shear stress > critical shear stress (between layer-bed)
       !
          !erosi = mers*dtau/tauers
 
          erosi = DTsed*erosion(erat,poro,frac,taum,tauers)
          erosi = erosi*flag
          eflux_avail = frac*srho*(1.0-poro) * dz_bed + dep
          erosi = min(erosi, eflux_avail)
       else
       !
       !              if shear stress < critical shear stress (between layer-bed)
       !
          erosi = 0.0
       endif
       !
       !           a.1.4. compute exchange bed-mud layer
       !                  (erosion+dewater)
       !
       excbed = dewat + erosi
         

       !
       !           *********************************************
       !           * a.2 computation of source term for effect *
       !           *     of settling/entrainment               *
       !           *     (interaction suspension/mud layer)    *
       !           *********************************************
       !
       !           a.2.1. calculate richardson number rib
       !

       us = Ub(i)
       vs = Vb(i)
       !
       uv2 = (us - um)**2 + (vs - vm)**2
       !
       !           a.2.2. calculate entrainment rate entr
       !
       !           fwat = friction coefficient suspension/mud layer
       !
       ustin2 = fwat*uv2/8.0
       tau    = rhosus*ustin2
       usty2  = taubng/rhosus
       ustbe2 = taum/rhosus
       usttot = (ustin2**1.5 + ustbe2**1.5)**(1./3.)
       !h0sus  = sepsus(nm) + real(dps(nm),sp)
       h0sus  = d(i)
       !buoyan = ag*max(h0sus, 0.01_sp)*rhofrac
       buoyan = grav*max(h0sus, 0.01_sp)*(rhomud-rhosus)/rhosus
       if (ustin2>usty2) then
          entr = DTsed*(0.5*(ustin2 - usty2)*sqrt(uv2) + 0.42*(usttot**2 -     &
                 & usty2)*usttot)*cmud/(buoyan + 0.25*uv2)
       else
          !entr = max(mers*(tau - tauers)/tauers, 0.0_sp)
          entr = DTsed*erosion(erat,poro,frac,tau,tauers)
          entr = entr*flag
          eflux_avail = frac*srho*(1.0-poro) * dz_bed
          entr = min(entr, eflux_avail)
       endif
       !
       !           a.2.3. calculate surface shear stress tau
       !
       !
       !           a.2.4. calculate settling
       !
       !           tauset : critical shear stress for settling
       !                   (between suspension/mud layer)
       !
       !            if ( tau.gt. tauset ) then
       !
       !              critical shear stress exceeded --> no settling
       !
       !               sett(nm) = 0.0
       !            else
       !
       !              shear stress < critical shear stress --> settling
       !
       !            wssus is the settling velocity including shear stress
       !            (computed in online sediment-module)
       !sett(nm) = wssus(nm)*rsed(nm, kmax)
       sett = Settling(i)
       !
       !     *                  * (1.- tau/tauset)
       !            endif
       !
       !           a.2.5. calculate source term for combined effect of
       !                  erosion/dewater (exchange) and settling/entrainment
       !
       soumud = (sett - entr + excbed)/cmud
       !
       !           a.2.6. check for drying/ flooding
       !
       !           check whether actual grid point will become dry on
       !           next time step.
       !
       !           approximate thickness of mud layer for next half time
       !           step. (=dnext)
       !
       !           if (soumud(nm)*hdt+dmud .lt. 0.0) then
       !
       if (entr>sett + excbed + (dmud2-min_depth_fml)*cmud) then
       !
       !              >>> drying
       !
       !              a.2.7. recalculate entrainment + source term
       !
       !              dmud = thickness mud layer
       !
          !if(ngid(i)==3247)print*,'xxx1:',soumud,sett,entr,excbed
          entr   = sett + excbed + (dmud2-min_depth_fml)*cmud
          !soumud = (sett - entr + excbed)/cmud
          !soumud = soumud*10.0
          !if(ngid(i)==3247)print*,'xxx2:',soumud,sett,entr,excbed
       endif

     !
     else
       !
       !           ************************************************
       !           * b. no mud layer present at actual grid point *
       !           *    so, only suspension layer                 *
       !           ************************************************
       !
       !           b.1 calculate bed shear stress for suspension
       !               layer
       !
       tau = taub(i)
       tauers = tau_ce(i)
       erat = erate(i)
       poro = bed_por(i)
       frac = bed_frac(i)
       flag = bed_flag(i)
       dz_bed = bed_thck(i)
       dep = Settling(i)
       !
       !           b.3. calculate entrainment rate entr and source term soumud
       !
       uv2=Ub(i)**2 + Vb(i)**2
       ustin2 = tau/rhosus
       usty2  = taubng/rhosus
       usttot = sqrt(ustin2)
       !h0sus  = sepsus(nm) + real(dps(nm),sp)
       h0sus  = d(i)
       !buoyan = ag*max(h0sus, 0.01_sp)*(rhomud - rhosus)/rhosus
       buoyan = grav*max(h0sus, 0.01_sp)*(rhomud-rhosus)/rhosus
       if (ustin2>usty2) then
          entr = DTsed*(0.5*(ustin2 - usty2)*sqrt(uv2) + 0.42*(usttot**2 -     &
                   & usty2)*usttot)*cmud/(buoyan + 0.25*uv2)
       else
          !entr(nm) = max(mers*(tau - tauers)/tauers, 0.0_sp)
          entr = DTsed*erosion(erat,poro,frac,tau,tauers)
          entr = entr*flag
          eflux_avail = frac*srho*(1.0-poro) * dz_bed + dep
          entr = min(entr, eflux_avail)
       endif
       !
       !           ***********************************************
       !           * b.2. calculate erosion + settling           *
       !           ***********************************************
       !
       !           calculate erosion in suspension layer
       !
       erosi = DTsed*erosion(erat,poro,frac,tau,tauers)
       erosi = erosi*flag
       eflux_avail = frac*srho*(1.0-poro) * dz_bed + dep
       erosi = min(erosi, eflux_avail)

       !if(entr > erosi) entr = erosi

       sett = Settling(i)
       excbed = erosi

       if(erosi >= sett)then
          entr = erosi + sett
          erosi = erosi !- sett
          soumud = -erosi/cmud
       else
          entr = erosi
          soumud = sett/cmud  !(sett - erosi)/cmud
       end if 
      ! soumud = (sett - entr + excbed)/cmud
      ! !
      ! !           a.2.6. check for drying/ flooding
      ! !
      ! !           check whether actual grid point will become dry on
      ! !           next time step.
      ! !
      ! !           approximate thickness of mud layer for next half time
      ! !           step. (=dnext)
      ! !
      ! !            if (soumud(nm)*hdt+dmud .lt. 0.0) then
      ! if(ngid(i)==3247)print*,'yyy1:',entr,soumud
      ! if (entr >sett + excbed ) then
      !    !
      !    !              >>> drying
      !    !
      !    !              a.2.7. recalculate entrainment + source term
      !    !
      !    !              dmud = thickness mud layer
      !    !
      !    entr   = sett + excbed
      !    soumud = (sett - entr + excbed)/cmud
      ! endif

     endif

     Entrainment(i) = entr
     BedErosion(i) = erosi
     dewatering(i) = dewat
     Source_Fluid_Mud(i) = soumud/DTsed

# if defined (WET_DRY)
   else
     Entrainment(i) = 0.0_sp
     Settling(i) = 0.0_sp
     BedErosion(i) = 0.0_sp
     dewatering(i) = 0.0_sp
     Source_Fluid_Mud(i) = 0.0_sp
   end if
# endif
   
   !if(ngid(i)==3247)print*,iswetn_fml(i),Entrainment(i),BedErosion(i)
   !if(ngid(i)==3247)print*, dewatering(i), Source_Fluid_Mud(i)&
   !     &,d_fml(i),settling(i)
end do
    

# if defined (MULTIPROCESSOR)
IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Entrainment,Settling,BedErosion)
IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,dewatering,Source_Fluid_Mud)
# endif

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: fluid_mud_source_sink"

return

end subroutine fluid_mud_source_sink

!==============================================================================|

SUBROUTINE FLUX_OBN_FML(K)
  USE ALL_VARS
  USE MOD_OBCS
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: K
  INTEGER              :: I,J,C1,C2
  REAL(SP)             :: FF,FLUX 

  FLUXOBN = 0.0_SP

  DO I = 1, IOBCN

     J = I_OBC_N(I)
     !Compute Boundary Flux From Continuity Flux Defect
     FLUX = -(ELF_FML(J)-ELRK_FML(J))*ART1(J)/(ALPHA_RK(K)*DTE_fml)-XFLUX_OBCN(I)

     !Set Flux In Adjacent Nonlinear BC Element 1 (If Exists)
     IF(NADJC_OBC(I) > 0) THEN
        C1 = ADJC_OBC(I,1) 
        FF = FLUXF_OBC(I,1)
        FLUXOBN(C1) = FLUXOBN(C1) + FF*FLUX 
     END IF

     !Set Flux In Adjacent Nonlinear BC Element 2 (If Exists)
     IF(NADJC_OBC(I) > 1) THEN
        C2 =  ADJC_OBC(I,2) 
        FF = FLUXF_OBC(I,2)
        FLUXOBN(C2) = FLUXOBN(C2) + FF*FLUX 
     END IF

  END DO

  RETURN
END SUBROUTINE FLUX_OBN_FML
!==============================================================================|

SUBROUTINE BCOND_GCN_2_FML
USE ALL_VARS
IMPLICIT NONE
INTEGER :: I

DO I=1,N

!
!--2 SOLID BOUNDARY EDGES------------------------------------------------------|
!
  IF(ISBCE(I) == 3) THEN
    UAF_FML(I)=0.0_SP
    VAF_FML(I)=0.0_SP
  END IF
END DO

RETURN

END SUBROUTINE BCOND_GCN_2_FML


!==========================================================================|
#  if defined (SPHERICAL)
#  if !defined (SEMI_IMPLICIT)
   SUBROUTINE EXTEL_EDGE_XY_FML(K,XFLUX)       
   use MOD_NORTHPOLE
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:MT)
   REAL(SP) :: DIJ,UIJ,VIJ
   INTEGER  :: I,J,K,I1,IA,IB,JJ,J1,J2,II

   REAL(SP) :: UIJ_TMP,VIJ_TMP,EXFLUX_TMP
   REAL(SP) :: DLTXE_TMP,DLTYE_TMP
   REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP

   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
   
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: extel_edge_XY_FML"
!----------INITIALIZE FLUX ARRAY ----------------------------------------------!

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA) = 0.0_SP
     END IF
     IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB) = 0.0_SP
     END IF  
   END DO  

!---------ACCUMULATE FLUX BY LOOPING OVER CONTROL VOLUME HALF EDGES------------!

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     DIJ = D1_FML(I1)

     UIJ = UA_FML(I1)
     VIJ = VA_FML(I1)

     IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
       UIJ_TMP = -VIJ*COS(XC(I1)*DEG2RAD)-UIJ*SIN(XC(I1)*DEG2RAD)
       VIJ_TMP = -VIJ*SIN(XC(I1)*DEG2RAD)+UIJ*COS(XC(I1)*DEG2RAD)
       
       VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
       VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

       VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))

       DLTXE_TMP = VX2_TMP-VX1_TMP
       DLTYE_TMP = VY2_TMP-VY1_TMP
       
       EXFLUX_TMP = DIJ*(-UIJ_TMP*DLTYE_TMP+VIJ_TMP*DLTXE_TMP)
     END IF  
     
     IF(IA == NODE_NORTHPOLE) THEN    
       XFLUX(IA) = XFLUX(IA)-EXFLUX_TMP
     ELSE IF(IB == NODE_NORTHPOLE)THEN
       XFLUX(IB) = XFLUX(IB)+EXFLUX_TMP
     END IF    

   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: extel_edge_XY_FML"

   RETURN
   END SUBROUTINE EXTEL_EDGE_XY_FML
#  endif
!==============================================================================|

!==============================================================================|
#  if !defined (SEMI_IMPLICIT)
   SUBROUTINE EXTUV_EDGE_XY_FML(K)       
   use MOD_NORTHPOLE
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:NT) :: RESX,RESY
   INTEGER  :: I,II

   REAL(SP) :: UARK_TMP,VARK_TMP,UAF_TMP,VAF_TMP,UA_TMP,VA_TMP

   REAL(SP) :: WUSURF2_TMP,WVSURF2_TMP,WUBOT_TMP,WVBOT_TMP
   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
   
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: extuv_edge_XY_FML"

!
!--ACCUMULATE RESIDUALS FOR EXTERNAL MODE EQUATIONS----------------------------|
!
   DO II=1,NP
     I=NP_LST(II)
     UA_TMP = -VA_FML(I)*COS(XC(I)*DEG2RAD)-UA_FML(I)*SIN(XC(I)*DEG2RAD)
     VA_TMP = -VA_FML(I)*SIN(XC(I)*DEG2RAD)+UA_FML(I)*COS(XC(I)*DEG2RAD)

     WUSURF2_TMP = -WVSURF2_FML(I)*COS(XC(I)*DEG2RAD)-WUSURF2_FML(I)*SIN(XC(I)*DEG2RAD)
     WVSURF2_TMP = -WVSURF2_FML(I)*SIN(XC(I)*DEG2RAD)+WUSURF2_FML(I)*COS(XC(I)*DEG2RAD)

     WUBOT_TMP = -WVBOT_FML(I)*COS(XC(I)*DEG2RAD)-WUBOT_FML(I)*SIN(XC(I)*DEG2RAD)
     WVBOT_TMP = -WVBOT_FML(I)*SIN(XC(I)*DEG2RAD)+WUBOT_FML(I)*COS(XC(I)*DEG2RAD)

     RESX(I) = ADX2D_FML(I)+ADVUA_FML(I)+DRX2D_FML(I)+PSTX_FML(I)                &
               -COR(I)*VA_TMP*D1_FML(I)*ART(I)                       &
               -(WUSURF2_TMP+WUBOT_TMP)*ART(I)
     RESY(I) = ADY2D_FML(I)+ADVVA_FML(I)+DRY2D_FML(I)+PSTY_FML(I)                &
               +COR(I)*UA_TMP*D1_FML(I)*ART(I)                       &
               -(WVSURF2_TMP+WVBOT_TMP)*ART(I)

!
!--UPDATE----------------------------------------------------------------------|
!
     UARK_TMP = -VARK_FML(I)*COS(XC(I)*DEG2RAD)-UARK_FML(I)*SIN(XC(I)*DEG2RAD)
     VARK_TMP = -VARK_FML(I)*SIN(XC(I)*DEG2RAD)+UARK_FML(I)*COS(XC(I)*DEG2RAD)

     UAF_TMP = (UARK_TMP*(H1_FML(I)+ELRK1_FML(I))              &
               -ALPHA_RK(K)*DTE_fml*RESX(I)/ART(I))/(H1_FML(I)+ELF1_FML(I))
     VAF_TMP = (VARK_TMP*(H1_FML(I)+ELRK1_FML(I))              &
               -ALPHA_RK(K)*DTE_fml*RESY(I)/ART(I))/(H1_FML(I)+ELF1_FML(I))
		       
     UAF_FML(I)  = VAF_TMP*COS(XC(I)*DEG2RAD)-UAF_TMP*SIN(XC(I)*DEG2RAD)
     VAF_FML(I)  = UAF_TMP*COS(XC(I)*DEG2RAD)+VAF_TMP*SIN(XC(I)*DEG2RAD)
     VAF_FML(I)  = -VAF_FML(I)			    

   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: extuv_edge_XY_FML"

   RETURN
   END SUBROUTINE EXTUV_EDGE_XY_FML
#  endif


   SUBROUTINE ADVAVE_EDGE_XY_FML(XFLUX,YFLUX,IFCETA)
   use MOD_NORTHPOLE
   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)
   REAL(SP) :: IFCETA

   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
 
   CALL FATAL_ERROR("PLease Update the ADVAVE Calculation with North P&
        &ole in Spherical Coordinate")

   END SUBROUTINE ADVAVE_EDGE_XY_FML

#  endif

!==============================================================================|

!==============================================================================|
!   READ IN STATIC WATER DEPTH AND CALCULATE RELATED QUANTITIES                |
!                                                                              |
!   INPUTS: H(NNODE) BATHYMETRIC DEPTH AT NODES				       |
!   INITIALIZES: D(NNODE) DEPTH AT NODES				       |
!   INITIALIZES: DT(NNODE) ???					               |
!   INITIALIZES: H1(NNODE) BATHYMETRIC DEPTH AT ELEMENTS		       |
!   INITIALIZES: D1(NNODE) DEPTH AT NODES                		       | 
!   INITIALIZES: DT1(NNODE) ??                                   	       |
!==============================================================================|

   SUBROUTINE SET_WATER_DEPTH_FML  
!------------------------------------------------------------------------------|
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: TEMP
   INTEGER  :: I,K,J1,J2
!------------------------------------------------------------------------------|
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_WATER_DEPTH_FML" 


!
!  ADJUST STATIC HEIGHT AND CALCULATE DYNAMIC DEPTHS (D) AND (DT)
!
   H_FML  = H_FML + STATIC_SSH_ADJ
   D_FML  = H_FML + EL_FML
   DT_FML = H_FML + ET_FML

!
!  ADJUST SIGMA VALUES ON OUTER BOUNDARY 
!

#  if !defined (PLBC)
   IF(IOBCN > 0) THEN
     DO I = 1,IOBCN
       J1 = I_OBC_N(I)
       J2 = NEXT_OBC(I)
       H_FML(J1) = H_FML(J2)
       D_FML(J1) = D_FML(J2)
       DT_FML(J1) = DT_FML(J2)	 
     END DO
   END IF
#  endif

! 
!  CALCULATE FACE-CENTERED VALUES OF BATHYMETRY AND DEPTH
!
   DO I=1,NT    
     H1_FML(I)  = (H_FML(NV(I,1))+H_FML(NV(I,2))+H_FML(NV(I,3)))/3.0_SP
     D1_FML(I)  = H1_FML(I)+EL1_FML(I)
     DT1_FML(I) = H1_FML(I)+ET1_FML(I)
      
   END DO

# if defined (MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,H_FML,D_FML,DT_FML)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,H1_FML,D1_FML,DT1_FML)
# endif
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_WATER_DEPTH_FML" 
   RETURN
 END SUBROUTINE SET_WATER_DEPTH_FML
!==============================================================================|

!==============================================================================|
   SUBROUTINE SET_WD_DATA_FML
!------------------------------------------------------------------------------|
!  INITIALIZE ARRAYS USED FOR WET/DRY TREATMENT                                |
!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE
   INTEGER :: I

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: SET_WD_DATA_FML"

   IF(STARTUP_TYPE == STARTUP_TYPE_COLDSTART) THEN

!-------- SET WET/DRY FLAGS AND MODIFY WATER SURFACE ELEVATION-----------------!
     CALL WET_JUDGE_FML
!-------- EXCHANGE MODIFIED FREE SURFACE ELEVATION ACROSS PROCESSOR BOUNDS-----!

#  if defined (MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,ELF_FML)
#  endif

!-------- TRANSFER ELEVATION FIELD TO DEPTH AND OLD TIME LEVELS----------------!

     EL1_FML = ELF1_FML
     D1_FML  = H1_FML + EL1_FML
     EL_FML = ELF_FML
     ET_FML = EL_FML
     D_FML  = EL_FML + H_FML
     DT_FML = D_FML
     DTFA_FML = D_FML
     ET1_FML = EL1_FML
     DT1_FML = D1_FML
   END IF 


   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SET_WD_DATA_FML"
   RETURN
   END SUBROUTINE SET_WD_DATA_FML

!==============================================================================|
!==============================================================================|

   SUBROUTINE ALLOC_WD_DATA_FML  

!------------------------------------------------------------------------------|
!  ALLOCATE AND INITIALIZE WET/DRY TREATMENT ARRAYS                            |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE
   INTEGER NCT,NDB

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: ALLOC_WD_DATA_FML"

# if !defined (DOUBLE_PRECISION)
   NDB = 1       !!GWC BASE THIS ON KIND
# else
   NDB = 2
# endif

!-----variables controlling porosities through wet/dry determination----------------!
                                                                                                                          
   ALLOCATE(ISWETN_FML(0:MT))        ; ISWETN_FML  = 1
   ALLOCATE(ISWETC_FML(0:NT))        ; ISWETC_FML  = 1
   ALLOCATE(ISWETNT_FML(0:MT))       ; ISWETNT_FML = 1
   ALLOCATE(ISWETCT_FML(0:NT))       ; ISWETCT_FML = 1
   ALLOCATE(ISWETCE_FML(0:NT))       ; ISWETCE_FML = 1

   ALLOCATE(WETNODES_FML(0:MT))        ; WETNODES_FML  = 1
   ALLOCATE(WETCELLS_FML(0:NT))        ; WETCELLS_FML  = 1
   ALLOCATE(WETNODES_CUR(0:MT))        ; WETNODES_CUR  = 1
   ALLOCATE(WETCELLS_CUR(0:NT))        ; WETCELLS_CUR  = 1

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: ALLOC_WD_DATA_FML"
   RETURN
   END SUBROUTINE ALLOC_WD_DATA_FML

!==============================================================================|
!==============================================================================|

   SUBROUTINE WET_JUDGE_FML

!------------------------------------------------------------------------------|
!  DETERMINE IF NODES/ELEMENTS ARE WET OR DRY                                  |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE
   REAL(SP) :: DTMP
   INTEGER  :: ITA_TEMP
   INTEGER  :: I,IL,IA,IB,K1,K2,K3,K4,K5,K6

   integer :: KT

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: WET_JUDGE_FML"
!
!--Determine If Node Points Are Wet/Dry Based on Depth Threshold---------------!
!
   ISWETN_FML = 1
   DO I = 1, M
     DTMP = H_FML(I) + ELF_FML(I)
     IF((DTMP - MIN_DEPTH_FML) < 1.0E-6_SP) ISWETN_FML(I) = 0
!     if(ngid(i)==3206)then
!       print*,'wd_judge1:',DTMP, H_FML(I) , ELF_FML(I), ISWETN_FML(I)
!     endif
   END DO



!
!--Determine if Cells are Wet/Dry Based on Depth Threshold---------------------!
!
   ISWETC_FML = 1
   DO I = 1, N
     DTMP =  MAX(ELF_FML(NV(I,1)),ELF_FML(NV(I,2)),ELF_FML(NV(I,3)))  + &
             MIN(  H_FML(NV(I,1)),  H_FML(NV(I,2)),  H_FML(NV(I,3)))
     IF((DTMP - MIN_DEPTH_FML) < 1.0E-6_SP) ISWETC_FML(I) = 0
   END DO

!
!--A Secondary Condition for Nodal Dryness-(All Elements Around Node Are Dry)--!
!
   DO I = 1, M
     IF(SUM(ISWETC_FML(NBVE(I,1:NTVE(I)))) == 0)  ISWETN_FML(I) = 0
!     if(ngid(i)==3206)then
!       DTMP = H_FML(I) + ELF_FML(I)
!       print*,'wd_judge2:',DTMP, H_FML(I) , ELF_FML(I), ISWETN_FML(I)
!     endif
   END DO

!
!--Adjust Water Surface So It Does Not Go Below Minimum Depth------------------!
!
   !ELF_FML = MAX(ELF_FML,-H_FML + MIN_DEPTH_FML)
   WHERE(ISWETN_FML>0)
     ELF_FML = MAX(ELF_FML,-H_FML + MIN_DEPTH_FML)
   ELSEWHERE
     ELF_FML = -H_FML + MIN_DEPTH_FML
   END WHERE
!
!--Recompute Element Based Depths----------------------------------------------!
!
   DO I = 1, N
     ELF1_FML(I) = ONE_THIRD*(ELF_FML(NV(I,1))+ELF_FML(NV(I,2))+ELF_FML(NV(I,3)))
   END DO

!
!--Extend Element/Node Based Wet/Dry Flags to Domain Halo----------------------!
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     CALL AEXCHANGE(EC,MYID,NPROCS,ISWETC_FML)
     CALL AEXCHANGE(NC,MYID,NPROCS,ISWETN_FML)
   END IF
#  endif

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: WET_JUDGE_FML"
   RETURN
   END SUBROUTINE WET_JUDGE_FML
!==============================================================================|


!==============================================================================|

   SUBROUTINE WD_UPDATE_FML(INCASE)

!------------------------------------------------------------------------------|
!  SHIFT WET/DRY VARIABLES TO NEW TIME LEVELS                                  |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: INCASE
   INTEGER :: I

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: WD_UPDATE_FML"

   SELECT CASE(INCASE)

!------------------------------------------------------------------------------!
   CASE(1)    !! SHIFT AT END OF EXTERNAL MODE
!------------------------------------------------------------------------------!
   ISWETCE_FML=ISWETC_FML
!------------------------------------------------------------------------------!
   CASE(2)    !! UPDATE NODE WET/DRY AFTER DEPTH ADJUSTMENT
!------------------------------------------------------------------------------!
   DO I = 1,M
     IF(DTFA_FML(I)-MIN_DEPTH_FML <= 1.0E-6_SP) THEN
       ISWETN_FML(I) = 0
     END IF 
#    if defined (WET_DRY)
     IF(WETTING_DRYING_ON.AND.ISWETN(I)==0)ISWETN_FML(I) = 0
#    endif
   END DO
!------------------------------------------------------------------------------!
   CASE(3)    !! SHIFT VARIABLES AT END OF INTERNAL MODE
!------------------------------------------------------------------------------!

   ISWETCT_FML=ISWETC_FML
   ISWETNT_FML=ISWETN_FML

   END SELECT

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: WD_UPDATE_FML"
   RETURN
   END SUBROUTINE WD_UPDATE_FML



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

# endif

End Module Mod_Fluid_Mud
