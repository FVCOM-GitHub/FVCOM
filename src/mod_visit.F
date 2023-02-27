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

!!==================================================================================
!!==================================================================================
!!
!! US FVCOM VISIT SIMULATION INTERFACE 
!! PROGRAMER: David StUebe 
!!
!!==================================================================================
!!==================================================================================


!  MODULE INTERFACE DOES NOT WORK WITH VISIT LIBRARY NAMING SCHEME
!
!==================================================================================
!==================================================================================

# ifdef VISIT

MODULE MOD_VISIT
  USE CONTROL, only : MSR, IPT, MPI_FVCOM_GROUP, Lag_particles_on
  USE MOD_PREC
  USE particle_class
  USE MOD_LAG, only : particle_list
  USE MOD_UTILS
  USE MOD_TIME
  USE MOD_NCDIO, ONLY : VISIT_CMD_DUMP
  IMPLICIT NONE

  ! Variables Used in us_fvcom.F and run_data.F go here

   CHARACTER(LEN=80) :: VISIT_OPT   !!(advanced/basic)  

   type(time) :: visit_time_ext
   type(time) :: visit_time_int
   integer :: visit_cycle

! Variables only used in visitsim.F and mod_visit.F

  LOGICAL :: VisitOneTStep

  LOGICAL :: VisitHalt

  LOGICAL :: VisitAnimate
  
  integer :: VisitRunFlag 
  
  integer :: VisitStep
  INTEGER :: VisitStepCount
    
  integer :: VisitParRank
  
!  integer :: loopcount
  
  INTEGER, parameter :: VISIT_COMMAND_PROCESS = 0
  INTEGER, parameter :: VISIT_COMMAND_SUCCESS = 1
  INTEGER, parameter :: VISIT_COMMAND_FAILURE = 2

  INTEGER, parameter :: VISIT_STEP_EXT = 0
  INTEGER, parameter :: VISIT_STEP_INT = 1
  INTEGER, parameter :: VISIT_STEP_10XINT = 2
  INTEGER, parameter :: VISIT_STEP_100XINT = 3

  
  INTEGER, parameter :: VISIT_TWODMESH = 1
  INTEGER, parameter :: VISIT_BATHYMESH = 2
  INTEGER, parameter :: VISIT_SSHMESH = 3
  INTEGER, parameter :: VISIT_LAYERMESH = 4
  INTEGER, parameter :: VISIT_LEVELMESH = 5
  
  
  integer :: BROADCASTINTCOUNT
  
  integer :: BROADCASTSTRCOUNT
  
  integer :: SLAVECALLBACKCOUNT
  
  TYPE VisitMeshType
     LOGICAL :: MESHALLOCATED
     TYPE(TIME) :: Updated_Time ! the time at which the mesh was last updated  
     INTEGER :: DATAOWNER
     INTEGER :: NDIMS
     INTEGER ::Nodes
     INTEGER :: Zones
     INTEGER :: LCONN 
     INTEGER :: MESH
     INTEGER :: LYRS
     INTEGER :: LVLS
     INTEGER, POINTER, DIMENSION(:,:) :: NV
     INTEGER, POINTER, DIMENSION(:,:) :: GHOSTZONES
     REAL*4, POINTER, DIMENSION(:,:) :: VX,VY,VZ
  END TYPE VisitMeshType
  
  TYPE LAGS
     REAL*4, POINTER, DIMENSION(:) :: S
  END TYPE LAGS

  TYPE VisitLag
     TYPE(TIME) :: Updated_Time ! the time at which the mesh was last updated  
     INTEGER :: NDIMS
     INTEGER ::Nodes
     REAL*4, POINTER, DIMENSION(:) :: VX,VY,VZ,VU,VW,VV,VD
     TYPE(LAGS), POINTER, DIMENSION(:) :: VS
  END TYPE VisitLag

  TYPE VisitSphereVel
     REAL*4, POINTER, DIMENSION(:,:) :: VelX,VelY,VelZ
     TYPE(TIME) :: Updated_Time ! the time at which the velocity was last updated       
  END TYPE VisitSphereVel



  TYPE(VisitMeshType), TARGET, DIMENSION(5) :: VISIT_MESH

  TYPE(VisitLag), TARGET :: VISIT_LAGDATA
  
# if defined (SPHERICAL)
  TYPE(VisitSphereVel), target :: VisitSphericalVel,&
       & VisitSphericalAVel, VisitSphericalWindVel
# if defined (ICE)
  TYPE(VisitSphereVel), target :: VisitSphericalIceVel
# endif

# endif

  SAVE
  
  
  !===================================================================================|
CONTAINS   !!INCLUDED SUBROUTINES FOLLOW
  !===================================================================================|
  !===================================================================================
  


  !============================================================================
  !============================================================================
  !============================================================================
  !
  ! Adding fvcom to visit data handeling functions:
  ! FVCOM data is in 2 dimensional arrays with rows and columns padding
  ! the idicies. These functions make passing the data to visit as
  ! transparent as possible to the user by making the mesh data complex.
  !
  !============================================================================
  !============================================================================
  !============================================================================
  
  SUBROUTINE PRINT_MESH(vmp)
    implicit none
    TYPE(VISITMESHTYPE), POINTER :: vmp
    
    if (.NOT. ASSOCIATED(VMP)) then
       WRITE(IPT,*)"Mesh data pointer passed to UPDATE MESH was not associa&
            &ted"
       return
    end if
    

    WRITE(IPT,*)"==========================================="
    WRITE(IPT,*)"VMP%MESH=",VMP%MESH  
    if(VMP%MESHALLOCATED) then
       WRITE(IPT,*)"VMP has benn allocated" 
    else
       WRITE(IPT,*)"VMP has NOT benn allocated" 
    end if
    WRITE(IPT,*)"VMP%NDIMS=",VMP%NDIMS
    WRITE(IPT,*)"VMP%LYRS=",VMP%LYRS
    WRITE(IPT,*)"VMP%LVLS=",VMP%LVLS
    WRITE(IPT,*)"VMP%Nodes=",VMP%Nodes
    WRITE(IPT,*)"VMP%Zones=",VMP%Zones
    WRITE(IPT,*)"VMP%LCONN=",VMP%LCONN
    WRITE(IPT,*)"VMP%DataOwner=",VMP%DataOwner
    WRITE(IPT,*)"VMP%Updated_Time=",VMP%Updated_Time
    WRITE(IPT,*)"==========================================="
    
    
  END SUBROUTINE PRINT_MESH
  
  SUBROUTINE UPDATE_MESH(vmp,error)
    USE ALL_VARS
    implicit none
    include "visitfortransiminterface.inc"
    integer, intent(out) :: error
    integer :: ind, VnCP, lind, uind
    TYPE(VISITMESHTYPE), POINTER :: vmp
    REAL*4, pointer, Dimension(:) :: PX,PY,PZ

    error =-1 ! BAD RESULT!
    
    if (.NOT. ASSOCIATED(VMP)) then
       if(MSR) WRITE(IPT,*) "Mesh data pointer passed to UPDATE MESH was not associa&
            &ted"
       return
    end if    

    if(dbg_set(dbg_scl)) then
       WRITE(IPT,*)"Starting update_mesh"
       CALL PRINT_MESH(VMP)
    end if


    IF (.NOT.VMP%MESHALLOCATED) then
       
       SELECT CASE(VMP%MESH)
          
       CASE ( VISIT_TWODMESH )
          VnCP=4  ! The number of integers to define one cell
          VMP%NDIMS=2 ! Spatial dimension of the mesh
          VMP%LYRS=1 ! The number of layers
          VMP%LVLS=1 ! The number of levels
          VMP%Nodes=VMP%LVLS * MT ! Total number of nodes
          VMP%Zones=VMP%LYRS * NT ! Total number of zones or cells
          VMP%LCONN=VnCP*VMP%Zones ! The number of integer in the
          ! connectivity array passed to visit
          VMP%DATAOWNER = VISIT_OWNER_SIM
         ! SINCE THE TWOD MESH IS STATIC SET THE TIME TO ZERO AND LEAVE IT
          VMP%Updated_Time=ZEROTIME
       case ( VISIT_BATHYMESH )
          VnCP=4  ! The number of integers to define one cell
          VMP%NDIMS=3
          VMP%LYRS=1
          VMP%LVLS=1
          VMP%Nodes= VMP%LVLS * MT
          VMP%Zones= VMP%LYRS * NT
          VMP%LCONN=VnCP*VMP%Zones
          VMP%DATAOWNER = VISIT_OWNER_SIM
           ! SINCE THE BATHYMETRY MESH IS STATIC SET THE TIME TO ZERO AND LEAVE IT
          VMP%Updated_Time=ZEROTIME
       CASE ( VISIT_SSHMESH )
          VnCP=4  ! The number of integers to define one cell
          VMP%NDIMS=3
          VMP%LYRS=1
          VMP%LVLS=1
          VMP%Nodes= VMP%LVLS * MT
          VMP%Zones= VMP%LYRS * NT
          VMP%LCONN=VnCP*VMP%Zones
          VMP%Nodes=MT
          VMP%Zones=NT
          VMP%LCONN=VnCP*NT
          VMP%DATAOWNER = VISIT_OWNER_SIM
          ! SINCE THE SSH MESH IS NOT STATIC SET THE TIME TO -1
          VMP%Updated_Time=days2time(-1000000)
       CASE ( VISIT_LAYERMESH )

          VnCP=7  ! The number of integers to define one cell
          VMP%NDIMS=3
          VMP%LYRS=KBM2
          VMP%LVLS=KBM1
          VMP%Nodes= VMP%LVLS * MT
          VMP%Zones= VMP%LYRS * NT
          VMP%LCONN=VnCP*VMP%Zones
          VMP%Nodes= KBM1 * MT
          VMP%Zones= KBM2 * NT
          VMP%LCONN=VnCP*NT*KBM2
          VMP%DATAOWNER = VISIT_OWNER_SIM
          ! SINCE THE LAYER MESH IS NOT STATIC SET THE TIME TO -1
          VMP%Updated_Time=days2time(-1000000)
       CASE ( VISIT_LEVELMESH )
          VnCP=7  ! The number of integers to define one cell
          VMP%NDIMS=3 
          VMP%LYRS=KBM1
          VMP%LVLS=KB
          VMP%Nodes= VMP%LVLS * MT
          VMP%Zones= VMP%LYRS * NT
          VMP%LCONN=VnCP*VMP%Zones
          VMP%DATAOWNER = VISIT_OWNER_SIM
          ! SINCE THE LEVEL MESH IS NOT STATIC SET THE TIME TO -1
          VMP%Updated_Time=days2time(-1000000)
          
       CASE DEFAULT
          CALL FATAL_ERROR( "VISIT UPDATE MESH",&
               & "CASE(DEFAULT): bad mesh!")

       END SELECT

       VMP%MESHALLOCATED=.TRUE.
       
       
       ! NV must be allocated by the total number of zones(NT * layers)
       ! by VnCP
       
       ! All other variables will be zones by layer or nodes by layers
       ALLOCATE(VMP%NV(VnCP,VMP%Zones),STAT=error)
       if ( error  .ne. 0 ) then
           if(MSR) WRITE(IPT,*)"Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
               &TE NV, error= ",error,"; MYID= ",MYID
          return
       end if

       if (VnCP == 4) then
          VMP%NV(1,:)=VISIT_CELL_TRI
          VMP%NV(2:4,:)=Transpose(NV(1:NT,1:3)) -1
          
       else if (VnCP == 7) then
          VMP%NV(1,:)=VISIT_CELL_WEDGE
          DO ind =0,(VMP%LYRS-1)

             lind=ind*NT+1
             uind=(ind+1)*NT
             !           1    NT
             !          NT+1 2*NT
             !
             ! (LYRS-1)*NT+1 LYRS*NT
             VMP%NV(2:4,lind:uind)=ind*MT + TRANSPOSE(NV(1:NT,1:3)) -1
             
             VMP%NV(5:7,lind:uind)=(ind+1)*MT + TRANSPOSE(NV(1:NT,1:3)) -1
          END DO
          
       else 
          call FATAL_ERROR("VISIT UPDATE MESH",&
               &"Unknown Cell type: VnCP") 
       end if
       
       ALLOCATE(VMP%GHOSTZONES(NT,VMP%LYRS),STAT=error)
       if ( error  .ne. 0 ) then
          if(MSR)WRITE(IPT,*) "Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
               &TE NV, error= ",error,"; MYID= ",MYID
          return
       end if
       
       DO ind = 1,VMP%LYRS
          VMP%GHOSTZONES(1:N,ind)=0
          if(NT .GT. N) VMP%GHOSTZONES((N+1):NT,ind)=1
       END DO
       
       
       ALLOCATE(VMP%VX(MT,VMP%LVLS),STAT=error)
       if (error  .ne. 0 ) then
          if(MSR)WRITE(IPT,*) "Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
             &TE VX, MYID=",MYID
          return
       end if
       
       
       
       ALLOCATE(VMP%VY(MT,VMP%LVLS),STAT=error)
       if (error  .ne. 0 ) then
          if(MSR)WRITE(IPT,*) "Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
               &TE VY, MYID=",MYID
          return
       end if
       
       
       
       ALLOCATE(VMP%VZ(MT,VMP%LVLS),STAT=error)
       if (error  .ne. 0 ) then
          if(MSR)WRITE(IPT,*) "Allocation error in UPDATE_VMP: CAN NOT ALLOCA&
               &TE VZ, MYID=",MYID
          return
       end if
       
       ! Add static data
       DO ind=1,VMP%LVLS
          VMP%VX(1:MT,ind)=VX(1:MT)
          VMP%VY(1:MT,ind)=VY(1:MT)
       END DO
       
       if (VMP%MESH == VISIT_TWODMESH) VMP%VZ(1:MT,1)=0
       
       if (VMP%MESH == VISIT_BATHYMESH) VMP%VZ(1:MT,1)=-H(1:MT)

# if defined (SPHERICAL)
       if (VMP%MESH == VISIT_TWODMESH .OR. &
            &VMP%MESH == VISIT_BATHYMESH) then

          PX=>VMP%VX(:,1)
          PY=>VMP%VY(:,1)
          PZ=>VMP%VZ(:,1)
          Call Sphere2Cart(PX,PY,PZ,MT)
          

       end if
# endif
      
    END IF  ! 
    
    
    SELECT CASE(VMP%MESH)
       
    CASE ( VISIT_TWODMESH )
       ! Do NOTHING
    CASE ( VISIT_BATHYMESH )
     ! Do NOTHING
    CASE ( VISIT_SSHMESH )
       
       if  (VMP%Updated_Time .NE. VISIT_TIME_EXT) then
          
# if defined (SPHERICAL)
          VMP%VX(1:MT,1)=VX(1:MT)
          VMP%VY(1:MT,1)=VY(1:MT)
          VMP%VZ(1:MT,1)=EL(1:MT)
          
          PX=>VMP%VX(:,1)
          PY=>VMP%VY(:,1)
          PZ=>VMP%VZ(:,1)
          Call Sphere2Cart(PX,PY,PZ,MT)
# else
          VMP%VZ(1:MT,1)=EL(1:MT)
          VMP%Updated_Time=VISIT_TIME_EXT
# endif
!          if (MSR) WRITE(IPT,*) "TIME CHANGED: Updated SSH"
       end if
       
    CASE ( VISIT_LAYERMESH )
       
       if  (VMP%Updated_Time .NE. VISIT_TIME_INT) then
# if defined (SPHERICAL)
          Do ind = 1,VMP%LVLS
             VMP%VX(1:MT,ind)=VX(1:MT)
             VMP%VY(1:MT,ind)=VY(1:MT)
             VMP%VZ(1:MT,ind)=EL(1:MT)+ ZZ(1:MT,ind)*D(1:MT)
             
             PX=>VMP%VX(:,ind)
             PY=>VMP%VY(:,ind)
             PZ=>VMP%VZ(:,ind)
             Call Sphere2Cart(PX,PY,PZ,MT)
          End Do
# else
          Do ind = 1,VMP%LVLS
             VMP%VZ(1:MT,ind)=EL(1:MT)+ ZZ(1:MT,ind)*D(1:MT)
          End Do
# endif       
!          if (MSR) WRITE(IPT,*) "TIME CHANGED: Updated SigmaLayers"
          VMP%Updated_Time=VISIT_TIME_INT
       end if
       
    CASE ( VISIT_LEVELMESH )

       if  (VMP%Updated_Time .NE. VISIT_TIME_INT) then
# if defined (SPHERICAL)
          Do ind = 1,VMP%LVLS
             VMP%VX(1:MT,ind)=VX(1:MT)
             VMP%VY(1:MT,ind)=VY(1:MT)
             VMP%VZ(1:MT,ind)=EL(1:MT)+ Z(1:MT,ind)*D(1:MT)

             PX=>VMP%VX(:,ind)
             PY=>VMP%VY(:,ind)
             PZ=>VMP%VZ(:,ind)
             Call Sphere2Cart(PX,PY,PZ,MT)
          End Do

# else
          Do ind = 1,VMP%LVLS
             VMP%VZ(1:MT,ind)=EL(1:MT)+ Z(1:MT,ind)*D(1:MT)
          End Do
# endif      
!          if (MSR) WRITE(IPT,*) "TIME CHANGED: Updated SigmaLevels"
          VMP%Updated_Time=VISIT_TIME_INT
       end if
       
    CASE DEFAULT
       CALL FATAL_ERROR("CASE(DEFAULT) in UPDATE_MESH",&
            &"ERROR while setting values: bad mesh!")

    END SELECT
    


!    if (MSR)   WRITE(IPT,*) "FINISHED MESH UPDATE"
!    if(MSR)    CALL  PRINT_MESH(VMP)

    
    error = 0
  END SUBROUTINE UPDATE_MESH

  
# if defined (SPHERICAL)
  SUBROUTINE SPHERE2CART(PX,PY,PZ,n)
    IMPLICIT NONE  
    integer, intent(in) :: n
    REAL*4, pointer, Dimension(:), intent(inout) :: PX,PY,PZ
    
    real*4, parameter :: rad=100000.0;
    real*4, parameter :: d2r=3.14159/180.0;
    
    real*4, dimension(n) :: tx, ty, tz
    
    tx = (rad+PZ)*cos(PX*d2r)*cos(PY*d2r);
    ty = (rad+PZ)*sin(PX*d2r)*cos(PY*d2r);
    tz = (rad+PZ)*sin(PY*d2r);
    
    PX=tx
    PY=ty
    PZ=tz
    
  END SUBROUTINE SPHERE2CART
  
  
  SUBROUTINE UpdateSphereVel(VSV)
    USE ALL_VARS, only : NT, KBM1, U, V, WW
    implicit none
    TYPE(VisitSphereVel), POINTER :: VSV
    REAL*4, pointer, Dimension(:) :: PU,PV,PW
    integer :: ind

    if (.NOT.associated(VSV%Velx)) allocate(VSV%Velx(NT,KBM1))
    if (.NOT.associated(VSV%Vely)) allocate(VSV%Vely(NT,KBM1))
    if (.NOT.associated(VSV%Velz)) allocate(VSV%Velz(NT,KBM1))
    
    VSV%Velx=U(1:NT,1:KBM1)
    VSV%Vely=V(1:NT,1:KBM1)
    VSV%Velz=WW(1:NT,1:KBM1)

    Do ind=1,KBM1
       PU=>VSV%VelX(:,ind)
       PV=>VSV%VelY(:,ind)
       PW=>VSV%VelZ(:,ind)
       CALL SphereVel2Cart(PU,PV,PW)
    End Do

    Nullify(PU)
    Nullify(PV)
    Nullify(PW)

    VSV%Updated_Time=Visit_Time_INT
    
  END SUBROUTINE UpdateSphereVel

  SUBROUTINE UpdateSphereAVel(VSV)
    USE ALL_VARS, only : NT, UA, VA
    implicit none
    TYPE(VisitSphereVel), POINTER :: VSV
    REAL*4, pointer, Dimension(:) :: PU,PV,PW
    
    if (.NOT.associated(VSV%Velx)) allocate(VSV%Velx(NT,1))
    if (.NOT.associated(VSV%Vely)) allocate(VSV%Vely(NT,1))
    if (.NOT.associated(VSV%Velz)) allocate(VSV%Velz(NT,1))
    
    VSV%Velx(1:NT,1)=UA(1:NT)
    VSV%Vely(1:NT,1)=VA(1:NT)
    VSV%Velz=0
    
    PU=>VSV%VelX(:,1)
    PV=>VSV%VelY(:,1)
    PW=>VSV%VelZ(:,1)
    CALL SphereVel2Cart(PU,PV,PW)
    
    Nullify(PU)
    Nullify(PV)
    Nullify(PW)
    
    VSV%Updated_Time=Visit_Time_EXT
    
  END SUBROUTINE UpdateSphereAVel

  SUBROUTINE UpdateSphereWindVel(VSV)
    USE ALL_VARS, only : NT, UUWIND, VVWIND
    implicit none
    TYPE(VisitSphereVel), POINTER :: VSV
    REAL*4, pointer, Dimension(:) :: PU,PV,PW
    
    if (.NOT.associated(VSV%Velx)) allocate(VSV%Velx(NT,1))
    if (.NOT.associated(VSV%Vely)) allocate(VSV%Vely(NT,1))
    if (.NOT.associated(VSV%Velz)) allocate(VSV%Velz(NT,1))
    
    VSV%Velx(1:NT,1)=UUWIND(1:NT)
    VSV%Vely(1:NT,1)=VVWIND(1:NT)
    VSV%Velz=0
    
    PU=>VSV%VelX(:,1)
    PV=>VSV%VelY(:,1)
    PW=>VSV%VelZ(:,1)
    CALL SphereVel2Cart(PU,PV,PW)
    
    Nullify(PU)
    Nullify(PV)
    Nullify(PW)
    
    VSV%Updated_Time=Visit_Time_INT
    
  END SUBROUTINE UpdateSphereWindVel
  
#if defined (ICE)
  SUBROUTINE UpdateSphereIceVel(VSV)
    USE mod_ice2d, only: UICE2, VICE2
    USE ALL_VARS, only : NT
    implicit none
    TYPE(VisitSphereVel), POINTER :: VSV
    REAL*4, pointer, Dimension(:) :: PU,PV,PW
    
    if (.NOT.associated(VSV%Velx)) allocate(VSV%Velx(NT,1))
    if (.NOT.associated(VSV%Vely)) allocate(VSV%Vely(NT,1))
    if (.NOT.associated(VSV%Velz)) allocate(VSV%Velz(NT,1))
    
    VSV%Velx(1:NT,1)=UICE2(1:NT)
    VSV%Vely(1:NT,1)=VICE2(1:NT)
    VSV%Velz=0
    
    PU=>VSV%VelX(:,1)
    PV=>VSV%VelY(:,1)
    PW=>VSV%VelZ(:,1)
    CALL SphereVel2Cart(PU,PV,PW)
    
    Nullify(PU)
    Nullify(PV)
    Nullify(PW)
    
    VSV%Updated_Time=Visit_Time_INT
    
  END SUBROUTINE UpdateSphereIceVel
# endif
  
  SUBROUTINE SPHEREVEL2CART(PU,PV,PW)
    USE ALL_VARS, only : XC,YC, NT
    IMPLICIT NONE  
    REAL*4, pointer, Dimension(:), intent(inout) :: PU,PV,PW
    real*4, parameter :: d2r=3.14159/180.0;
    real*4, dimension(NT) :: tu, tv, tw


    tu=-sin(XC(1:NT)*d2r)*PU &
         & -sin(YC(1:NT)*d2r)*cos(XC(1:NT)*d2r)*PV &
         & + cos(XC(1:NT)*d2r)*cos(YC(1:NT)*d2r)*PW

    tv=cos(XC(1:NT)*d2r)*PU &
         & -sin(YC(1:NT)*d2r)*sin(XC(1:NT)*d2r)*PV &
         & + sin(XC(1:NT)*d2r)*cos(YC(1:NT)*d2r)*PW

    tw= 0.0 &
         & + cos(YC(1:NT)*d2r)*PV &
         & + sin(YC(1:NT)*d2r)*PW


    PU = tu
    PV = tv
    PW = tw

  END SUBROUTINE SPHEREVEL2CART

  SUBROUTINE DeAllocate_SphereVel(VSV,error)
    implicit none
    integer, intent(out) :: error
    TYPE(VisitSphereVel), POINTER :: VSV
    error =-1 ! BAD RESULT!
    
    if (associated(VSV%VELX)) deallocate(VSV%VELX)
    if (associated(VSV%VELY)) deallocate(VSV%VELY)
    if (associated(VSV%VELZ)) deallocate(VSV%VELZ)
    

    error = 0
    
  END SUBROUTINE DeAllocate_SphereVel


# endif

  SUBROUTINE DeAllocate_MESH(vmp,error)
    implicit none
    integer, intent(out) :: error
    TYPE(VISITMESHTYPE), POINTER :: vmp
    error =-1 ! BAD RESULT!
    
    if (associated(VMP%VX)) deallocate(VMP%VX)
    if (associated(VMP%VY)) deallocate(VMP%VY)
    if (associated(VMP%VZ)) deallocate(VMP%VZ)
    
    if (associated(VMP%NV)) deallocate(VMP%NV)
    
    if (associated(VMP%GHOSTZONES)) deallocate(VMP%GHOSTZONES)
  
    VMP%MESHALLOCATED=.FALSE.
    error = 0
    
  END SUBROUTINE DeAllocate_MESH
  

! UPDATE FOR LAGRANGIAN DATA
  SUBROUTINE UPDATE_LAG(LMP,error)
    USE LINKED_LIST
    USE MOD_LAG
    implicit none
    integer :: ind, error
    ! FOR LAG TRACKING MESH DATA
    type(link_node), pointer :: lp
    type(VISITLAG), pointer  :: LMP
    
    error = -1

    if  (LMP%Updated_Time .LE. visit_time_int .OR. .NOT. associated(LMP%VX)) then
       
       LMP%NDIMS=3

       if (associated(LMP%VX)) deallocate(LMP%VX)
       if (associated(LMP%VY)) deallocate(LMP%VY)
       if (associated(LMP%VZ)) deallocate(LMP%VZ)

       if (associated(LMP%VU)) deallocate(LMP%VU)
       if (associated(LMP%VV)) deallocate(LMP%VV)
       if (associated(LMP%VW)) deallocate(LMP%VW)

       if (associated(LMP%VS)) deallocate(LMP%VS)
       if (associated(LMP%VD)) deallocate(LMP%VD)

       
       
       LMP%NODES = listsize(particle_list)
       allocate(LMP%VX(LMP%NODES)); LMP%VX=0.0_sp
       allocate(LMP%VY(LMP%NODES)); LMP%VY=0.0_sp
       allocate(LMP%VZ(LMP%NODES)); LMP%VZ=0.0_sp

       allocate(LMP%VU(LMP%NODES)); LMP%VU=0.0_sp
       allocate(LMP%VV(LMP%NODES)); LMP%VV=0.0_sp
       allocate(LMP%VW(LMP%NODES)); LMP%VW=0.0_sp
       
       allocate(LMP%VS(LMP%NODES)) ! IT IS A ARRAY OF POINTERS....
       allocate(LMP%VD(LMP%NODES)); LMP%VD=0.0_sp
       
       call add_scalars

       ind = 1
       lp  => particle_list%first%next
       do
          if(.not. associated(lp) ) exit  !end of list, exit
          


          LMP%VX(ind)=lp%v%x(1)
          LMP%VY(ind)=lp%v%x(2)
          LMP%VZ(ind)=lp%v%zloc

          LMP%VU(ind)=lp%v%U
          LMP%VV(ind)=lp%v%V
          LMP%VW(ind)=lp%v%W

          LMP%VS(ind)%S=>lp%v%S

          LMP%VD(ind)=lp%v%pathlength
          
          ind=ind+1
          lp => lp%next                          !set object
       end do
       

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------

       CALL SPHERE2CART(LMP%VX,LMP%VY,LMP%VZ,LMP%NODES)

# endif

       
       if (MSR) WRITE(IPT,*) "TIME CHANGED: Updated LAGRANGIAN MESH"
       LMP%Updated_Time=visit_time_int
       error = 0

    else
       error = 0
    end if
    
  END SUBROUTINE UPDATE_LAG

  SUBROUTINE DeAllocate_LAG(LMP,error)
    USE LINKED_LIST
    implicit none
    integer, intent(out) :: error
    ! FOR LAG TRACKING MESH DATA
    type(VISITLAG), pointer  :: LMP
    
    error = -1

    if (associated(LMP%VX)) deallocate(LMP%VX)
    if (associated(LMP%VY)) deallocate(LMP%VY)
    if (associated(LMP%VZ)) deallocate(LMP%VZ)
    
    if (associated(LMP%VU)) deallocate(LMP%VU)
    if (associated(LMP%VV)) deallocate(LMP%VV)
    if (associated(LMP%VW)) deallocate(LMP%VW)
    
    if (associated(LMP%VS)) deallocate(LMP%VS)
    if (associated(LMP%VD)) deallocate(LMP%VD)
    
    error = 0
        
  END SUBROUTINE DeAllocate_LAG


END MODULE MOD_VISIT


#else
! if visit is not defined compile a dummy subroutine!

subroutine visit_dummy_mod

implicit none

end subroutine visit_dummy_mod

#endif

!end MODULE MOD_VISIT
