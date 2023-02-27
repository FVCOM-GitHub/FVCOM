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

!  PETSC 2.3.2 / PETSC 2.3.3 COMPATIBILITY ISSUE FOR VecScatter:
! 2.3.2 : CALL VecScatterBegin(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
!
! 2.3.3 : CALL VecScatterBegin(L2G_EL,BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
!


MODULE MOD_PETSc
# if defined (SEMI_IMPLICIT) || defined (NH) || (defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT))
  IMPLICIT NONE

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

# if defined (SEMI_IMPLICIT) || (NH)
  REAL*8, POINTER          :: XVALS_EL(:)
# endif
# if defined (NH)
  REAL*8, POINTER          :: XVALS(:)
# endif
# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  REAL*8, POINTER          :: XVALS_WAVE(:)
# endif

!======== PETSc VARIABLE BLOCK ==============================================
!petsc 2d data structures
# if defined (SEMI_IMPLICIT) || (NH)
  Mat         :: A_EL
  Vec         :: X_EL
  Vec         :: B_EL
  Vec         :: XL_EL
  Vec         :: BL_EL
  PC          :: Pc_EL
  KSP         :: Ksp_EL
  PetscInt    :: PSIZE_EL_HALO
  PetscReal   :: NORM_EL
  IS          :: ISLOCAL_EL
  IS          :: ISGLOBAL_EL
  VecScatter  :: L2G_EL
  VecScatter  :: G2L_EL
  PetscInt    :: ITS_EL
  ISLocalToGlobalMapping :: ISL2G_EL
# endif

# if defined (NH)
!petsc 3d data structure
  Mat         :: A
  Vec         :: X
  Vec         :: B
  Vec         :: U
  Vec         :: BLOCAL
  Vec         :: XLOCAL
  PC          :: Pc
  KSP         :: Ksp
  PetscInt    :: PSIZE
  PetscInt    :: PsizeGL
  PetscInt    :: PSIZE_WHALO
  PetscReal   :: NORM
  IS          :: ISLOCAL 
  IS          :: ISGLOBAL
  VecScatter  :: L2G
  VecScatter  :: G2L
  PetscInt    :: ITS
  ISLocalToGlobalMapping :: ISL2G
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  Mat         :: A_WAVE
  Vec         :: X_WAVE
  Vec         :: B_WAVE
  Vec         :: XL_WAVE
  Vec         :: BL_WAVE
  PC          :: Pc_WAVE
  KSP         :: Ksp_WAVE
  PetscInt    :: PSIZE_WAVE_HALO
  PetscReal   :: NORM_WAVE
  IS          :: ISLOCAL_WAVE
  IS          :: ISGLOBAL_WAVE
  VecScatter  :: L2G_WAVE
  VecScatter  :: G2L_WAVE
  PetscInt    :: ITS_WAVE
  ISLocalToGlobalMapping :: ISL2G_WAVE
# endif

  PetscInt    :: N_VERTS
  PetscInt    :: N_LAY
  PetscInt    :: N_HALO
  PetscScalar :: NEG_ONE = -1.0D0
  PetscScalar :: ONE     =  1.0D0
  PetscScalar :: ZERO    =  0.0D0

  PetscViewer :: viewer
  
!integer mappings
  INTEGER,ALLOCATABLE  :: PLO_2_ALO_NODE(:)  !maps PLO to ALO on a node (surface vertex) basis
  INTEGER,ALLOCATABLE  :: ALO_2_PLO_NODE(:)  !maps ALO to PLO on a node (surface vertex) basis

  INTEGER,ALLOCATABLE,SAVE  :: PLO_2_PGO_NODE(:)

  INTEGER,PARAMETER    :: MAX_NZ_ROW = 30    !max nonzero values in any row

  INTEGER :: PUNIT
  INTEGER :: PSCRN
  LOGICAL :: PMSR  = .TRUE.                  !petsc master processor

!-------------------------------------------------------------------------------------
! The relative residual reduction required for convergence. this can be overriden
! at runtime from the command line.
!-------------------------------------------------------------------------------------
  PetscReal :: RTOL = 1.0E-10

!-------------------------------------------------------------------------------------
! AR Aspect ratio. In the fake poisson discretization, the d^2/dz^2 terms are set to
! this number instead of order 1.
!-------------------------------------------------------------------------------------
  PetscInt  :: AR = 1

!-------------------------------------------------------------------------------------
! CHECK_EXACT If this is  true, we set up a RHS using a solution vector u = 1 (Au=b)
! then we solve Ax=b and make sure x and u are the same (within tol).
!-------------------------------------------------------------------------------------
  LOGICAL   :: CHECK_EXACT = .FALSE.

!-------------------------------------------------------------------------------------
! If this is true. we set the x vector from qnh and this is used by petsc as an
! initial condition for the Solve.  If this is not true we skip setting the initial
! condition (to save time) and Petsc zeroes out the x vector on its own.
!-------------------------------------------------------------------------------------
  LOGICAL   :: USE_LAST = .TRUE.

!-------------------------------------------------------------------------------------
! If this is true, extra debug statements are printed to file unit numbers myid+800
! where myid is the processor id.
!-------------------------------------------------------------------------------------
  LOGICAL   :: PDEBUG   = .FALSE.

!============ END PETSc BLOCK ========================================================

  !----------------------------------------------------------------------------
  CONTAINS
  ! PETSc_SET                :
  ! PETSc_AllocA             :
  ! PETSc_Solve              :
  ! PETSc_Solve_EL           :
  ! PETSc_Solve_WAVE         :
  ! PETSc_CleanUp            : 
  ! And two functions
  ! PLO_2_ALO                : map [i] PLO to [i,k] ALO
  ! ALO_2_PLO                : map [i,k] ALO to [i] PLO
  !----------------------------------------------------------------------------

  SUBROUTINE PETSc_SET
  !============================================================================
  !  ALO:  Application Local Node Ordering
  !  PLO:  Petsc Local Node Ordering
  !  AGO:  Application Global Node Ordering
  !  PGO:  Petsc Global Node Ordering
  !============================================================================
  USE MOD_PAR,  ONLY : EL_PID,NLID,NLID_X
  USE CONTROL,  ONLY : MPI_FVCOM_GROUP
  USE ALL_VARS, ONLY : NVG, MSR, NGL, MGL, MT, NPROCS, MYID, KBM1, M, IPT
  IMPLICIT NONE

  PetscInt:: IERR
  INTEGER :: I,L,J,K,ICNT,PID,N1,PGL_NODE
  INTEGER :: MINDOMSIZE, MINDOM
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: CNT
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: ACCUM
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PARTID
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PART_NODES
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: OWNER_COUNT
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PLO_2_AGO_NODE
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: AGO_2_PGO_NODE
!  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PLO_2_PGO_NODE
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: OWNERS
  CHARACTER(LEN=20)                  :: SUBNAME = 'PETSc_SETUP'
  LOGICAL                            :: NEWOWN

  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PTMP
  
# if defined (NH)
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: PLO_2_PGO
# endif
  !-----------------------------------------------------------------------------
  !setup Petsc vectors and matrices (first call only)
  !we will assume that the non-zero structure of matrix DOES NOT CHANGE
  !during iteration process.  This will allow us to preallocate properly
  !-----------------------------------------------------------------------------

  PSCRN = IPT             ! set petsc unit number for dumping regular information
  PUNIT = 800+MYID        ! set petsc unit number for dumping debug information
  PMSR = MSR              ! set petsc master (primarily for printing) to FVCOM master

  IF(PDEBUG) WRITE(PUNIT,*) 'STARTING: ', TRIM(SUBNAME)
  IF(PMSR)   WRITE(PSCRN,*) 'PETSc SETUP        : BEGINNING'

  CALL PETScInitialize(PETSC_NULL_CHARACTER,IERR) ; CHKERRQ(IERR)
  IF(PMSR)   WRITE(PSCRN,*) 'PETSc INITIALIZE   : COMPLETE'

  !---------------------------------------------------------
  ! set up AO<-->PO mappings, sizes, etc.
  !---------------------------------------------------------

  !-----> simple partitioning does not give great load balance partition the domain by
  ! nodes such that each local node in Petsc Partition is also local in Application partition
  ! But, we do not have overlapping domains What were originally boundary nodes in FVCOM decomp
  ! and had multiple partition owners are assigned to the adjacent partition of lower number
  ! petsc partition id for each node in global application order ==> [partid(1:mgl)]
  ! we should really consider a load balancing procedure where we first assign node that
  ! arent boundary nodes and then assign boundary nodes to the adjacent partition which is smaller.
  ! Ultimately, we want N_verts to be as equal as possible in each domain.  For right now,
  ! use this simple partitioning strategy

!  allocate(partid(mgl)) ; partid = nprocs+1
!  do i=1,ngl
!  do l=1,3
!    partid(nvg(i,l)) = min( partid(nvg(i,l)) , el_pid(i))
!  end do
!  end do

  !-----> more complex partitioning for better load balance partition the domain by nodes
  ! such that each local node in Petsc Partition is also local in Application partition
  ! we will do two sweeps. In a first sweep, any node that was not an interprocessor boundary
  ! node will be assigned to the owner that had the node in the FVCOM partition.
  ! then we will count nodes in each Petsc Partition then we will loop over boundary nodes
  ! (which have multiple owners in the FVCOM partition) and decide who gets it based
  ! on who has the lowest current node count this should result in a nearly optimal loadbalance

  !--sweep 1, determine owners, owner count, assume max # owners = 10
  ALLOCATE(OWNER_COUNT(MGL)) ; OWNER_COUNT = 0
  ALLOCATE(OWNERS(MGL,0:10)) ; OWNERS      = 0
  DO I =1, NGL
    DO L =1, 3
      !check if owner is already in list --> can use loc!
      NEWOWN = .TRUE.
      DO J=1,OWNER_COUNT(NVG(I,L))
        IF(OWNERS(NVG(I,L),J) == EL_PID(I) ) NEWOWN = .FALSE.
      ENDDO
      IF(NEWOWN)THEN
        OWNER_COUNT(NVG(I,L)) = OWNER_COUNT(NVG(I,L)) + 1
        IF(OWNER_COUNT(NVG(I,L)) > 10) THEN
          WRITE(PSCRN,*) 'ERROR IN PETSc PARTITIONING'
          WRITE(PSCRN,*) 'INCREASE OWNER COUNT FROM 10'
          CALL MPI_FINALIZE(IERR)
          STOP
        ENDIF
        OWNERS( NVG(I,L),OWNER_COUNT(NVG(I,L)) ) = EL_PID(I)
      ENDIF
    ENDDO
  ENDDO

  !--now sweep and petsc partition the non-interprocessor boundary nodes
  ALLOCATE(PART_NODES(0:NPROCS)) ; PART_NODES = 0
  ALLOCATE(PARTID(MGL)) ; PARTID = -1
  DO I= 1, MGL
    IF(OWNER_COUNT(I) == 1) THEN
      PART_NODES(OWNERS(I,1)) = PART_NODES(OWNERS(I,1)) + 1
      PARTID(I) = OWNERS(I,1)
    ENDIF
  ENDDO

  !--now sweep boundary nodes, assigning a node to owner with least count
  DO I= 1, MGL
    IF(OWNER_COUNT(I) > 1) THEN
      MINDOMSIZE = HUGE(I)
      DO J= 1, OWNER_COUNT(I)
        !find which adjacent dom has min count
        IF (PART_NODES(OWNERS(I,J)) < MINDOMSIZE) THEN
          MINDOMSIZE = PART_NODES(OWNERS(I,J) )
          MINDOM     = OWNERS(I,J)
        ENDIF
      ENDDO
      !assign boundary node to that partition
      PARTID(I) = MINDOM
      PART_NODES(MINDOM) = PART_NODES(MINDOM) + 1
    ENDIF
  ENDDO

  !--make sure everybody got assigned
  DO I = 1, MGL
    IF(PARTID(I) == -1) THEN
      WRITE(PSCRN,*) 'ERROR IN PETSc PARTITIONING'
      WRITE(PSCRN,*) 'GLOBAL NODE: ',I,' NEVER GOT ASSIGNED A PETSc PARTITION!'
      WRITE(PSCRN,*) 'FIX IT!'
      CALL MPI_FINALIZE(IERR)
      STOP
    ENDIF
  ENDDO

  !--deallocate the data we used for 'fancy' partition
  DEALLOCATE(OWNERS)
  DEALLOCATE(OWNER_COUNT)

  ! lets count (or recount) nodes in each Petsc partition ==> [part_nodes(0:nprocs)]
  ! we need to start from zero so we can access partition (myid-1) later the total will equal MGL,
  ! unlike the sum of nodes over FVCOM partitions which is > MGL
  IF(.NOT.ALLOCATED(PART_NODES) ) ALLOCATE(PART_NODES(0:NPROCS))
  PART_NODES = 0
  DO I= 1, MGL
    PART_NODES( PARTID(I)) = PART_NODES( PARTID(I) ) + 1
  ENDDO

  !count number of nodes in my local Petsc partition ==> [N_verts] halt if partition is empty
  N_VERTS = PART_NODES(MYID)

  print*,'N_VERTS and M=',N_VERTS,M,MT,MYID
  IF(N_VERTS == 0) THEN
    WRITE(PUNIT,*) 'MY PETSc PARTITION DOSE NOT CONTAIN ANY UNKNOWNS'
    CALL MPI_FINALIZE(IERR)
    STOP
  ENDIF
  IF(PDEBUG) WRITE(PUNIT,*) 'NODES IN PETSc PARTION: ',N_VERTS,' APPLICATION',M

# if defined (NH)
  !set number of layers and total problem size ==> [N_lay,Psize]
  N_LAY   = KBM1
  PSIZE   = N_LAY*N_VERTS
  PSIZEGL = MGL*N_LAY
# endif

  !report all sizes(master only)
  IF(PMSR) WRITE(PSCRN,*)'PARTITIONING       : COMPLETE'
  IF(PMSR) THEN
    DO I= 1, NPROCS
      WRITE(PSCRN,*)'PETSc PARTITION: ',I,'  SIZE  ',PART_NODES(I)
    ENDDO
  ENDIF

  !make sure total sum of nodes in all petsc partitions is equal to mgl
  IF(SUM(PART_NODES) /= MGL) THEN
    WRITE(PSCRN,*)'TOTAL SUM OF NODES IN PETSc PARTIONS DOES NOT EQUAL GLOBAL NUMBER OF NODES'
    WRITE(PSCRN,*)'EXITING!'
    CALL MPI_FINALIZE(IERR)
    STOP
  ENDIF

  !determine which global nodes (in application order) are in my Petsc partition
  !gives us a mapping (for interior nodes) of PLO to AGO
  !allocate it to size MGL because we do not know size of petsc partition halo
  !gives us ==> [plo_2_ago_node(1:N_verts)]
  ALLOCATE(PLO_2_AGO_NODE(MGL)) ; PLO_2_AGO_NODE = 0
  ICNT = 0
  DO I= 1, MGL
    IF(PARTID(I) == MYID) THEN
      ICNT = ICNT + 1
      PLO_2_AGO_NODE(ICNT) = I
    ENDIF
  ENDDO

  IF(PDEBUG)THEN
    WRITE(PUNIT,*) 'PLO_2_AGO_NODE'
    DO I= 1, N_VERTS
      WRITE(PUNIT,*) I, PLO_2_AGO_NODE(I)
    ENDDO
  ENDIF

  !get a mapping from petsc local nodes to fvcom local nodes ==> plo_2_alo_node(1:N_verts)
  !we already have the global application order number of each node, we will
  !use nlid to convert this to local application order of each node
  !this will be used to in our plo_2_alo function which is ultimately used to
  !load and unload FVCOM vectors (q and source term)
  !gives us ==> [plo_2_alo_node(N_verts)]
  ALLOCATE(PLO_2_ALO_NODE(N_VERTS)); PLO_2_ALO_NODE = 0
  IF(NPROCS > 1) THEN
    DO I=1,N_VERTS
      PLO_2_ALO_NODE(I) = NLID(PLO_2_AGO_NODE(I))
    ENDDO
  ELSE
    PLO_2_ALO_NODE = PLO_2_AGO_NODE
  ENDIF
  IF(PDEBUG)THEN
    WRITE(PUNIT,*) 'PLO_2_ALO_NODE:'
    DO I= 1, N_VERTS
      WRITE(PUNIT,*) I, PLO_2_ALO_NODE(I)
    ENDDO
  ENDIF

  !get a mapping from application local nodes to petsc local nodes
  !first we map without setting petsc halos
  !gives us ==> [alo_2_plo_node(1:m)]
  ALLOCATE(ALO_2_PLO_NODE(MT)) ; ALO_2_PLO_NODE = 0
  IF(NPROCS > 1)THEN
    ICNT = 0
    DO I= 1, MGL
      IF(NLID_X(I) /= 0 .AND. PARTID(I) == MYID) THEN
        ICNT = ICNT + 1
        ALO_2_PLO_NODE(NLID_X(I)) = ICNT
      ENDIF
    ENDDO
  ELSE
    ALO_2_PLO_NODE = PLO_2_ALO_NODE
  ENDIF
  IF(PDEBUG)THEN
    WRITE(PUNIT,*) 'ALO_2_PLO_NODE:'
    DO I= 1, N_VERTS
      WRITE(PUNIT,*) I, ALO_2_PLO_NODE(I)
    ENDDO
  ENDIF

  IF(ICNT /= N_VERTS)THEN
    WRITE(PUNIT,*) 'PROBLEM IN ALO to PLO MAPPING!'
    CALL MPI_FINALIZE(IERR)
    STOP
  ENDIF

  !we continue mapping ALO to PLO but now we accumulate halo nodes in PLO
  !we also extend our plo_2_ago map to include a mapping from petsc local
  !order to application global order in the petsc partition halos
  !gives us ==> [n_halo]
  !gives us ==> [alo_2_plo_node(m+1:mt)]
  !gives us ==> [plo_2_alo_node(n_verts+1:n_verts+n_halo)]
  N_HALO = 0
  IF(NPROCS > 1)THEN
    DO I= 1, MGL
      IF(NLID_X(I) /= 0 .AND. PARTID(I) /= MYID) THEN
        ICNT = ICNT + 1
        N_HALO = N_HALO + 1
        ALO_2_PLO_NODE(NLID_X(I)) = ICNT
        PLO_2_AGO_NODE(ICNT     ) = I
      ENDIF
    ENDDO
  ENDIF
  IF(PDEBUG) THEN
    WRITE(PUNIT,*) 'ALO_2_PLO_NODE WITH HALO:'
    DO I= 1, MT
      WRITE(PUNIT,*) I, ALO_2_PLO_NODE(I)
    ENDDO
  ENDIF

  !we now have a map from all all local application order nodes, including
  !halos to petsc local nodes, including  halos.  We also have a map from all local
  !petsc partition nodes, including the halo to application global order.
  !what we need is a map from petsc partition nodes including halo to petsc
  !global order, so that we can insert matrix values in local order and they will
  !get mapped accordingly.
  !we will construct a map of application global order to petsc global order and
  !use that to map our palo_2_ago into what we really need, a plo_2_pgo
  !gives us ==> [ago_2_pgo_node(1:MGL)]
  ALLOCATE(AGO_2_PGO_NODE(MGL)) ; AGO_2_PGO_NODE = 0
  ALLOCATE(CNT(NPROCS))         ; CNT   = 0
  ALLOCATE(ACCUM(NPROCS) )      ; ACCUM = 0
  DO I= 2, NPROCS
    ACCUM(I)  = ACCUM(I-1) + PART_NODES(I-1)
  ENDDO
  DO I= 1, MGL
    PID = PARTID(I)
    CNT(PID) = CNT(PID) + 1
    AGO_2_PGO_NODE(I) = ACCUM(PID) + CNT(PID)
  ENDDO
  IF(PDEBUG)THEN
    WRITE(PUNIT,*)'AGO_2_PGO_NODE:'
    DO I= 1, MGL
      WRITE(PUNIT,*) I, AGO_2_PGO_NODE(I)
    ENDDO
  ENDIF

  !now, use our ago_2_pgo_node to construct a plo_2_pgo_node
  !gives us ==> [plo_2_pgo_node(1:n_verts+n_halo)]
  ALLOCATE(PLO_2_PGO_NODE(N_VERTS+N_HALO)) ; PLO_2_PGO_NODE = 0
  DO I= 1, N_VERTS + N_HALO
    PLO_2_PGO_NODE(I) = AGO_2_PGO_NODE( PLO_2_AGO_NODE(I) )
  ENDDO
  IF(PDEBUG)THEN
    WRITE(PUNIT,*) 'COUNT PLO_2_PGO_NODE:'
    DO I=1, N_VERTS+N_HALO
      WRITE(PUNIT,*) I, PLO_2_PGO_NODE(I)
    ENDDO
  ENDIF

  !set a total problem size including petsc partition halo nodes
  ! ==> [Psize_whalo]
# if defined (SEMI_IMPLICIT) || (NH)
  PSIZE_EL_HALO = N_VERTS + N_HALO
# endif

# if defined (NH)
  PSIZE_WHALO   = (N_VERTS + N_HALO)*N_LAY

  !now create a mapping between petsc local order and petsc global
  !order for each unknown in the 3-D mesh
  !we will always assume k (layers) is fastest running index
  ! ==> [plo_2_pgo(1:Psize_whalo)] 
  ! should be 1 to 1 for single processor case
  ALLOCATE(PLO_2_PGO( PSIZE_WHALO )) ; PLO_2_PGO = 0 
  ICNT = 0
  DO I= 1, N_VERTS + N_HALO
    PGL_NODE = PLO_2_PGO_NODE(I)
    DO K= 1, N_LAY
      ICNT = ICNT + 1
      PLO_2_PGO(ICNT) = K + (PGL_NODE-1)*N_LAY
    ENDDO
  ENDDO
  IF(PDEBUG)THEN
    WRITE(PUNIT,*) 'COUNT  PLO_2_PGO'
    DO I= 1, PSIZE_WHALO
      WRITE(PUNIT,*) I, PLO_2_PGO(I)
    ENDDO
  ENDIF
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  PSIZE_WAVE_HALO = N_VERTS + N_HALO
# endif

  IF(PMSR) WRITE(PSCRN,*) 'MAPPINGS           : COMPLETE!'

  !---------------------------------------------------------
  ! set up the petsc matrix A
  !---------------------------------------------------------
  IF(PDEBUG) WRITE(PUNIT,*) 'SETTING UP THE MATRIX'

# if defined (NH)
  !setup petsc local to global mapping (isl2g) to all for local indexing in assembly of matrix
  ALLOCATE(PTMP(PSIZE_WHALO)) ; PTMP = 0
  ICNT = 0
  DO I= 1, N_VERTS + N_HALO
    DO K=1, N_LAY
      ICNT = ICNT + 1
      PTMP(ICNT) = PLO_2_PGO(ICNT)-1
    ENDDO
  ENDDO
  CALL ISLocalToGlobalMappingCreate(PETSC_COMM_SELF,PSIZE_WHALO,PTMP,ISL2G,IERR);CHKERRQ(IERR)
  DEALLOCATE(PTMP)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  ALLOCATE(PTMP(PSIZE_EL_HALO)) ; PTMP = 0
  DO I= 1, N_VERTS + N_HALO
    PTMP(I) = PLO_2_PGO_NODE(I)-1
  ENDDO
  CALL ISLocalToGlobalMappingCreate(PETSC_COMM_SELF,PSIZE_EL_HALO,PTMP,ISL2G_EL,IERR);CHKERRQ(IERR)
  DEALLOCATE(PTMP)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  !setup petsc local to global mapping (isl2g) to all for local indexing in assembly of matrix
  ALLOCATE(PTMP(PSIZE_WAVE_HALO)) ; PTMP = 0
  DO I= 1, N_VERTS + N_HALO
    PTMP(I) = PLO_2_PGO_NODE(I)-1
  END DO
  CALL ISLocalToGlobalMappingCreate(PETSC_COMM_SELF,PSIZE_WAVE_HALO,PTMP,ISL2G_WAVE,IERR);CHKERRQ(IERR)
  DEALLOCATE(PTMP)
# endif

  !preallocate data
  CALL PETSc_AllocA

  !allow user to set options from command line
# if defined (NH)
  CALL MatSetFromOptions(A,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL MatSetFromOptions(A_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL MatSetFromOptions(A_WAVE,IERR);CHKERRQ(IERR)
# endif

  !allowing local indexing for assembly of matrix
# if defined (NH)
  CALL MatSetLocalToGlobalMapping(A,ISL2G,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL MatSetLocalToGlobalMapping(A_EL,ISL2G_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL MatSetLocalToGlobalMapping(A_WAVE,ISL2G_WAVE,IERR);CHKERRQ(IERR)
# endif

  IF(PDEBUG) WRITE(PUNIT,*) 'MatSetLocalToGlobalMapping COMPLETE!'
  IF(PMSR)   WRITE(PSCRN,*) 'MATRIX PREALLOC    : COMPLETE!'

  !---------------------------------------------------------
  ! set up vector index sets used for scatter ops
  !---------------------------------------------------------

  !setup local vector index sets  (Petsc Ordering) -->  these are zero based
# if defined (NH)
  CALL ISCreateStride(PETSC_COMM_SELF,PSIZE,0,1,ISLOCAL,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL ISCreateStride(PETSC_COMM_SELF,N_VERTS,0,1,ISLOCAL_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL ISCreateStride(PETSC_COMM_SELF,N_VERTS,0,1,ISLOCAL_WAVE,IERR);CHKERRQ(IERR)
# endif

  !setup global vector index set (Petsc Ordering) --> zero based
# if defined (NH)
  N1 = ACCUM(MYID)*N_LAY
  CALL ISCreateStride(PETSC_COMM_SELF,PSIZE,N1,1,ISGLOBAL,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL ISCreateStride(PETSC_COMM_SELF,N_VERTS,ACCUM(MYID),1,ISGLOBAL_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL ISCreateStride(PETSC_COMM_SELF,N_VERTS,ACCUM(MYID),1,ISGLOBAL_WAVE,IERR);CHKERRQ(IERR)
# endif

  IF(PMSR) WRITE(PSCRN,*)'VECTOR INDEX SETS  : COMPLETE!'

  !---------------------------------------------------------
  ! set up vector contexts: local vecs (xlocal,blocal)
  !---------------------------------------------------------
# if defined (NH)
  CALL VecCreateSeq(PETSC_COMM_SELF,PSIZE,BLOCAL,IERR);CHKERRQ(IERR)
  CALL VecCreateSeq(PETSC_COMM_SELF,PSIZE,XLOCAL,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL VecCreateSeq(PETSC_COMM_SELF,N_VERTS,BL_EL,IERR);CHKERRQ(IERR)
  CALL VecCreateSeq(PETSC_COMM_SELF,N_VERTS,XL_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL VecCreateSeq(PETSC_COMM_SELF,N_VERTS,BL_WAVE,IERR);CHKERRQ(IERR)
  CALL VecCreateSeq(PETSC_COMM_SELF,N_VERTS,XL_WAVE,IERR);CHKERRQ(IERR)
# endif

  !---------------------------------------------------------
  ! set up vector contexts: global parallel vecs (x,u,b)
  !---------------------------------------------------------
# if defined (NH)
  CALL VecCreate(MPI_FVCOM_GROUP,X,IERR);CHKERRQ(IERR)
  CALL VecSetSizes(X,PSIZE,PETSC_DECIDE,IERR);CHKERRQ(IERR)
  CALL VecSetFromOptions(X,IERR);CHKERRQ(IERR)
  CALL VecDuplicate(X,B,IERR);CHKERRQ(IERR)
  CALL VecDuplicate(X,U,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL VecCreate(MPI_FVCOM_GROUP,X_EL,IERR);CHKERRQ(IERR)
  CALL VecSetSizes(X_EL,N_VERTS,PETSC_DECIDE,IERR);CHKERRQ(IERR)
  CALL VecSetFromOptions(X_EL,IERR);CHKERRQ(IERR)
  CALL VecDuplicate(X_EL,B_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL VecCreate(MPI_FVCOM_GROUP,X_WAVE,IERR);CHKERRQ(IERR)
  CALL VecSetSizes(X_WAVE,N_VERTS,PETSC_DECIDE,IERR);CHKERRQ(IERR)
  CALL VecSetFromOptions(X_WAVE,IERR);CHKERRQ(IERR)
  CALL VecDuplicate(X_WAVE,B_WAVE,IERR);CHKERRQ(IERR)
# endif

  IF(PMSR) WRITE(PSCRN,*) 'VECTOR SETUP       : COMPLETE!'

  !---------------------------------------------------------
  ! set up vector scatter operations
  !---------------------------------------------------------
  !for insertion (gather):  petsc local --> petsc global
# if defined (NH)
  CALL VecScatterCreate(XLOCAL,ISLOCAL,X,ISGLOBAL,L2G,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL VecScatterCreate(XL_EL,ISLOCAL_EL,X_EL,ISGLOBAL_EL,L2G_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL VecScatterCreate(XL_WAVE,ISLOCAL_WAVE,X_WAVE,ISGLOBAL_WAVE,L2G_WAVE,IERR);CHKERRQ(IERR)
# endif

  !for a true scatter: petsc global --> petsc local
# if defined (NH)
  CALL VecScatterCreate(X,ISGLOBAL,XLOCAL,ISLOCAL,G2L,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL VecScatterCreate(X_EL,ISGLOBAL_EL,XL_EL,ISLOCAL_EL,G2L_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL VecScatterCreate(X_WAVE,ISGLOBAL_WAVE,XL_WAVE,ISLOCAL_WAVE,G2L_WAVE,IERR);CHKERRQ(IERR)
# endif

  IF(PMSR) WRITE(PSCRN,*) 'VECTOR SCATTER OPS : COMPLETE!'

  !---------------------------------------------------------
  ! create the linear solver and setup the options
  !   set type
  !   set initial guess is zero or nonzero
  !   set convergence tolerance
  !   use set from options to allow override of these
  !   and defaults from the command line at runtime
  !---------------------------------------------------------
# if defined (NH)
  CALL KSPCreate(MPI_FVCOM_GROUP,Ksp,IERR);CHKERRQ(IERR)
  CALL KSPSetOperators(Ksp,A,A,SAME_NONZERO_PATTERN,IERR);CHKERRQ(IERR)
  CALL KSPSetType(Ksp,KSPGMRES,IERR);CHKERRQ(IERR)
  CALL KSPGetPC(Ksp,Pc,IERR);CHKERRQ(IERR)
  CALL PCSetType(Pc,PCHYPRE,IERR);CHKERRQ(IERR)
  CALL PCHYPRESetType(Pc,"boomeramg",IERR);CHKERRQ(IERR)
  IF(USE_LAST) CALL KSPSetInitialGuessNonzero(Ksp,PETSC_TRUE,IERR);CHKERRQ(IERR)
  CALL KSPSetTolerances(Ksp,RTOL,PETSC_DEFAULT_DOUBLE_PRECISION,  &
         PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,IERR);CHKERRQ(IERR) 
  CALL KSPSetFromOptions(Ksp,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL KSPCreate(MPI_FVCOM_GROUP,Ksp_EL,IERR);CHKERRQ(IERR)
  CALL KSPSetOperators(Ksp_EL,A_EL,A_EL,SAME_NONZERO_PATTERN,IERR);CHKERRQ(IERR)
  CALL KSPSetType(Ksp_EL,KSPGMRES,IERR);CHKERRQ(IERR)
  CALL KSPGetPC(Ksp_EL,Pc_EL,IERR);CHKERRQ(IERR)
  CALL PCSetType(Pc_EL,PCHYPRE,IERR);CHKERRQ(IERR)
  CALL PCHYPRESetType(Pc_EL,"boomeramg",IERR);CHKERRQ(IERR)
  IF(USE_LAST) CALL KSPSetInitialGuessNonzero(Ksp_EL,PETSC_TRUE,IERR);CHKERRQ(IERR)
  CALL KSPSetTolerances(Ksp_EL,RTOL,PETSC_DEFAULT_DOUBLE_PRECISION,  &
         PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,IERR);CHKERRQ(IERR)
  CALL KSPSetFromOptions(Ksp_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL KSPCreate(MPI_FVCOM_GROUP,Ksp_WAVE,IERR);CHKERRQ(IERR)
!  CALL KSPSetOperators(Ksp_WAVE,A_WAVE,A_WAVE,SAME_NONZERO_PATTERN,IERR);CHKERRQ(IERR)
  CALL KSPSetOperators(Ksp_WAVE,A_WAVE,A_WAVE,DIFFERENT_NONZERO_PATTERN,IERR);CHKERRQ(IERR)
  CALL KSPSetType(Ksp_WAVE,KSPGMRES,IERR);CHKERRQ(IERR)
!  CALL KSPSetType(Ksp_WAVE,KSPBCGS,IERR);CHKERRQ(IERR)
!  CALL KSPGetPC(Ksp_WAVE,Pc_WAVE,IERR);CHKERRQ(IERR)
!  CALL PCSetType(Pc_WAVE,PCHYPRE,IERR);CHKERRQ(IERR)
!  CALL PCHYPRESetType(Pc_WAVE,"boomeramg",IERR);CHKERRQ(IERR)
  IF(USE_LAST) CALL KSPSetInitialGuessNonzero(Ksp_WAVE,PETSC_TRUE,IERR);CHKERRQ(IERR)
  CALL KSPSetTolerances(Ksp_WAVE,RTOL,PETSC_DEFAULT_DOUBLE_PRECISION,  &
         PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,IERR);CHKERRQ(IERR)
  CALL KSPSetFromOptions(Ksp_WAVE,IERR);CHKERRQ(IERR)
# endif

  IF(PMSR) WRITE(PSCRN,*) 'SOLVER CONTEXT     : COMPLETE!'

  !---------------------------------------------------------
  ! deallocate unneeded arrays
  ! we need to keep the scatters (l2g,g2l)
  ! we need to keep the local to global matrix map (isl2g)
  ! we need to keep both local-local maps: alo_2_plo_node
  !                                        plo_2_alo_node
  ! we may or may need the index sets (islocal,isglobal)
  !---------------------------------------------------------
  DEALLOCATE(ACCUM)
  DEALLOCATE(PARTID)
  DEALLOCATE(PART_NODES)
  DEALLOCATE(PLO_2_AGO_NODE)
  DEALLOCATE(AGO_2_PGO_NODE)
  DEALLOCATE(CNT)
!  DEALLOCATE(PLO_2_PGO_NODE)
# if defined (NH)
  DEALLOCATE(PLO_2_PGO)
# endif

  !report and reset logical first_entry
  IF(PMSR) WRITE(PSCRN,*) 'PETSc SETUP        : COMPLETE!'
  IF(PDEBUG) WRITE(PUNIT,*) 'FINISHING: ',TRIM(SUBNAME)

  RETURN
  END SUBROUTINE PETSc_SET

!=========================================================================================
  SUBROUTINE PETSc_AllocA
!============================================================================
!                         PETSc_AllocA
!     preallocate matrix A
!     assume our nonzero pattern will  not change
!     this will make the implementation much more efficient
!     malloc-ing on the fly can be very slow
!
!     currently assumes stencil is:
!       -node       + surrounding nodes
!       -node above + surrounding nodes
!       -node below + surrounding nodes
!
!     typical number of nonzero entries:  21
!
!     dependencies:  all internal to this module
!          except:  m,kbm1
!============================================================================
  USE ALL_VARS, ONLY : M, KBM1, NTSN, NBSN, MGL
  USE CONTROL,  ONLY : MPI_FVCOM_GROUP
  IMPLICIT NONE
  INTEGER :: NODE,LAY,PETSc_POS,NCOL,NCOL2,NEYNUM,NEY,COLCOUNT,I
  INTEGER :: COL2(MAX_NZ_ROW)
  INTEGER,ALLOCATABLE  :: ONNZ2(:)
  INTEGER,ALLOCATABLE  :: DNNZ2(:)
  CHARACTER(LEN=20)    :: SUBNAME = 'PETSc_AllocA'

# if defined (NH)
  INTEGER :: COL(MAX_NZ_ROW)
  INTEGER,ALLOCATABLE  :: ONNZ(:), DNNZ(:)
# endif

  PetscInt :: IERR

  IF(PDEBUG) WRITE(PUNIT,*) 'STARTING: ', TRIM(SUBNAME)
  !-----------------------------------------------------------------------------
  ! assemble matrix, use local vars (requires that MatSetLocalToGlobalMapping is setup)
  ! note: row and column are C order (starts from 0)
  ! see page 51 of petsc manual, bottom paragraph
  ! we will set entries row by row
  !-----------------------------------------------------------------------------
# if defined (NH)
  ALLOCATE(ONNZ(PSIZE))    ; ONNZ = 0
  ALLOCATE(DNNZ(PSIZE))    ; DNNZ = 0
# endif
  ALLOCATE(ONNZ2(N_VERTS)) ; ONNZ2 = 0
  ALLOCATE(DNNZ2(N_VERTS)) ; DNNZ2 = 0

  DO NODE=1, M

    COL2  = 0
    NCOL2 = 0
    IF( ALO_2_PLO_NODE(NODE) <= N_VERTS) THEN
      NCOL2 = 1
      COL2(1) = ALO_2_PLO_NODE(NODE)

      DO NEYNUM =1, NTSN(NODE)-1
        NEY = NBSN(NODE,NEYNUM)
        IF(NEY == NODE) CYCLE
        NCOL2 = NCOL2 + 1
        !make sure we havent exceeded max_nz_row
        IF(NCOL2 > 10) THEN
          WRITE(PUNIT,*) 'WE HAVE EXCEEDED MAX_NZ_ROW==10'
          WRITE(PUNIT,*) 'INCREASE AND RECOMPILE'
          CALL MPI_FINALIZE(IERR)
          STOP
        ENDIF
        COL2(NCOL2) = ALO_2_PLO_NODE(NEY)
      ENDDO

      DO COLCOUNT =1, NCOL2
        IF(COL2(COLCOUNT) <= N_VERTS) THEN
          ONNZ2(ALO_2_PLO_NODE(NODE)) = ONNZ2(ALO_2_PLO_NODE(NODE)) + 1
        ELSE
          DNNZ2(ALO_2_PLO_NODE(NODE)) = DNNZ2(ALO_2_PLO_NODE(NODE)) + 1
        ENDIF
      ENDDO

    ENDIF

#   if defined (NH)
    DO LAY =1, KBM1

      !get node number in petsc local domain, cycle if not in domain
      PETSc_POS = ALO_2_PLO(NODE,LAY)
      IF(PETSc_POS > PSIZE) CYCLE

      !initialize column numbers and row values
      COL = 0
      NCOL = 1 !diagonal entry

      !set row number in PLO.  If it is not in PLO interior, just cycle
      !this means it is a boundar node not in the PLO partition and 
      !another partition will be setting that row
      COL(1) = ALO_2_PLO(NODE,LAY)
!      IF(COL(1) > PSIZE) CYCLE

      !set discretization in current layer" 
      DO NEYNUM =1, NTSN(NODE)-1
        NEY = NBSN(NODE,NEYNUM)
        IF(NEY == NODE) CYCLE
        NCOL = NCOL + 1
        !make sure we havent exceeded max_nz_row
        IF(NCOL > MAX_NZ_ROW) THEN
          WRITE(PUNIT,*) 'WE HAVE EXCEEDED MAX_NZ_ROW'
          WRITE(PUNIT,*) 'INCREASE AND RECOMPILE'
          CALL MPI_FINALIZE(IERR) 
          STOP
        ENDIF
        COL(NCOL) = ALO_2_PLO(NEY,LAY)
      ENDDO

      !set discretization in layer above
      IF(LAY > 1) THEN
        DO NEYNUM = 1, NTSN(NODE)-1
          NEY = NBSN(NODE,NEYNUM)
          IF(NEY == NODE) CYCLE
          NCOL = NCOL + 1
          !make sure we havent exceeded max_nz_row
          IF(NCOL > MAX_NZ_ROW) THEN
            WRITE(PUNIT,*) 'WE HAVE EXCEEDED MAX_NZ_ROW'
            WRITE(PUNIT,*) 'INCREASE AND RECOMPILE'
            CALL MPI_FINALIZE(IERR)
            STOP
          ENDIF
          COL(NCOL) = ALO_2_PLO(NEY,LAY-1)
        ENDDO
      ENDIF

      !set discretization in layer below
      IF(LAY < N_LAY) THEN
        DO NEYNUM = 1, NTSN(NODE)-1
          NEY = NBSN(NODE,NEYNUM)
          IF (NEY == NODE) CYCLE
          NCOL = NCOL + 1
          !make sure we havent exceeded max_nz_row
          IF(NCOL > MAX_NZ_ROW) THEN
            WRITE(PUNIT,*) 'WE HAVE EXCEEDED MAX_NZ_ROW'
            WRITE(PUNIT,*) 'INCREASE AND RECOMPILE'
            CALL MPI_FINALIZE(IERR)
            STOP
          ENDIF
          COL(NCOL) = ALO_2_PLO(NEY,LAY+1)
        ENDDO
      ENDIF

      !set vertical "discretization"
      IF(LAY > 1) THEN
        NCOL = NCOL + 1
        !make sure we havent exceeded max_nz_row
        IF(NCOL > MAX_NZ_ROW) THEN
          WRITE(PUNIT,*) 'WE HAVE EXCEEDED MAX_NZ_ROW'
          WRITE(PUNIT,*) 'INCREASE AND RECOMPILE'
          CALL MPI_FINALIZE(IERR) 
          STOP
        ENDIF
        COL(NCOL) = ALO_2_PLO(NODE,LAY-1)
      ENDIF

      IF(LAY < KBM1) THEN
        NCOL = NCOL + 1
        !make sure we havent exceeded max_nz_row
        IF(NCOL > MAX_NZ_ROW) THEN
          WRITE(PUNIT,*) 'WE HAVE EXCEEDED MAX_NZ_ROW'
          WRITE(PUNIT,*) 'INCREASE AND RECOMPILE'
          CALL MPI_FINALIZE(IERR) 
          STOP
        ENDIF
        COL(NCOL) = ALO_2_PLO(NODE,LAY+1)
      ENDIF

      !update diagonal block nonzero counter for this row if it is in petsc local domain
      DO COLCOUNT = 1, NCOL
        IF(COL(COLCOUNT) <= PSIZE) THEN
          DNNZ(PETSc_POS) = DNNZ(PETSc_POS) + 1
        ELSE
          ONNZ(PETSc_POS) = ONNZ(PETSc_POS) + 1
        ENDIF
      ENDDO
      
    ENDDO
#   endif

  ENDDO

# if defined (NH)
  !preallocate the matrix 
  IF(PDEBUG) THEN
    WRITE(PUNIT,*)'PREALLOCATING MATRIX'
    WRITE(PUNIT,*)'N_VERTS: ',N_VERTS,' MGL: ',MGL
    WRITE(PUNIT,*)' COUNT   DNNZ   ONNZ'
    DO I=1, PSIZE   
      WRITE(PUNIT,*) I,DNNZ(I),ONNZ(I)
    ENDDO
  ENDIF
# endif

# if defined (NH)
  CALL MatCreateMPIAIJ(MPI_FVCOM_GROUP,PSIZE,PSIZE,PsizeGL,PsizeGL,0,DNNZ,0,ONNZ,A,IERR);CHKERRQ(IERR)
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  CALL MatCreateMPIAIJ(MPI_FVCOM_GROUP,N_VERTS,N_VERTS,MGL,MGL,0,DNNZ2,0,ONNZ2,A_EL,IERR);CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL MatCreateMPIAIJ(MPI_FVCOM_GROUP,N_VERTS,N_VERTS,MGL,MGL,0,DNNZ2,0,ONNZ2,A_WAVE,IERR);CHKERRQ(IERR)
# endif

  !deallocate non-zero counters
# if defined (NH)
  DEALLOCATE(ONNZ,DNNZ)
# endif
  DEALLOCATE(ONNZ2,DNNZ2)

  IF(PDEBUG) WRITE(PUNIT,*)'FINISHING: ',TRIM(SUBNAME)

  RETURN
  END SUBROUTINE PETSc_AllocA
!=============================================================================

# if defined (NH)
  SUBROUTINE PETSc_SOLVER
  !============================================================================
  !                         petsc_solve
  ! solve Ax=b, report convergence success, iteration count
  ! option to initialize b using exact solution x to check solution
  !
  ! basic function:  solve for a new x
  !
  ! dependencies:  pmsr,u,A,b,one,ierr,ksb,x,neg_one,check_exact from this module
  !                norm,its from this module
  !
  !============================================================================
  IMPLICIT NONE
  CHARACTER(LEN=20) :: SUBNAME = 'PETSc_SOLVER'
  PetscInt :: REASON, IERR

  IF(PDEBUG) WRITE(PUNIT,*) 'STARTING: ', TRIM(SUBNAME)
  !---------------------------------------------------------
  ! option: use an exact solution of u=1 to create a RHS
  !---------------------------------------------------------
  IF(CHECK_EXACT) THEN
    CALL VecSet(U,ONE,IERR);CHKERRQ(IERR)
    CALL VecSet(X,ZERO,IERR);CHKERRQ(IERR)
    CALL MatMult(A,U,B,IERR);CHKERRQ(IERR)
  ENDIF

  !---------------------------------------------------------
  ! solve the linear system Ax=b and gather info
  ! about convergence rates, iterations, etc.
  !---------------------------------------------------------
  CALL KSPSolve(Ksp,B,X,IERR);CHKERRQ(IERR)

  !check what happened (divergence/convergence) and report
  CALL KSPGetConvergedReason(Ksp,REASON,IERR)
  IF(REASON < 0) THEN
    IF(PMSR) THEN
      WRITE (PSCRN,*) 'PETSc SOLVER HAS DIVERGED: STOPPING...'
    ENDIF
    CALL MPI_FINALIZE(IERR)
    STOP
  ELSE
    CALL KSPGetIterationNumber(Ksp,ITS,IERR);CHKERRQ(IERR)
    !qxu IF(PMSR) WRITE(PSCRN,*)'PETSc CONVERGED IN : ',ITS,'AT RATE ',RTOL**(1.0d0/FLOAT(ITS))
  ENDIF

  !option to check solution for exact case
  IF(CHECK_EXACT) THEN
    CALL VecAXPY(X,NEG_ONE,U,IERR);CHKERRQ(IERR)
    CALL VecNorm(X,NORM_2,NORM,IERR);CHKERRQ(IERR)
    IF(PMSR) WRITE(PUNIT,*) 'ERROR NORM FROM EXACT TEST: ',NORM/FLOAT(PSIZEGL)
    IF(PMSR) WRITE(PSCRN,*) 'ERROR NORM         : ',NORM/FLOAT(PSIZEGL)
  ELSE
    CALL VecNorm(X,NORM_2,NORM,IERR);CHKERRQ(IERR)
    IF(PMSR) WRITE(PSCRN,*) 'SOLUTION NORM      : ',NORM
  ENDIF

  IF(PDEBUG) WRITE(PUNIT,*) 'FINISHING: ',TRIM(SUBNAME)

  RETURN
  END SUBROUTINE PETSc_SOLVER
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  SUBROUTINE PETSc_SOLVER_EL
  IMPLICIT NONE
  PetscInt :: REASON, IERR

  CALL KSPSolve(Ksp_EL,B_EL,X_EL,IERR);CHKERRQ(IERR)

  CALL KSPGetConvergedReason(Ksp_EL,REASON,IERR)
  IF(REASON < 0) THEN
    IF(PMSR) THEN
      WRITE (PSCRN,*) 'PETSc SOLVER HAS DIVERGED IN EL CAL: STOPPING...'
    ENDIF
    CALL MPI_FINALIZE(IERR)
    STOP
  ELSE
    CALL KSPGetIterationNumber(Ksp_EL,ITS_EL,IERR);CHKERRQ(IERR)
!qxu    IF(PMSR) WRITE(PSCRN,*)'PETSc CONVERGED IN : ',ITS_EL,'AT RATE ',RTOL**(1.0d0/FLOAT(ITS_EL))
  ENDIF

  CALL VecNorm(X_EL,NORM_2,NORM_EL,IERR);CHKERRQ(IERR)
!qxu  IF(PMSR) WRITE(PSCRN,*) 'SOLUTION OF EL NORM      : ',NORM_EL

  RETURN
  END SUBROUTINE PETSc_SOLVER_EL
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  SUBROUTINE PETSc_SOLVER_WAVE
  IMPLICIT NONE
  PetscInt :: REASON
  PetscInt :: IERR

  CALL KSPSolve(Ksp_WAVE,B_WAVE,X_WAVE,IERR);CHKERRQ(IERR)

  CALL KSPGetConvergedReason(Ksp_WAVE,REASON,IERR)
  IF(REASON < 0)THEN
    IF(PMSR)THEN
      WRITE(PSCRN,*) 'PETSc SOLVER HAS DIVERGED IN WAVE CAL: STOPPING...'
    END IF
    CALL MPI_FINALIZE(IERR)
    STOP
  ELSE
    CALL KSPGetIterationNumber(Ksp_WAVE,ITS_WAVE,IERR);CHKERRQ(IERR)
    ! LWU IF(PMSR) WRITE(PSCRN,*)'PETSc CONVERGED IN : ',ITS_WAVE,'AT RATE ',RTOL**(1.0d0/FLOAT(ITS_WAVE))
  END IF

  CALL VecNorm(X_WAVE,NORM_2,NORM_WAVE,IERR);CHKERRQ(IERR)
  ! LWU IF(PMSR) WRITE(PSCRN,*) 'SOLUTION OF WAVE NORM      : ',NORM_WAVE

  RETURN
  END SUBROUTINE PETSc_SOLVER_WAVE
# endif

  SUBROUTINE PETSc_CLEANUP
  !============================================================================
  !                         petsc_cleanup
  !      cleanup memory from petsc vars (called at end of FVCOM run)
  !
  ! external deps:  ksp,u,x,b,blocal,xlocal,A,ierr,pmsr from this module
  !============================================================================
  IMPLICIT NONE
  PetscInt :: IERR

# if defined (NH)
  CALL KSPDestroy(Ksp,IERR)   ;CHKERRQ(IERR)
  CALL VecDestroy(U,IERR)     ;CHKERRQ(IERR)
  CALL VecDestroy(X,IERR)     ;CHKERRQ(IERR)
  CALL VecDestroy(B,IERR)     ;CHKERRQ(IERR)
  CALL VecDestroy(BLOCAL,IERR);CHKERRQ(IERR)
  CALL VecDestroy(XLOCAL,IERR);CHKERRQ(IERR)
  CALL MatDestroy(A,IERR)     ;CHKERRQ(IERR)
# endif  

# if defined (SEMI_IMPLICIT) || (NH) 
  CALL KSPDestroy(Ksp_EL,IERR)   ;CHKERRQ(IERR)
  CALL VecDestroy(X_EL,IERR)     ;CHKERRQ(IERR)
  CALL VecDestroy(B_EL,IERR)     ;CHKERRQ(IERR)
  CALL VecDestroy(BL_EL,IERR)    ;CHKERRQ(IERR)
  CALL VecDestroy(XL_EL,IERR)    ;CHKERRQ(IERR)
  CALL MatDestroy(A_EL,IERR)     ;CHKERRQ(IERR)
# endif

# if defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT)
  CALL KSPDestroy(Ksp_WAVE,IERR)   ;CHKERRQ(IERR)
  CALL VecDestroy(X_WAVE,IERR)     ;CHKERRQ(IERR)
  CALL VecDestroy(B_WAVE,IERR)     ;CHKERRQ(IERR)
  CALL VecDestroy(BL_WAVE,IERR)    ;CHKERRQ(IERR)
  CALL VecDestroy(XL_WAVE,IERR)    ;CHKERRQ(IERR)
  CALL MatDestroy(A_WAVE,IERR)     ;CHKERRQ(IERR)
# endif

  CALL petscfinalize(IERR)       ;CHKERRQ(IERR)

  IF(PMSR) WRITE(PSCRN,*) '    '
  IF(PMSR) WRITE(PSCRN,*) 'PETSc CLEANUP      : COMPLETE!'

  RETURN
  END SUBROUTINE PETSc_CLEANUP

  FUNCTION ALO_2_PLO(NODE,LAY) RESULT(PLO_LOC)
  !============================================================================
  !                      alo_2_plo
  ! map an unknown at node [node] and layer [lay] in application local
  ! ordering to a point in petsc local ordering in the full set of petsc
  ! local unknowns [1 --> (n_verts+n_halo)*(N_lay)].
  ! we will assume that layers run fastest in our vector
  !
  ! external deps:  N_lay from this module
  !============================================================================
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NODE
  INTEGER, INTENT(IN) :: LAY
  INTEGER             :: PLO_LOC,PLO_NODE
  
  PLO_NODE = ALO_2_PLO_NODE(NODE)
  PLO_LOC  = LAY + (PLO_NODE-1)*N_LAY

  RETURN
  END FUNCTION ALO_2_PLO

  FUNCTION PLO_2_ALO(PLO_LOC) RESULT(ALO_LOC)
  !============================================================================
  !                      plo_2_alo
  ! map an unknown at petsc local unknown plo_loc to a particular node and
  ! layer in the application local ordering
  ! ordering to a point in petsc local ordering in the full set of petsc
  ! local unknowns [1 --> (n_verts+n_halo)*(N_lay)].
  ! we will assume that layers run fastest in our vector
  !
  ! external deps: N_lay from this module
  !============================================================================
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: PLO_LOC
  INTEGER             :: ALO_LOC(2)
  INTEGER             :: LAY,NODE,PLO_NODE

  !get the local plo node and layer
  LAY = N_LAY
  IF(MOD(PLO_LOC,N_LAY) /= 0) LAY = MOD(PLO_LOC,N_LAY)
  PLO_NODE = (PLO_LOC-LAY)/N_LAY+1
  NODE     = PLO_2_ALO_NODE(PLO_NODE)

  !set result
  ALO_LOC(1) = NODE
  ALO_LOC(2) = LAY

  RETURN
  END FUNCTION PLO_2_ALO

# endif

!=========================================================================================
!
!=========================================================================================
END MODULE MOD_PETSc

