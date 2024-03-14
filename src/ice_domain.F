!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!
!BOP
!
! !MODULE: ice_domain - sets array sizes for local subdomain and parallel info
!
! !DESCRIPTION:
!
! Sets array sizes for local subdomain and related parallel 
! processing info. Code originally based on domain.F in POP
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      module ice_domain
!
! !USES:
!
        use ice_kinds_mod
        use ice_model_size

!      connect to FVCOM

# if defined (MULTIPROCESSOR)
       use lims   !, only :MT
       USE MOD_PREC
       USE MOD_PAR
#endif

!
!EOP
!
      implicit none
      save

      ! ghost cells define boundary data to avoid frequent 
      ! off-processor communication 
      integer (kind=int_kind), parameter ::  &
       num_ghost_cells = 1

      ! imt_local, jmt_local define total size of subdomain - 
      ! the formula below automatically pads arrays if NPROC_X
      ! or NPROC_Y do not evenly divide global problem size

      !!!!NOTE:  This version of the ice model is not set up to !!!!
      !!!!       use the padding properly.  The number of       !!!!
      !!!!       processors MUST divide evenly into the grid.   !!!!

!      integer (kind=int_kind), parameter ::  &
!       imt_local =  1   !, & 
!       jmt_local =   100  ! MT  ggao need change to MT

!     &  imt_local = (imt_global-1)/NPROC_X + 1 + 2*num_ghost_cells,
!     &  jmt_local = (jmt_global-1)/NPROC_Y + 1 + 2*num_ghost_cells

!      integer (kind=int_kind), parameter :: &
!       ilo = 1, &  ! beg index of actual physical subdomain
!       ihi = 1, &  ! end index 
!       jlo = 1   ! beg index
!       ,jhi = MT     !  ggao need change to MT   ! end index 

!    &  ilo = num_ghost_cells+1, ! beg index of actual physical subdomain
!    &  ihi = imt_local - num_ghost_cells, ! end index
!    &  jlo = num_ghost_cells+1,           ! beg index
!    &  jhi = jmt_local - num_ghost_cells  ! end index

!      integer (kind=int_kind), dimension (ilo:ihi,jlo:jhi) :: &
       integer (kind=int_kind), dimension (:,:),allocatable,save :: &
       index_global          ! global index for local subdomain (int)


!      integer (kind=int_kind), allocatable :: &
!       index_global(:,:)           ! global index for local subdomain (int)

!       real (kind=dbl_kind), allocatable ::  &
       integer (kind=int_kind), dimension (:,:),allocatable,save :: &
       rndex_global           ! global index for local subdomain (dbl)

       integer (kind=int_kind),save  :: &
       my_task,      &          ! task id for local process
       master_task            ! task id for master process



      integer (kind=int_kind) :: ilo,ihi,jlo,jhi,imt_local,jmt_local

      end module ice_domain
 
!=======================================================================
