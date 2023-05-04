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
! !MODULE: ice_work - globally accessible, temporary work arrays
!
! !DESCRIPTION:
!
! The intent is to save memory by allocating global arrays only when
! necessary.  Globally accessible, local (i.e., on-processor) work 
! arrays are also available to conserve memory.  These arrays should 
! be used only within a single subroutine!
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      module ice_work
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
!
!EOP
!
      implicit none

      ! global
      real (kind=dbl_kind), dimension(:,:), allocatable :: work_g1
      real (kind=dbl_kind), dimension(:,:), allocatable :: work_g2
      real (kind=real_kind),dimension(:,:), allocatable :: work_gr

      ! local
!      real (kind=dbl_kind) :: 
!     &   work_l1(imt_local,jmt_local)
!     &,  work_l2(imt_local,jmt_local)
!     &,  worka(ilo:ihi,jlo:jhi)
!     &,  workb(ilo:ihi,jlo:jhi)

      real (kind=dbl_kind), dimension(:,:), allocatable  :: &
         work_l1  ,  work_l2 ,  worka ,  workb


!=======================================================================

      end module ice_work

!=======================================================================
