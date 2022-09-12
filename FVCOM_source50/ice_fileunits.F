!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_fileunits 
!
! !DESCRIPTION:
!
! Defines unit numbers for files opened for reading or writing
!
! !REVISION HISTORY:
!
!  author: Elizabeth C. Hunke
!          Fluid Dynamics Group, Los Alamos National Laboratory
!
! !INTERFACE:
!
      module ice_fileunits
!
! !USES:
      use ice_kinds_mod
      USE CONTROL, only: nu_diag => ipt
!
!EOP
!=======================================================================

      implicit none
      save

      ! DAVID Changed to used the same log output file as the fvcom
      ! main code!

      integer (kind=int_kind), parameter :: &
         nu_grid        = 81  &! grid file
      ,  nu_kmt         = 82  &! land mask file
      ,  nu_nml         = 83  &! namelist input file
!      ,  nu_diag        = 6   &! diagnostics output file
      ,  nu_forcing     = 84  &! forcing data file
      ,  nu_dump        = 85  &! dump file for restarting
      ,  nu_restart     = 85  &! restart input file
      ,  nu_rst_pointer = 86   ! pointer to latest restart file

!=======================================================================

      end module ice_fileunits

!=======================================================================
