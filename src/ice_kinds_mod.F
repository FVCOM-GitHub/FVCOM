!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!=======================================================================
!BOP
!
! !MODULE: ice_kinds_mod - defines variable precision
!
! !DESCRIPTION:
!
! Defines variable precision for all common data types \\
! Code originally based on kinds_mod.F in POP
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !INTERFACE:
!
      module ice_kinds_mod
!
! !USES:
!
!EOP
!=======================================================================
!!  ggao 0104/2008

     use MOD_PREC

      implicit none
      save

!      integer, parameter :: char_len  = 80,  &
!                           char_len_long  = 128, & 
!                           int_kind  = kind(1),  &
!                           log_kind  = kind(.true.), &
!     !                    real_kind = selected_real_kind(6),&
!     !                    dbl_kind  = selected_real_kind(13)
!
!                         real_kind = SP,&
!                         dbl_kind  = DP

#  if !defined (DOUBLE_PRECISION)

      integer, parameter :: char_len  = 80,  &
                           char_len_long  = 128, &
                           int_kind  = kind(1),  &
                           log_kind  = kind(.true.), &
!     !                    real_kind = selected_real_kind(6),&
!     !                    dbl_kind  = selected_real_kind(13)

                         real_kind = SP,&
                         dbl_kind  = SP
#  else
!--Double Precision Coding------------------------------------------------------!

      integer, parameter :: char_len  = 80,  &
                           char_len_long  = 128, &
                           int_kind  = kind(1),  &
                           log_kind  = kind(.true.), &
!     !                    real_kind = selected_real_kind(6),&
!     !                    dbl_kind  = selected_real_kind(13)

                         real_kind = DP,&
                         dbl_kind  = DP


#  endif
 


!=======================================================================

      end module ice_kinds_mod

!=======================================================================
