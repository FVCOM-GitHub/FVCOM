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

!===============================================================================!
! DEFINE FLOATING POINT PRECISION USING KIND                                    !
!===============================================================================!
MODULE MOD_PREC
#  if defined (MULTIPROCESSOR)
   use mpi
#  endif

   IMPLICIT NONE

!--Single Precision Coding------------------------------------------------------!
#  if !defined (DOUBLE_PRECISION)
   INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,30)
#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_F = MPI_REAL
# endif
#  else 
!--Double Precision Coding------------------------------------------------------!
   INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(12,300)
#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_F = MPI_DOUBLE_PRECISION
# endif
#  endif

  
   INTEGER, PARAMETER :: SPA = SELECTED_REAL_KIND(6,30)
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(12,300)

#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_DP = MPI_DOUBLE_PRECISION
# endif


END MODULE MOD_PREC

