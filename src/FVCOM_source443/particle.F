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

!==============================================================================!
!                LAGRANGIAN PARTICLE CLASS                                     !
!==============================================================================!

! MODFIED BY DAVID FOR THE NEW INPUT: NEEDS TESTING

! ALL TIMES ARE NOW EITHER TYPE(TIME) OR IN DAYS!
  
module particle_class
  use mod_prec
  use mod_time
  use control, only : zerotime
  USE MOD_UTILS

  ! INTEL 9.1.X does not like initialized types
  integer, parameter :: number_of_floats = 26
  integer, parameter :: number_of_scalars = 5

!!$  type particle 
!!$    integer  :: id          !!global particle number
!!$    integer  :: pid         !!processor id number        
!!$    integer  :: elem        !!element containing particle
!!$    integer  :: group       !!element group id
!!$    TYPE(TIME) :: tbeg      !! START TIME
!!$    TYPE(TIME) :: tend      !! END TIME
!!$    real(sp) :: x(3)           !!particle position
!!$    real(sp) :: xn(3)          !!last time step position 
!!$    real(sp) :: el          !!surface elevation at particle point
!!$    real(sp) :: h           !!bathymetry at particle point
!!$    real(sp) :: u           !!x-velocity at particle location
!!$    real(sp) :: v           !!y-velocity at particle location
!!$    real(sp) :: w           !!sigma-velocity at particle location
!!$    real(sp) :: zloc        !!z position of particle [m]            
!!$    real(sp) :: s(number_of_scalars) !!particle scalar
!!$    real(sp) :: chi(3,4)    !!Runge-Kutta stage contributions
!!$    real(sp) :: deltat      !!particle time step
!!$    real(sp) :: pathlength  !!particle integrated pathlength [m]
!!$    logical  :: found      
!!$  end type particle 

  ! USE THIS FOR INTEL 10.1+ 
  type particle 
    integer  :: id            =0 !!global particle number
    integer  :: pid           =0 !!processor id number        
    integer  :: elem          =0 !!element containing particle
    integer  :: group         =0 !!element group id
    TYPE(TIME) :: tbeg        !! START TIME
    TYPE(TIME) :: tend        !! END TIME
    real(sp) :: x(3)          =0.0_sp !!particle position
    real(sp) :: xn(3)         =0.0_sp !!last time step position 
    real(sp) :: el            =0._SP !!surface elevation at particle point
    real(sp) :: h             =0._SP !!bathymetry at particle point
    real(sp) :: u             =0._SP !!x-velocity at particle location
    real(sp) :: v             =0._SP !!y-velocity at particle location
    real(sp) :: w             =0._SP !!sigma-velocity at particle location
    real(sp) :: zloc          =0._SP !!z position of particle [m]            
    real(sp) :: s(number_of_scalars) =0._SP !!particle scalar
    real(sp) :: chi(3,4)      =0._SP !!Runge-Kutta stage contributions
    real(sp) :: deltat        =0._SP !!particle time step
    real(sp) :: pathlength    =0._SP !!particle integrated pathlength [m]
    logical  :: found         = .false.
  end type particle 

  INTEGER MPI_PARTICLE

  interface operator (<)    !for sorting, insert
    module procedure less_than_particle ; end interface
  interface operator (==)   !for sorting or delete
    module procedure equal_to_particle  ; end interface

  contains

  function less_than_particle (particle_1,particle_2) result(Boolean)
    type(particle), intent(in) :: particle_1
    type(particle), intent(in) :: particle_2
    logical                    :: Boolean
    Boolean = particle_1%id < particle_2%id
  end function less_than_particle

  function equal_to_particle (particle_1,particle_2) result(Boolean)
    type(particle), intent(in) :: particle_1
    type(particle), intent(in) :: particle_2
    logical                    :: Boolean
    Boolean = particle_1%id == particle_2%id
  end function equal_to_particle

  subroutine screen_write(p,hprint)
    type(particle), intent(in) :: p
    logical,        intent(in) :: hprint
    if(hprint)then
      write(ipt,*)'id   group    x           y         z      elem pid'
    endif
    write(ipt,'(2I5,3F10.2,I8,2I5)')p%id,p%group,p%x,p%elem,p%pid
  end subroutine screen_write

  subroutine shift_pos(p)
    type(particle), intent(inout) :: p
    p%xn = p%x
  end subroutine shift_pos 

!!$  subroutine set_pathlength(p)
!!$    type(particle), intent(inout) :: p
!!$    p%pathlength = sqrt( (p%x-p%xn)**2)
!!$  end subroutine set_pathlength

  !dump particle info to screen
  subroutine particle_print(p)
    type(particle), intent(in) :: p
!    type(particle), pointer :: p
    write(ipt,*)
    write(ipt,*)'id:        ',p%id
    write(ipt,*)'processor: ',p%pid
    write(ipt,*)'element:   ',p%elem
    write(ipt,*)'group:     ',p%group
    write(ipt,*)'tbeg:      ',p%tbeg
    write(ipt,*)'tend:      ',p%tend
    write(ipt,*)'x:         ',p%x
    write(ipt,*)'xn:        ',p%xn
    write(ipt,*)'el:        ',p%el
    write(ipt,*)'h:         ',p%h
    write(ipt,*)'u:         ',p%u
    write(ipt,*)'v:         ',p%v
    write(ipt,*)'w:         ',p%w
    write(ipt,*)'zloc:      ',p%zloc
    write(ipt,*)'s:         ',p%s
    write(ipt,*)'chi1       ',p%chi(:,1)
    write(ipt,*)'chi2       ',p%chi(:,2)
    write(ipt,*)'chi3       ',p%chi(:,3)
    write(ipt,*)'chi4       ',p%chi(:,4)
    write(ipt,*)'deltat     ',p%deltat   
    write(ipt,*)'pathlength ',p%pathlength
  end subroutine particle_print

  !initialize values
  subroutine zero_out(p)
    implicit none
    type(particle), intent(inout) :: p

    ! SET VALUES
    p%id     = 0
    p%pid    = 0 
    p%group  = 0
    p%elem   = 0
    p%tbeg   = ZEROTIME
    p%tend   = ZEROTIME
    p%x      = 0.0_SP  
    p%xn     = 0.0_SP  
    p%el     = 0.0_SP 
    p%u      = 0.0_SP 
    p%v      = 0.0_SP 
    p%w      = 0.0_SP 
    p%zloc   = 0.0_SP 
    p%s      = 0.0_SP 
    p%chi    = 0.0_SP 
    p%deltat = 0.0
    p%pathlength = 0.0
    p%found  = .false.
  end subroutine zero_out

  SUBROUTINE NEW_PARTICLES(P,i)
    implicit none
    integer, intent(in) :: i
    type(particle), allocatable, intent(inout) :: p(:)
    integer :: status, J
    
    allocate(p(i),stat=status)
    if(status/=0) CALL FATAL_ERROR("NEW_PARTICLE: COULD NOT ALLOCATE?")


    DO J=1,i
       call zero_out(P(j))
    END DO

  END SUBROUTINE NEW_PARTICLES

  
  SUBROUTINE CREATE_MPI_PARTICLE
    USE MOD_PAR
    IMPLICIT NONE
# if defined (MULTIPROCESSOR)
    
    TYPE(TYPE_DEF) :: PRTCL

    ! ZERO THE TYPE
    PRTCL = INIT_TYPE_DEF()

    !ADD EACH COMPONENT
    CALL ADD_TO_MPI_TYPE(PRTCL,MPI_INTEGER,4) ! 5 integers

    CALL ADD_TO_MPI_TYPE(PRTCL,MPI_DOUBLE_PRECISION,4) ! 2 time types(2XLongLong)
    ! CHEATING - LONG-LONG is the same size as DOUBLE

    ! MPI_F is really MPI_DOUBLE_PRECISION in when if defined double_precision
    CALL ADD_TO_MPI_TYPE(PRTCL,MPI_F,number_of_floats + number_of_scalars) ! 31 FLOATING POINT NUMBERS

    CALL ADD_TO_MPI_TYPE(PRTCL,MPI_LOGICAL,1)
    
    ! DEFINE AND COMMIT THE TYPE
    MPI_PARTICLE = CREATE_MPI_TYPE(PRTCL)

  
# endif  
  END SUBROUTINE CREATE_MPI_PARTICLE

end module particle_class
