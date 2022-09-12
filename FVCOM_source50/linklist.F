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

!---------------------------------------------------------
! A linklist module
!---------------------------------------------------------
module linked_list
  use lims, only : myid
  use particle_class
  implicit none

  type link_node
!    private
    type(particle)                 :: v
    type(link_node), pointer       :: next
  end type link_node 
   
  type link_list
!    private 
    type(link_node), pointer :: first
  end type link_list

  contains
   
  subroutine node_delete (links,obj,found)
   implicit none
   type(link_list), intent(inout) :: links 
   type(particle),  intent(in)    :: obj
   logical        , intent(out)   :: found
   type(link_node), pointer       :: previous,current

   !find location of obj
   previous => links%first
   current  => previous%next
   found = .false.

   do
     if(found .or. (.not. associated (current))) return
     if(obj == current%v)then
       found = .true. ; exit
     else
       previous => previous%next  
       current  => current%next
     endif
  end do !find location of node with obj

  if (found) then
    previous%next => current%next
    deallocate(current)
  endif

  end subroutine node_delete

  subroutine delete_NOT_MINE(links,ME)
    implicit none
    type(link_list), intent(inout) :: links
    integer, intent(in)            :: ME
    type(link_node), pointer       :: previous,current
    
    !find location of obj
    previous => links%first
    current  => previous%next
    

    do
       if(.not. associated (current)) return
       if(MYID /= current%v%PID)then
          previous%next => current%next
          deallocate(current)
          current  => previous%next
       else
          previous => previous%next  
          current  => current%next
       endif
    end do !find location of node with obj
    
  end subroutine delete_NOT_MINE

  subroutine delete_not_found (links)
    implicit none
    type(link_list), intent(inout) :: links 
    type(link_node), pointer       :: previous,current
    
    !find location of obj
    previous => links%first
    current  => previous%next
    
    do
       if(.not. associated (current)) return
       if(.not. current%v%found)then
          previous%next => current%next
          deallocate(current)
       else
          previous => previous%next  
          current  => current%next
       endif
    end do !find location of node with obj
    
  end subroutine delete_not_found

  subroutine node_insert( links,obj )
   type(link_list), intent(inout) :: links 
   type(particle),  intent(in)    :: obj
   type(link_node), pointer       :: previous,current

   previous => links%first
   current  => previous%next

   do 
    if( .not. associated (current) )exit
      if( obj < current%v ) exit
        previous => current
        current  => current%next
   end do

   !insert before current 
   allocate(previous%next)       !new node space
   previous%next%v = obj     !new object inserted
   previous%next%next => current !new next pointer
   
  end subroutine node_insert

  function empty_list(links)  result(t_or_f)
   type(link_list), intent(in) :: links
   logical                     :: t_or_f
   t_or_f = .not. associated (links%first%next)
  end function empty_list
  
  function new_list () result (OBJ)
   type(link_list) :: OBJ
   integer :: status
   allocate ( OBJ%first, stat=status)
   if(status/=0) CALL FATAL_ERROR("LinkList: Could not allocate new linklist")
   nullify(OBJ%first%next)
  end function new_list

  subroutine print_list (links)
   type(link_list), intent(in) :: links
   type(link_node), pointer    :: current
   logical                     :: headprint
   integer                     :: count
   current => links%first%next
   headprint = .true.
   count = 0
   do 
    if(.not. associated(current) ) exit 
    call screen_write(current%v,headprint)
    current => current%next
    headprint = .false.
    count = count +1
   end do
   write(ipt,*) "! PROC:",myid,"; # of local particles:", count
  end subroutine print_list

  subroutine print_data (links)
   type(link_list), intent(in) :: links
   type(link_node), pointer    :: current
   current => links%first%next
   do
    if(.not. associated(current) ) exit
    call particle_print(current%v)
    current => current%next
   end do
  end subroutine print_data

!!$  subroutine update_pathlength (links)
!!$   type(link_list), intent(in) :: links
!!$   type(link_node), pointer    :: current
!!$   current => links%first%next
!!$   do
!!$    if(.not. associated(current) ) exit
!!$    call set_pathlength(current%v)
!!$    current => current%next
!!$   end do
!!$  end subroutine update_pathlength


  subroutine print_id_list (links)
   type(link_list), intent(in) :: links
   type(link_node), pointer    :: current
   current => links%first%next
   do 
    if(.not. associated(current) ) exit 
    write(*,*)current%v%id 
    current => current%next
   end do
  end subroutine print_id_list

  subroutine shift_pos_list (links)
   type(link_list), intent(in) :: links
   type(link_node), pointer    :: current
   logical                     :: headprint
   current => links%first%next
   headprint = .true.
   do
    if(.not. associated(current) ) exit
    call shift_pos(current%v)
    current => current%next
   end do
  end subroutine shift_pos_list 
 

  function listsize (links) result (counter)
    type(link_list), intent(in) :: links
    type(link_node), pointer    :: current
    integer                     :: counter 
    counter = 0
    current => links%first%next
    
    do 
     if(.not. associated(current) ) exit 
     counter = counter + 1
     current => current%next
    end do
  end function listsize

  subroutine set_not_found (links) 
    type(link_list), intent(in) :: links
    type(link_node), pointer    :: current
    
   current => links%first%next
   do
    if(.not. associated(current) ) exit
    current%v%found = .false. 
    current => current%next
   end do
  end subroutine set_not_found 
 
  
end module linked_list
 
      
