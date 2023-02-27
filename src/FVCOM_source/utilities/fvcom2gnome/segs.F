Module Segment_Class
!  Use Gparms
  Use Mod_Prec
  Implicit None
  
  Type Edge_Type
    integer  :: v(2)
    integer  :: etype
    integer  :: ecntr
    integer  :: e(2)    
    integer  :: elem
    logical  :: used = .false.
    real(dp) :: length
    real(dp) :: x(2)
    real(dp) :: y(2)
    real(dp) :: h(2)
    real(dp) :: xc
    real(dp) :: yc
    integer  :: lney = 0
    integer  :: rney = 0
  End Type Edge_Type

  Type Segment_Type 
    integer          :: n_edges
    integer          :: beg_point
    integer          :: end_point
    integer          :: stype
    logical          :: closed = .false.
    real(dp)         :: length
    integer, pointer :: e(:)  
    integer          :: extra_pts = 0
    real(dp),pointer :: xextra(:)
    real(dp),pointer :: yextra(:)
  End Type Segment_Type


  contains

  integer function countfree(n,edgelist) result(sumfree)
  implicit none
  integer                      , intent(in) :: n 
  type(edge_type), dimension(n), intent(in) :: edgelist
  integer :: i

  sumfree = 0
  do i=1,n
    if(.not. edgelist(i)%used) sumfree = sumfree + 1
  end do
  
  end function countfree

  !---------------------------------------------
  !reorder edges so that node1 to node2 goes
  !ccw around the interior cell
  !---------------------------------------------
  subroutine order_edge(e)
  type(edge_type), intent(inout) :: e
  real(dp) :: dx1,dx2,dy1,dy2,prod,ftmp
  integer  :: itmp

  !cross product of v1-->v2  X v1-->(xc,yc)
  dx1 = e%x(2)-e%x(1)
  dx2 = e%xc  -e%x(1)
  dy1 = e%y(2)-e%y(1)
  dy2 = e%yc  -e%y(1)

  prod = dx1*dy2 - dx2*dy1
 
  !reverse edge node order if necessary
  if(prod < 0)then
    itmp   = e%v(1)
    e%v(1) = e%v(2)
    e%v(2) = itmp

    ftmp   = e%x(1)
    e%x(1) = e%x(2)
    e%x(2) = ftmp

    ftmp   = e%y(1)
    e%y(1) = e%y(2)
    e%y(2) = ftmp

    ftmp   = e%h(1)
    e%h(1) = e%h(2)
    e%h(2) = ftmp
  endif

  end subroutine order_edge

  subroutine edge_print(ee)
  type(edge_type), intent(in) :: ee
  write(*,*)'edge info'
  write(*,*)'node 1: ',ee%v(1)
  write(*,*)'coords: ',ee%x(1),ee%y(1),ee%h(1)
  write(*,*)'node 2: ',ee%v(2)
  write(*,*)'coords: ',ee%x(2),ee%y(2),ee%h(2)
  write(*,*)'lney: ',ee%lney,' rney ',ee%rney
  write(*,*)'elem: ',ee%elem
  end subroutine edge_print
    

  Subroutine Set_Neighbors(n,elist)
  Implicit None

  Integer        ,  intent(in)  :: n
  Type(edge_type),  intent(inout)  :: elist(n)
  Logical, allocatable :: mark(:)
  Integer  :: i

  allocate(mark(n)) ; mark = .false.

  
  do i=1,n
    if(.not.elist(i)%used)then
      elist(i)%used = .true.
      call leftney(n,elist,i)
      call rightney(n,elist,i)
    endif
  end do


  End Subroutine 
 
  recursive subroutine leftney(n,elist,istart) 
  Implicit None
  Integer        ,  intent(in)  :: n
  Type(edge_type),  intent(inout)  :: elist(n)
  Integer,          intent(in)  :: istart
  Integer :: i,v1,v2

  v1 = elist(istart)%v(1)
 
  do i=1,n
    v2 = elist(i)%v(2)
    if(v1 == v2 .and. .not. elist(i)%used) then
      elist(istart)%lney = i
      elist(i)%rney = istart
      call leftney(n,elist,i)
      elist(i)%used = .true.
    endif
  end do

  End Subroutine leftney


  recursive subroutine rightney(n,elist,istart)
  Implicit None
  Integer        ,  intent(in)  :: n
  Type(edge_type),  intent(inout)  :: elist(n)
  Integer,          intent(in)  :: istart
  Integer :: i,v1,v2

  v1 = elist(istart)%v(2)

  do i=1,n
    v2 = elist(i)%v(1)
    if(v1 == v2 .and. .not. elist(i)%used) then
      elist(istart)%lney = i
      elist(i)%rney = istart
      call rightney(n,elist,i)
      elist(i)%used = .true.
    endif
  end do

  End Subroutine rightney


End Module Segment_Class
