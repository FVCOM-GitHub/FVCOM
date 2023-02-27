!=======================================================================
! OSCAR Input Utility Module         
!
! Copyright:    2005(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED VERSIONS
!
! Authors:      G. Cowles 
!               School for Marine Science and Technology     
!               University of Massachusetts-Dartmouth
!
! Comments:     OSCAR Input Module  
!=======================================================================

Module Input_Util

Use Mod_Prec 
Use Mod_Utils
Implicit None

character(len=1), parameter :: com_char = '!'



!----------------------------------------
! interface definition: Overload Get_Val
!----------------------------------------
Interface Get_Val
  Module Procedure Get_Float_Val
  Module Procedure Get_String_Val
  Module Procedure Get_Logical_Val
  Module Procedure Get_Integer_Val
End Interface

Interface Get_Val_Array
  Module Procedure Get_Float_Val_Array
  Module Procedure Get_string_Val_Array
  Module Procedure Get_Logical_Val_Array
  Module Procedure Get_Integer_Val_Array
End Interface

!----------------------------------------
! functions    
!----------------------------------------

contains
  
  !---------------------------------------------------------------------
  !Get_Integer_Val -> Return Integer Value                    
  !---------------------------------------------------------------------
  Subroutine Get_Integer_Val(ival,filename,argname,line,echo,start,finish)  
  implicit none
  integer,          intent(out) :: ival 
  character(len=*), intent(in ) :: filename
  character(len=*), intent(in ) :: argname
  integer,          intent(out) :: line
  logical, optional,intent(in)  :: echo
  integer, optional,intent(in)  :: start
  integer, optional,intent(in)  :: finish 

  character(len=80) :: argval
  character(len=80) :: argtype
  integer           :: ierr

  !make sure file exists
  if(.not.check_exist(filename))then
    call error('Get_Val','halt','error reading: '//trim(argname), &
              'file: '//trim(filename)//' does not exist') 
  endif

  !check start and finish lines
  if(present(start) .and. present(finish))then
    if(finish < start)then
    call error('Get_Val','halt','finish must be greater than or equal to start')
    endif
  endif

  !parse file for argument line 
  call Extract_Val_String(filename,argname,argval,ierr,line,start,finish)
  if(ierr == 1)then
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
              'does not exist in file: '//trim(filename)) 
  endif

  !check argument type
  call Check_Arg_Type(argval,argtype)
  if(trim(argtype) /= 'integer')then
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
              'in file: '//trim(filename), &        
              'should be integer but is: '//trim(argtype))
  endif

  ! set return value
  read(argval,*)ival   

  ! echo to screen
  if(present(echo))then
    if(echo)then 
      write(*,'(A20,I10)')trim(argname)//': ',ival
    endif
  endif
  
  End Subroutine Get_Integer_Val


  !---------------------------------------------------------------------
  !Get_Float_Val -> Return Float
  !---------------------------------------------------------------------
  Subroutine Get_Float_Val(fval,filename,argname,line,echo,start,finish)  
  implicit none
  real(sp),         intent(out) :: fval
  character(len=*), intent(in ) :: filename
  character(len=*), intent(in ) :: argname
  integer,          intent(out) :: line
  logical, optional,intent(in)  :: echo
  integer, optional,intent(in)  :: start
  integer, optional,intent(in)  :: finish 

  character(len=80) :: argval
  character(len=80) :: argtype
  integer           :: ierr


  !make sure file exists
  if(.not.check_exist(filename))then
    call error('Get_Val','halt','error reading: '//trim(argname), &
              'file: '//trim(filename)//' does not exist') 
  endif

  !check start and finish lines
  if(present(start) .and. present(finish))then
    if(finish < start)then
    call error('Get_Val','halt','finish must be greater than or equal to start')
    endif
  endif

  !parse file for argument line 
  call Extract_Val_String(filename,argname,argval,ierr,line,start,finish)
  if(ierr == 1)then
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
              'does not exist in file: '//trim(filename)) 
  endif

  !check argument type
  call Check_Arg_Type(argval,argtype)
  if(trim(argtype) /= 'float')then
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
              'in file: '//trim(filename), &        
              'should be float but is: '//trim(argtype))
  endif

  ! set return value
  read(argval,*)fval  

  ! echo to screen 
  if(present(echo))then
    if(echo)then 
      write(*,'(A20,F10.4)')trim(argname)//': ',fval
    endif
  end if
    
  
  End Subroutine Get_Float_Val


  !---------------------------------------------------------------------
  !Get_Logical_Val -> Return Logical 
  !---------------------------------------------------------------------
  Subroutine Get_Logical_Val(lval,filename,argname,line,echo,start,finish)  
  implicit none
  logical,          intent(out) :: lval
  character(len=*), intent(in ) :: filename
  character(len=*), intent(in ) :: argname
  integer,          intent(out) :: line
  logical, optional,intent(in)  :: echo
  integer, optional,intent(in)  :: start
  integer, optional,intent(in)  :: finish 

  character(len=80) :: argval
  character(len=80) :: argtype
  integer           :: ierr

  !make sure file exists
  if(.not.check_exist(filename))then
    call error('Get_Val','halt','error reading: '//trim(argname), &
              'file: '//trim(filename)//' does not exist') 
  endif

  !check start and finish lines
  if(present(start) .and. present(finish))then
    if(finish < start)then
    call error('Get_Val','halt','finish must be greater than or equal to start')
    endif
  endif

  !parse file for argument line 
  call Extract_Val_String(filename,argname,argval,ierr,line,start,finish)
  if(ierr == 1)then
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
              'does not exist in file: '//trim(filename)) 
  endif

  !check argument type
  call Check_Arg_Type(argval,argtype)
  if(trim(argtype) /= 'logical')then
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
              'in file: '//trim(filename), &        
              'should be logical but is: '//trim(argtype))
  endif

  ! set return value
  lval = .false.
  if(argval(1:1) == 'T') lval = .true.

  ! echo to screen 
  if(present(echo))then
    if(echo)then 
      write(*,'(A20,L10)')trim(argname)//': ',lval  
    endif
  endif
  
  End Subroutine Get_Logical_Val


  !---------------------------------------------------------------------
  !Get_String_Val -> Return String 
  !---------------------------------------------------------------------
  Subroutine Get_String_Val(sval,filename,argname,line,echo,start,finish)  
  implicit none
  character(len=*), intent(out) :: sval
  character(len=*), intent(in ) :: filename
  character(len=*), intent(in ) :: argname
  integer,          intent(out) :: line
  logical, optional,intent(in)  :: echo
  integer, optional,intent(in)  :: start
  integer, optional,intent(in)  :: finish 

  character(len=80) :: argval
  character(len=80) :: argtype
  integer           :: ierr

  !make sure file exists
  if(.not.check_exist(filename))then
    call error('Get_Val','halt','error reading: '//trim(argname), &
              'file: '//trim(filename)//' does not exist') 
  endif

  !check start and finish lines
  if(present(start) .and. present(finish))then
    if(finish < start)then
    call error('Get_Val','halt','finish must be greater than or equal to start')
    endif
  endif

  !parse file for argument line 
  call Extract_Val_String(filename,argname,argval,ierr,line,start,finish)
  if(ierr == 1)then
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
              'does not exist in file: '//trim(filename)) 
  endif

  !check argument type
  call Check_Arg_Type(argval,argtype)
  if(trim(argtype) /= 'string')then
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
              'in file: '//trim(filename), &        
              'should be string but is: '//trim(argtype))
  endif

  ! set return value
  sval = adjustl(trim(argval))
  
  ! echo to screen 
  if(present(echo))then
    if(echo)then 
      write(*,'(A20,A20)') trim(argname)//': ',trim(sval)
    endif
  endif
  End Subroutine Get_String_Val


  !---------------------------------------------------------------------
  !Get_Float_Val_Array -> Return Float Array
  !---------------------------------------------------------------------
  Subroutine Get_Float_Val_Array(fval,filename,argname,size,echo)  
  implicit none

  character(len=*), intent(in ) :: filename
  character(len=*), intent(in ) :: argname
  integer,intent(in)            :: SIZE
  real(sp),dimension(size), intent(out) :: fval
  logical, optional,intent(in)  :: echo
 
 

  character(len=80) :: argval
  character(len=80) :: argtype
  integer           :: ierr
  integer           :: iscan
  integer           :: ntemp
  real(sp)          :: ftemp
  real(sp)          :: REALVEC(150)
  CHARACTER(LEN=30) :: FFMT

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start:Get_Float_Val_Array "

  !make sure file exists
  if(.not.check_exist(filename))then
    call error('Get_Val','halt','error reading: '//trim(argname), &
              'file: '//trim(filename)//' does not exist') 
  endif

  !read in 
  ISCAN = SCAN_FILE3(filename,argname,FVEC = REALVEC,NSZE = NTEMP)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING: ',trim(argname),': ',ISCAN
     call error('Get_Val','halt','error reading variable: '//trim(argname), &
        'does not exist in file: '//trim(filename)) 
  END IF

  IF(NTEMP > SIZE)THEN
    WRITE(*,*)'EXPECTED SIZE:',SIZE
    WRITE(*,*)'AND READ: ',REALVEC(1:ntemp)
    call error('Get_Val','warning','warning reading variable: '//trim(argname), &
        'number of specified size in : '//trim(filename)//' in not eq&
        &ual to expected size,specified size is used in the model  ')
  ELSEIF(NTEMP < SIZE)THEN
    WRITE(*,*)'EXPECTED SIZE:',SIZE
    WRITE(*,*)'AND READ: ',REALVEC
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
        'number of specified size in : '//trim(filename)//' in not eq&
        &ual to expected size ')
  END IF
  
  FVAL(1:SIZE)= REALVEC(1:SIZE)

  ! echo to screen 
  if(present(echo))then
    if(echo)then 
      WRITE(FFMT,'("(A20,",I0,"F10.4)")') SIZE
      write(*,FFMT) trim(argname)//': ',fval(1:SIZE)
!!$      write(*,'(A20,<size>F10.4)')trim(argname)//': ',fval(1:SIZE)
    endif
  end if
    
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Get_Float_Val_Array "  

  End Subroutine Get_Float_Val_Array





  !---------------------------------------------------------------------
  !Get_Float_Val_Array -> Return Float Array
  !---------------------------------------------------------------------
  Subroutine Get_Integer_Val_Array(ival,filename,argname,size,echo)  
  implicit none

  character(len=*), intent(in ) :: filename
  character(len=*), intent(in ) :: argname
  integer,intent(in)            :: SIZE
  integer,dimension(size), intent(out) :: ival
  logical, optional,intent(in)  :: echo
 
 

  character(len=80) :: argval
  character(len=80) :: argtype
  integer           :: ierr
  integer           :: iscan
  integer           :: ntemp
  real(sp)          :: ftemp
  real(sp)          :: REALVEC(150)
  CHARACTER(LEN=30) :: FFMT

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Get_Integer_Val_Array"

  !make sure file exists
  if(.not.check_exist(filename))then
    call error('Get_Val','halt','error reading: '//trim(argname), &
              'file: '//trim(filename)//' does not exist') 
  endif

  !read in 
  ISCAN = SCAN_FILE3(filename,argname,FVEC = REALVEC,NSZE = NTEMP)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING',trim(argname),': ',ISCAN
     call error('Get_Val','halt','error reading variable: '//trim(argname), &
        'does not exist in file: '//trim(filename)) 
  END IF
  print*,ntemp, size
  IF(NTEMP > SIZE)THEN
    WRITE(*,*)'EXPECTED SIZE:',SIZE
    WRITE(*,*)'AND READ: ',REALVEC(1:ntemp)
    call error('Get_Val','warning','warning reading variable: '//trim(argname), &
        'number of specified size in : '//trim(filename)//' in not eq&
        &ual to expected size,specified size is used in the model  ')
  ELSEIF(NTEMP < SIZE)THEN
    WRITE(*,*)'EXPECTED SIZE:',SIZE
    WRITE(*,*)'AND READ: ',REALVEC
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
        'number of specified size in : '//trim(filename)//' in not eq&
        &ual to expected size ')
  END IF

  iVAL(1:SIZE)= REALVEC(1:SIZE)

  ! echo to screen 
  if(present(echo))then
    if(echo)then 
      WRITE(FFMT,'("(A20,",I0,"I10)")') SIZE
      write(*,FFMT) trim(argname)//': ',ival(1:SIZE)
!!$      write(*,'(A20,<size>I10)')trim(argname)//': ',ival(1:SIZE)
    endif
  end if

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Get_Integer_Val_Array"    

  End Subroutine Get_Integer_Val_Array




  !---------------------------------------------------------------------
  !Get_Float_Val_Array -> Return Float Array
  !---------------------------------------------------------------------
  Subroutine Get_Logical_Val_Array(cval,filename,argname,size,echo)  
  implicit none

  character(len=*), intent(in ) :: filename
  character(len=*), intent(in ) :: argname
  integer,intent(in)            :: SIZE
  logical,dimension(size), intent(out) :: cval
  logical, optional,intent(in)  :: echo
 
 

  character(len=80) :: argval
  character(len=80) :: argtype
  integer           :: ierr,i
  integer           :: iscan
  integer           :: ntemp
  real(sp)          :: ftemp
  character(len=80), dimension(100):: CHARVEC
  CHARACTER(LEN=30) :: FFMT

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Get_Logical_Val_Array" 

  !make sure file exists
  if(.not.check_exist(filename))then
    call error('Get_Val','halt','error reading: '//trim(argname), &
              'file: '//trim(filename)//' does not exist') 
  endif

  !read in 
  ISCAN = SCAN_FILE3(filename,argname,CVEC = CHARVEC,NSZE = NTEMP)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING:',trim(argname),': ',ISCAN
     call error('Get_Val','halt','error reading variable: '//trim(argname), &
        'does not exist in file: '//trim(filename)) 
  END IF
 
  IF(NTEMP > SIZE)THEN
    WRITE(*,*)'EXPECTED SIZE:',SIZE
    WRITE(*,*)'AND READ: ',CHARVEC(1:ntemp)
    call error('Get_Val','warning','warning reading variable: '//trim(argname), &
        'number of specified size in : '//trim(filename)//' in not eq&
        &ual to expected size,specified size is used in the model  ')
  ELSEIF(NTEMP < SIZE)THEN
    WRITE(*,*)'EXPECTED SIZE:',SIZE
    WRITE(*,*)'AND READ: ',CHARVEC
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
        'number of specified size in : '//trim(filename)//' in not eq&
        &ual to expected size ')
  END IF

  do i=1,size
     IF((trim(CHARVEC(i))== "T" .OR. trim(CHARVEC(i))=="F").and.(len(trim(CHARVEC(i)))==1))THEN
        if(trim(CHARVEC(i)) == "T")cVAL(i)=.true.
        if(trim(CHARVEC(i)) == "F")cVAL(i)=.false.
     end if
  end do

  ! echo to screen 
  if(present(echo))then
    if(echo)then 
      WRITE(FFMT,'("(A20,",I0,"L10)")') SIZE
      write(*,FFMT) trim(argname)//': ',cval(1:SIZE)
!!$      write(*,'(A20,<size>L10)')trim(argname)//': ',cval(1:SIZE)
    endif
  end if
    
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Get_Logical_Val_Array"   

  End Subroutine Get_Logical_Val_Array








  !---------------------------------------------------------------------
  !Get_String_Val_Array -> Return String Array
  !---------------------------------------------------------------------
  Subroutine Get_String_Val_Array(sval,filename,argname,size,echo)  
  implicit none

  character(len=*), intent(in ) :: filename
  character(len=*), intent(in ) :: argname
  integer,intent(in)            :: SIZE
  character(len=*),dimension(size), intent(out) :: sval
  logical, optional,intent(in)  :: echo
 
 

  character(len=80) :: argval
  character(len=80) :: argtype
  integer           :: ierr,i
  integer           :: iscan
  integer           :: ntemp
  real(sp)          :: ftemp
  character(len=80), dimension(100):: CHARVEC
  CHARACTER(LEN=30) :: FFMT

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Get_String_Val_Array" 

  !make sure file exists
  if(.not.check_exist(filename))then
    call error('Get_Val','halt','error reading: '//trim(argname), &
              'file: '//trim(filename)//' does not exist') 
  endif

  !read in 
  ISCAN = SCAN_FILE3(filename,argname,CVEC = CHARVEC,NSZE = NTEMP)
  IF(ISCAN /= 0)THEN
     WRITE(IPT,*)'ERROR READING:',trim(argname),': ',ISCAN
     call error('Get_Val','halt','error reading variable: '//trim(argname), &
        'does not exist in file: '//trim(filename)) 
  END IF
 
  IF(NTEMP > SIZE)THEN
    WRITE(*,*)'EXPECTED SIZE:',SIZE
    WRITE(*,*)'AND READ: ',CHARVEC(1:ntemp)
    call error('Get_Val','warning','warning reading variable: '//trim(argname), &
        'number of specified size in : '//trim(filename)//' in not eq&
        &ual to expected size,specified size is used in the model  ')
  ELSEIF(NTEMP < SIZE)THEN
    WRITE(*,*)'EXPECTED SIZE:',SIZE
    WRITE(*,*)'AND READ: ',CHARVEC
    call error('Get_Val','halt','error reading variable: '//trim(argname), &
        'number of specified size in : '//trim(filename)//' in not eq&
        &ual to expected size ')
  END IF


  sVAL(1:SIZE)= CHARVEC(1:SIZE)

  ! echo to screen 
  if(present(echo))then
    if(echo)then 
      WRITE(FFMT,'("(A20,",I0,"A10)")') SIZE
      write(*,FFMT) trim(argname)//': ',sval(1:SIZE)
!!$      write(*,'(A20,<size>A10)')trim(argname)//': ',sval(1:SIZE)
    endif
  end if
    
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Get_String_Val_Array"   

  End Subroutine Get_String_Val_Array



  !---------------------------------------------------------------------
  !Extract_Val_String: 
  !   Parse Input File for Value String Assigned to Desired Variable
  !---------------------------------------------------------------------
  Subroutine Extract_Val_String(fname,argname,argval,ierr,line,start,finish)
  implicit none
  character(len=* ), intent(in ) :: fname 
  character(len=* ), intent(in ) :: argname 
  character(len=80), intent(out) :: argval 
  integer,           intent(out) :: ierr
  integer,           intent(out) :: line 
  integer, optional, intent(in ) :: start
  integer, optional, intent(in ) :: finish 
  !------------------------------------------
  character(len=80)  :: temp,temp2,text_line,inpline
  integer            :: locex,eqloc,length,i,linecnt

  ierr   = 0
  open(10,file=trim(fname)) ; rewind(10) 

  linecnt = 0
  !move forward to start line if present
  if(present(start))then
    do i=1,start-1
      read(10,*)
      linecnt = linecnt + 1
    end do
  endif

   do while(.true.)
     read(10,'(A)',end=20) inpline 
     linecnt = linecnt + 1

    ! determine location of comment character and remove blanks
    locex = index(inpline,com_char)
    if(locex == 0)locex = len_trim(inpline)+1
    temp = inpline(1:locex-1)
    text_line = temp

    ! set line length
    length = len_trim(text_line)
   
    ! ensure "=" exist and determine location
    eqloc = index(text_line,"=")


    ! split off variable name and value strings
    temp2 = text_line(1:eqloc-1)
    if(adjustl(trim(temp2)) == trim(argname))then
      argval  = adjustl(text_line(eqloc+1:length))
      line    = linecnt
      close(10)
      return
    endif

    ! exit if finish is present
    if(present(finish))then
      if(linecnt==finish)exit
    endif
  end do

  ! not found, return error = 1
  20 continue
  ierr = 1
  close(10)

  End Subroutine Extract_Val_String  

  !---------------------------------------------------------------------
  !check argument type (string,float,integer,logical)         
  !---------------------------------------------------------------------
  Subroutine Check_Arg_Type(argval,argtype)
  implicit none
  character(len=80), intent(out) :: argtype 
  character(len=80), intent(in ) :: argval  

  character(len=16) :: numchars 
  integer           :: largval
  integer           :: dotloc
  integer           :: i
  numchars = "0123456789+-Ee. " 

  ! argument length
  largval = len_trim(adjustl(argval))

  ! check if logical
  if(largval == 1)then
  if(argval(1:1) == "T" .or. argval(1:1) == "F")then 
    argtype = 'logical' 
    return
  endif
  endif
    
  ! check if it is a string  (contains characters other than 0-9,+,-,E,e,.)
  do i=1,largval
    if(index(numchars,argval(i:i)) == 0)then
      argtype = "string" 
      return
    endif
  end do

  ! check if it is a float
  dotloc = index(argval,".")
  if(dotloc /= 0) then
    argtype = "float"
  else
    argtype = "integer"
  end if

  End Subroutine Check_Arg_Type


  !-------------------------------------------------------
  !Error
  !Write Error Message to Screen
  !-------------------------------------------------------
  Subroutine Error(subname,errtype,err1,err2,err3)
    character(len=*)           :: subname
    character(len=*)           :: errtype
    character(len=*)           :: err1    
    character(len=*), optional :: err2 
    character(len=*), optional :: err3 

    if(errtype == 'halt')then
      write(*,*)
      write(*,*)'FATAL ERROR in subroutine: ',trim(subname)
      write(*,*)trim(err1)
      if(present(err2)) write(*,*)trim(err2)
      if(present(err2)) write(*,*)trim(err3)
      stop
    else if(errtype == 'warning')then
      write(*,*)
      write(*,*)'WARNING from subroutine: ',trim(subname)
      write(*,*)trim(err1)
      if(present(err2)) write(*,*)trim(err2)
      if(present(err2)) write(*,*)trim(err3)
    else
      write(*,*)'incorrect usage of subroutine Error'
      write(*,*)'called from subroutine: ',trim(subname)
      write(*,*)'errtype must be halt or warning'
      stop
    endif

  End Subroutine Error

  !-------------------------------------------------------
  !Check_Exist
  !Check for Existence of a File
  !-------------------------------------------------------
  Function Check_Exist(filein) Result(fexist) 
    implicit none
    character(len=*)           :: filein 
    logical fexist 
    inquire(exist=fexist,file=filein)
  End Function Check_Exist  

  !-------------------------------------------------------
  !Check_Alloc
  !Check For Successfull Allocation of Data
  !-------------------------------------------------------
  Subroutine Check_Alloc(errorin,dsize,varname,subname) 
    implicit none
    integer,          intent(in) :: errorin
    integer,          intent(in) :: dsize 
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: subname

    if(errorin /= 0)then
      write(*,*)'insufficient space to allocate array: '//trim(varname)
      write(*,*)'in subroutine: '//trim(subname)
      write(*,*)'tried to allocate ',dsize,' units'
      stop
    endif
   
  End Subroutine Check_Alloc  


  !-------------------------------------------------------
  !Check_IOError
  !Check For Successfull R/W 
  !-------------------------------------------------------
  Subroutine Check_IOerror(errorin,filename,varname) 
    implicit none
    integer,          intent(in) :: errorin
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname

    if(errorin /= 0)then
      write(*,*)'failed read/write of: '//trim(varname)  
      write(*,*)'from file: '//trim(filename)    
      stop
    endif
   
  End Subroutine Check_IOerror

End Module Input_Util
