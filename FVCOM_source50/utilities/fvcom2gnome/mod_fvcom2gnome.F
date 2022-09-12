Module mod_fvcom2gnome
  use mod_nctools
  use mod_ncll
  use control,   only : ipt, msr, par, serial
  use mod_utils, only : pshutdown, fatal_error
  use segment_class
  implicit none


  !constants
!  real(sp) :: zero  = 0.0
  real(sp) :: ahalf = 0.5
  real(sp) :: one   = 0.0

  character(len=80) :: inputfile  !no default, required
  character(len=80) :: bndryfile  !no default, required
  character(len=80) :: outputfile = 'gnome.nc'
  character(len=80) :: tstring 
  character(len=80) :: basedate 
  character(len=80) :: user_date 
  logical :: get_time_from_model = .true.
  logical :: have_bndry_file = .true.
  logical :: with_basedate = .false.
  real(sp):: t0   = 0.0
  integer :: fbeg = 1
  integer :: fend = huge(fbeg) 
  integer :: fint = 1
  integer :: fvcom_dim = 3

  logical :: bndryfile_arg = .false.
  logical :: inputfile_arg = .false.


  !problem dimensions
  integer :: N_frames
  integer :: N_frames_out
  integer :: N_elems
  integer :: N_verts
  integer :: N_bsegs
  integer :: N_siglay
  integer :: N_siglev
  integer :: N_MaxNode
  integer :: N_MaxElem

  !mesh variables
  real(sp), allocatable, target :: x(:) 
  real(sp), allocatable, target :: y(:) 
  real(sp), allocatable, target :: lon(:) 
  real(sp), allocatable, target :: lat(:) 
  integer , allocatable :: nv(:,:)
  real(sp), allocatable, target :: art1(:)
  real(sp), allocatable :: art(:)
  integer , allocatable :: ntve(:)
  integer , allocatable :: nbve(:,:)
  integer , allocatable :: ntsn(:)
  integer , allocatable :: nbsn(:,:)
  real(sp), allocatable :: sigma(:)
  real(sp), allocatable, target :: xc(:)
  real(sp), allocatable, target :: yc(:)
  integer , allocatable :: ntype(:)

  integer :: N_segs
  type(segment_type), allocatable :: seglist(:)
  !integer*2, allocatable :: GNOME_bndry(:,:) 
  integer, allocatable :: GNOME_bndry(:,:) 

  integer :: N_Edges
  type(edge_type),    allocatable :: edgelist(:)

  integer :: N_obc
  integer, allocatable :: obc_nodes(:)

  !time-dependent data from input
  real(sp), target, allocatable :: u(:,:)
  real(sp), target, allocatable :: v(:,:)
  real(sp), target, allocatable :: ua(:)
  real(sp), target, allocatable :: va(:)

  !time-dependent data for output
  real(sp), target :: gtime 
  real(sp), target, allocatable :: un(:)
  real(sp), target, allocatable :: vn(:)
  real(sp), target, allocatable :: junk(:)
  real(sp), target, allocatable :: z(:)

  !input netcdf vars
  type(ncfile), pointer :: nc_infile 
  type(ncdim),  pointer :: DIM_nele
  type(ncdim),  pointer :: DIM_node
  type(ncdim),  pointer :: DIM_three
  type(ncdim),  pointer :: DIM_four  
  type(ncdim),  pointer :: DIM_siglay
  type(ncdim),  pointer :: DIM_siglev
  type(ncdim),  pointer :: DIM_time  
  type(ncdim),  pointer :: DIM_MaxNode
  type(ncdim),  pointer :: DIM_MaxElem

  !output file netcdf vars
  type(ncfile), pointer :: nc_outfile
  type(ncdim),  pointer :: DIM_nele_out
  type(ncdim),  pointer :: DIM_node_out
  type(ncdim),  pointer :: DIM_nbnd_out  
  type(ncdim),  pointer :: DIM_nbi_out 
  type(ncdim),  pointer :: DIM_sig_out 
  type(ncdim),  pointer :: DIM_time_out  
  integer :: ofid
  integer :: node_did,nele_did,nbnd_did,nbi_did,time_did
  integer :: bnd_vid,time_vid,lon_vid,lat_vid,u_vid,v_vid,z_vid

  contains

  !=====================================================================
  ! Write Help Message to Screen for Command Line Usage of fvcom2gnome
  !=====================================================================
  Subroutine HelpMessage
    Implicit None
             
    write(IPT,*) 'fvcom2gnome Usage:'
    write(IPT,*) 'fvcom2gnome --inputfile=XXX --bndryfile=ZZZ --outputfile=YYY -fbeg=i -fend=j -fint=k'
    write(IPT,*) 'required:'
    write(IPT,*) '  inputfile:  FVCOM data file'
    write(IPT,*) '              XXX is the full path to the input file'
    write(IPT,*) '  bndryfile:  FVCOM open boundary node list file'
    write(IPT,*) '              ZZZ is the full path to the input file'
    write(IPT,*) '              if no open boundry, specify bndryfile=none'
    write(IPT,*) 'optional: '
    write(IPT,*) ' outputfile: GNOME-compatible netcdf file'
    write(IPT,*) '              default = gnome.nc'
    write(IPT,*) '              YYY is the full path to the output file'
    write(IPT,*) '  fbeg     :  beginning frame [positive integer] to pull from inputfile' 
    write(IPT,*) '              default = 1'
    write(IPT,*) '  fend     :  last frame [positive integer] to pull from inputfile' 
    write(IPT,*) '              default = last frame in file'
    write(IPT,*) '  fint     :  frame interval [positive integer] to pull from inputfile' 
    write(IPT,*) '              default = 1'
    write(IPT,*) '  tstring  :  string description for start time: e.g.: 2008-09-15 23:00:00 UTC' 
    write(IPT,*) '              this is the gregorian date for the first frame in FVCOM output'
    write(IPT,*) '              default = will calculate days since 1970-01-01 00:00:00 00:00'
    write(IPT,*) '                        from model time.' 
  End Subroutine HelpMessage
          

  !=====================================================================
  ! Parse commandline and set control variables for fvcom2gnome 
  !=====================================================================
  Subroutine parse_commandline(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)
    use mod_sng
    use mod_time

    character(len=*), INTENT(IN)::CVS_Id  ! [sng] CVS Identification
    character(len=*), INTENT(IN)::CVS_Date ! [sng] Date string
    character(len=*), INTENT(IN)::CVS_Name ! [sng] File name string
    character(len=*), INTENT(IN)::CVS_Revision ! [sng] File revision string

    
    character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
    ! Command-line parsing
    character(80)::arg_val ! [sng] command-line argument value
    character(200)::cmd_ln ! [sng] command-line
    character(80)::opt_sng ! [sng] Option string
    character(2)::dsh_key  ! [sng] command-line dash and switch
    character(200)::prg_ID ! [sng] Program ID

    type(ncfile), pointer :: ncf
    type(ncvar),  pointer :: var
    !integer :: itime,itime2
    logical :: FOUND
    type(time) :: start_time
    character(len=80) :: gregtime

    logical :: fexist
    integer ::arg_idx 
    integer ::arg_nbr 
    integer ::opt_lng 

    ! Main code
    call ftn_strini(cmd_ln) ! [sng] sng(1:len)=NUL

    call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
    call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID

    arg_nbr = command_argument_count() 

    !sanity on number of arguments
    if (arg_nbr == 0  .or. arg_nbr > 5) then
     if(msr) call HelpMessage
     call PSHUTDOWN
    end if

    arg_idx=1 ! initialize Counting index

    !=> main loop over command line arguments and parse
    do while (arg_idx <= arg_nbr)
      call ftn_getarg_wrp(arg_idx,arg_val) ! call getarg, increment arg_idx
      dsh_key=arg_val(1:2) 
      if (dsh_key == "--") then
        opt_lng=ftn_opt_lng_get(arg_val) ! Length of option
        if (opt_lng <= 0) then
           if(MSR) write(IPT,*) "Long option has no name"
           call PSHUTDOWN
        end if

        !grab the actual argument from the string
        opt_sng=arg_val(3:2+opt_lng) 

        select case(opt_sng)

        !---------------------------------------------
        case("inputfile") !set the inputfile 
        !---------------------------------------------
          call ftn_arg_get(arg_idx,arg_val,inputfile) 
          inputfile_arg = .true.
        !---------------------------------------------
        case("bndryfile") !set the bndryfile 
        !---------------------------------------------
          call ftn_arg_get(arg_idx,arg_val,bndryfile) 
          if(bndryfile(1:4) == "NONE" .or. bndryfile(1:4) == "none")then
            have_bndry_file = .false.
          endif
          bndryfile_arg = .true.
        !---------------------------------------------
        case("outputfile") !set the outputfile name
        !---------------------------------------------
          outputfile = ''
          call ftn_arg_get(arg_idx,arg_val,outputfile) 
        !---------------------------------------------
        case("fbeg")       !set the begin frame
        !---------------------------------------------
          call ftn_arg_get(arg_idx,arg_val,fbeg) 
        !---------------------------------------------
        case("fend")       !set the end frame
        !---------------------------------------------
          call ftn_arg_get(arg_idx,arg_val,fend) 
        !---------------------------------------------
        case("fint")       !set the frame interval 
        !---------------------------------------------
          call ftn_arg_get(arg_idx,arg_val,fint) 
        !---------------------------------------------
        case("tstring")    !set the time string 
        !---------------------------------------------
          get_time_from_model = .false.
          call ftn_arg_get(arg_idx,arg_val,user_date) 
        !---------------------------------------------
        case default       !option not valid
        !---------------------------------------------
          arg_idx=arg_idx-1 ! [idx] Counting index
          if(msr) call HelpMessage
          call PSHUTDOWN
  
        end select

        cycle !parse next arg (awkward)

       endif ! endif long option

       if (dsh_key == "-V" .or.dsh_key == "-v" ) then

          if(MSR) write(IPT,*) prg_id
          call PSHUTDOWN

       else if (dsh_key == "-H" .or.dsh_key == "-h" ) then
          if(MSR) call HelpMessage
          Call PSHUTDOWN

       else ! Option not recognized
          arg_idx=arg_idx-1 ! [idx] Counting index
          if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
       endif ! endif arg_val


    end do 
    !<= end main loop over command line arguments and parse

    !check the arguments
    if(.not.inputfile_arg)then  
       Call warning("YOU MUST SPECIFY AN INPUTFILE") 
       call HelpMessage
       call pshutdown
    endif
    if(.not.bndryfile_arg)then  
       Call warning("YOU MUST SPECIFY A BOUNDARYFILE") 
       call HelpMessage
       call pshutdown
    endif
 
    inquire(exist=fexist,file=trim(inputfile)) 
    if(.not.fexist)then
       Call fatal_error("THE INPUTFILE DOES NOT EXIST",trim(inputfile)) 
    endif
    if(have_bndry_file)then 
      inquire(exist=fexist,file=trim(inputfile)) 
      if(.not.fexist)then
         Call fatal_error("THE INPUTFILE DOES NOT EXIST",trim(inputfile)) 
      endif
    endif


    if(fbeg < 0)then 
       call fatal_error("Your begin frame is non-positive") 
    endif

    if(fend < 0)then 
       call fatal_error("Your end frame is non-positive") 
    endif

    if(fend < fbeg)then 
       call fatal_error("Your begin frame is after your end frame") 
    endif

    if(fint < 1)then
       call fatal_error("Your frame interval must be a positive number") 
    endif

    !report
    write(IPT,*)'====  command line summary ====='
    write(IPT,*)'inputfile:   ',trim(inputfile)
    write(IPT,*)'bndryfile:   ',trim(bndryfile)
    write(IPT,*)'outputfile:  ',trim(outputfile)
    write(IPT,*)'begin frame: ',fbeg
    write(IPT,*)'frame intvl: ',fint
    if(fend == huge(fbeg))then
      write(IPT,*)'end frame  : ',"last frame in file"
      fend = -1
    else
      write(IPT,*)'end frame  : ',fend
    endif

    !get model first time.
    ncf => new_file()
    ncf%fname = trim(inputfile)
    ncf%ftime=>new_ftime()
    ncf%ftime%next_stkcnt = 0
    call nc_open(ncf)
    call nc_load(ncf)
    start_time = get_file_time(ncf,fbeg)      
    nullify(ncf)
    T0 = 0.0 !float(start_time%mjd) + float(start_time%musod)/(3600.*24.*1e6)
    if(get_time_from_model)then
      T0 = 40587.  !MJD of Jan 1, 1970
      gregtime = write_datetime(start_time,6,"UTC")
      !!tstring = "days since "//gregtime(1:10)//" "//gregtime(12:19)//" UTC"
      tstring = "days since 1970-01-01 00:00:00 00:00"  
      with_basedate = .false. !.true.
    else
      tstring = "days since "//trim(user_date)
      with_basedate = .false.
    endif

    write(*,*)'gnome file time represents:'
    write(*,*)trim(tstring)
 



  End Subroutine parse_commandline

  !=====================================================================
  ! open the input file, allocate and connect variables 
  !=====================================================================
  Subroutine Setup_Input 
    logical :: found
    type(ncdim),  pointer :: dim
    type(ncfile), pointer :: ncf
    type(ncvar),  pointer :: var
    integer  :: i,itmp,iscan
    real(sp) :: xvts(3),yvts(3)
  
  
    !-----------------------------------------
    !set file pointer and make sure it exists 
    !-----------------------------------------
    nullify(nc_infile)
    ncf => new_file()
    ncf%fname = trim(inputfile)
    ncf%ftime=>new_ftime()
    ncf%ftime%next_stkcnt = 0
    call nc_open(ncf)
    call nc_load(ncf)
    nc_infile => ncf

    !-----------------------------------------
    !read problem dimensions from FVCOM file
    !-----------------------------------------

    !time
    dim => find_dim(nc_infile,'time',found)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN FILENAME",&
            & "FILE NAME: "//TRIM(inputfile),&
            &"COULD NOT FIND DIMENSION 'time'")

    N_frames = DIM%DIM

    !number of elements
    dim => find_dim(nc_infile,'nele',found)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN FILENAME",&
            & "FILE NAME: "//TRIM(inputfile),&
            &"COULD NOT FIND DIMENSION 'nele'")

    N_elems = DIM%DIM

    !number of vertices 
    dim => find_dim(nc_infile,'node',found)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN FILENAME",&
            & "FILE NAME: "//TRIM(inputfile),&
            &"COULD NOT FIND DIMENSION 'node'")

    N_verts = DIM%DIM

    !number of vertices 
    dim => find_dim(nc_infile,'siglay',found)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN FILENAME",&
            & "FILE NAME: "//TRIM(inputfile),&
            &"COULD NOT FIND DIMENSION 'siglay'")

    N_siglay = DIM%DIM

    !number of maximum vertices surrounding a node
    dim => find_dim(nc_infile,'maxnode',found)
       IF(.not. FOUND) then 
         CALL warning &
            & ("IN FILENAME",&
            & "FILE NAME: "//TRIM(inputfile),&
            &"COULD NOT FIND DIMENSION 'maxnode'")
        call fatal_error("need grid metrics in the FVCOM output file", &
            & 'Set "NC_GRID_METRICS = T" in your NETCDF namelist', &
            & "in your runtime control file and rerun") 
        ENDIF

    N_maxnode = DIM%DIM

    !number of maximum elements surrounding a node
    dim => find_dim(nc_infile,'maxelem',found)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN FILENAME",&
            & "FILE NAME: "//TRIM(inputfile),&
            &"COULD NOT FIND DIMENSION 'maxelem'")

    N_maxelem = DIM%DIM

    !-----------------------------------------
    !load mesh into memory
    !-----------------------------------------

    !x
    allocate(x(N_verts)) ; x = zero
    var => find_var(nc_infile,'x',found)  
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'x'")
    call nc_connect_avar(var,x)
    call nc_read_var(var)
    nullify(var)

    !y
    allocate(y(N_verts)) ; y = zero
    var => find_var(nc_infile,'y',found)  
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'y'")
    call nc_connect_avar(var,y)
    call nc_read_var(var)
    nullify(var)

    !xc
    allocate(xc(N_elems)) ; xc = zero
    var => find_var(nc_infile,'xc',found)
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'xc'")
    call nc_connect_avar(var,xc)
    call nc_read_var(var)
    nullify(var)

    !yc
    allocate(yc(N_elems)) ; yc = zero
    var => find_var(nc_infile,'yc',found)
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'yc'")
    call nc_connect_avar(var,yc)
    call nc_read_var(var)
    nullify(var)

    !longitude
    allocate(lon(N_verts)) ; lon = zero
    var => find_var(nc_infile,'lon',found)  
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'lon'")
    call nc_connect_avar(var,lon)
    call nc_read_var(var)
    nullify(var)

    !latitude  
    allocate(lat(N_verts)) ; lat = zero
    var => find_var(nc_infile,'lat',found)  
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'lat'")
    call nc_connect_avar(var,lat)
    call nc_read_var(var)
    nullify(var)

    !nv        
    allocate(nv(N_elems,3)) ; nv = zero
    var => find_var(nc_infile,'nv',found)  
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'nv'")
    call nc_connect_avar(var,nv)
    call nc_read_var(var)
    if(maxval(nv) /= N_verts) call fatal_error &
            & ("maxval(nv) /= N_verts") 
    nullify(var)

    !open boundary nodes
    !gwc debug, add to mesh metrics output 
    if(have_bndry_file)then
      open(unit=33,file= bndryfile)
      iscan = scan_file(33,"OBC Node Number",iscal = N_obc)
      if(iscan /= 0) then
         call fatal_error&
              &('Improper formatting of open_boundary_file', & 
              & 'The header must contain: "OBC Node Number = "', &
              & 'Followed by and integer number of boundary nodes')
      end if
      close(33)
      open(unit=33,file= bndryfile)
      allocate(obc_nodes( max(N_obc,1) )) ; obc_nodes = 0
      read(33,*)
      do i=1,N_obc
        read(33,*)itmp,obc_nodes(i)
      end do
      close(33)
    endif
 

    !nbve       
    allocate(nbve(N_verts,N_MaxElem)) ; nbve = zero
    var => find_var(nc_infile,'nbve',found)  
    if(.NOT. found)then
        call warning &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'nbve'")
        call fatal_error('Set "NC_GRID_METRICS = T" in your NETCDF namelist', &
            & "in your runtime control file") 
    endif
    call nc_connect_avar(var,nbve)
    call nc_read_var(var)
    nullify(var)

    !ntve       
    allocate(ntve(N_verts)) ; ntve = zero
    var => find_var(nc_infile,'ntve',found)  
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'ntve'")
    call nc_connect_avar(var,ntve)
    call nc_read_var(var)
    nullify(var)

   !nbsn       
    allocate(nbsn(N_verts,N_MaxNode)) ; nbsn = zero
    var => find_var(nc_infile,'nbsn',found)
    if(.NOT. found)then
        call warning &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'nbsn'")
        call fatal_error('Set "NC_GRID_METRICS = T" in your NETCDF namelist', &
            & "in your runtime control file")
    endif
    call nc_connect_avar(var,nbsn)
    call nc_read_var(var)
    nullify(var)

    !ntsn       
    allocate(ntsn(N_verts)) ; ntsn = zero
    var => find_var(nc_infile,'ntsn',found)
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'ntsn'")
    call nc_connect_avar(var,ntsn)
    call nc_read_var(var)
    nullify(var)

    !art1 (area of elements surrounding a node)
    allocate(art1(N_verts)) ; art1 = zero
    var => find_var(nc_infile,'art1',found)  
    if(.NOT. found) call fatal_error &
            & ("in filename",&
            & "file name: "//trim(inputfile),&
            &"could not find variable 'art1'")
    call nc_connect_avar(var,art1)
    call nc_read_var(var)
    nullify(var)

    !determine if output is 3D or 2D
    fvcom_dim = 3
    var => find_var(nc_infile,'u',found)
    if(.not.found)then
      fvcom_dim = 2
      var => find_var(nc_infile,'ua',found)
      if(.not.found)then
        if(.NOT. found) call fatal_error &
             & ("in filename",&
             & "file name: "//trim(inputfile),&
             &"could not find variable 'u' or 'ua'")
      endif
   endif
    
   !allocate dataspace
   if(fvcom_dim ==3)then  !3D
    allocate(u(N_elems,N_siglay)); u = zero
    allocate(v(N_elems,N_siglay)); v = zero
   else !2D
    allocate(ua(N_elems)); ua = zero
    allocate(va(N_elems)); va = zero
   endif


    !report
    write(ipt,*) '============ input dimensions ================='
    write(ipt,*) 'N_frames  :',N_frames
    write(ipt,*) 'N_elems   :',N_elems
    write(ipt,*) 'N_verts   :',N_verts
    write(ipt,*) 'N_siglay  :',N_frames
    write(ipt,*) 'N_MaxElem :',N_maxelem
    write(ipt,*) 'N_MaxNode :',N_maxnode
    write(ipt,*) 'N_obc     :',N_obc    
    if(fvcom_dim ==3)then
       write(ipt,*) 'ProbDim   : 3D - using surface u/v'
    else
       write(ipt,*) 'ProbDim   : 2D - using ua/va' 
    endif

    !set/reset begin and end frame
    if(fbeg > N_frames) call fatal_error &
            & ("fbeg exceeds the number of frames in the file"  ) 

    if(fend == -1)fend = N_frames

    if(fend > N_frames) call fatal_error &
            & ("fend exceeds the number of frames in the file"  ) 

    !determine the number of frames to output
    N_frames_out = 0
    do i=fbeg,fend,fint
      N_frames_out = N_frames_out + 1 
    end do
    write(ipt,*) '#OutFrames:',N_frames_out

    !calculate cell areas
    allocate(art(N_elems)) ; art = zero
    do i=1,N_elems
      xvts = x(nv(i,1:3))
      yvts = y(nv(i,1:3))
      art(i) = cell_area(xvts,yvts)
    end do

    !mark open boundary nodes
    allocate(ntype(N_verts)) ; ntype = 0
    do i=1,N_obc
      ntype(obc_nodes(i)) = 1
    end do

    
    !dump basic mesh data
    write(ipt,*) '============ basic mesh info  ================='
    write(ipt,*) 'Avg Longitude: ',sum(lon)/N_verts
    write(ipt,*) 'Avg Latitude : ',sum(lat)/N_verts
    write(ipt,*) 'Min Cell Area: ',minval(art)
    write(ipt,*) 'Max Cell Area: ',maxval(art)
    write(ipt,*) 'Min Edge Lnth: ',sqrt(2.*minval(art))
    


    
  End Subroutine Setup_Input

  !=====================================================================
  ! Setup the Boundary 
  !=====================================================================
  Subroutine Get_Edgelist 
    Integer :: e1,icell,n1,n2,ecnt,jnode,ee,jcell,i,j
    Type(edge_type),    allocatable :: etmp(:)
    Logical, allocatable :: eset(:)

    write(ipt,*) '========= setting up the boundary ============='

    !---------------------------------------------------------
    !set up a list of boundary edges 
    !---------------------------------------------------------
    allocate(etmp(N_elems*3)) ; etmp%ecntr = 0 ; etmp%etype = 0
    allocate(eset(N_verts))

    !set vertices on either side of edge (v(1),v(2))
    ecnt = 0
    do i=1,N_verts 
      do j=1,ntsn(i)
        jnode = nbsn(i,j)
        if(.not.eset(jnode))then
          ecnt = ecnt + 1
          etmp(ecnt)%v(1) = i
          etmp(ecnt)%v(2) = jnode
        endif
      end do
      eset(i) = .true.
    end do

    !set elements on either side of edge 
    do ee=1,ecnt
      etmp(ee)%e(1) = 0
      etmp(ee)%e(2) = 0
      n1 = etmp(ee)%v(1) ; n2 = etmp(ee)%v(2)
      do i=1,ntve(n1)
        icell = nbve(n1,i)
        do j=1,ntve(n2)
          jcell = nbve(n2,j)
          if(icell == jcell)then
            etmp(ee)%ecntr = etmp(ee)%ecntr + 1
            etmp(ee)%e( etmp(ee)%ecntr ) = icell
          end if
        end do
      end do
    end do

    !count and flag boundary edges
    e1 = 0
    do ee=1,ecnt
      etmp(ee)%etype = 0
      n1 = etmp(ee)%v(1)
      n2 = etmp(ee)%v(2)
      if(etmp(ee)%ecntr == 1) then 
        if (ntype(n1)+ntype(n2) == 2)then
          etmp(ee)%etype = 2   !edge is open boundary edge
        else
          etmp(ee)%etype = 1   !edge is solid boundary edge
        endif
        e1 = e1 + 1
      end if
    end do

    !extract a real list containing only boundary edges
    N_bsegs = e1
    write(ipt,*) '#Bndry Edges : ',N_bsegs          
    allocate(edgelist(N_bsegs))
  
    e1 = 0
    do ee=1,ecnt
      if(etmp(ee)%etype > 0)then
        e1 = e1 + 1
        edgelist(e1) = etmp(ee)
        if(edgelist(e1)%e(1) /= 0)edgelist(e1)%elem = edgelist(e1)%e(1)
        if(edgelist(e1)%e(2) /= 0)edgelist(e1)%elem = edgelist(e1)%e(2)
      endif
    end do
    deallocate(etmp)

    !set xc,yc, coordinates of neighbor cell center
    !this will be used to reorient edge later
    do e1=1,N_bsegs
      icell = edgelist(e1)%elem
      edgelist(e1)%xc = xc(icell)
      edgelist(e1)%yc = yc(icell)
      n1 = edgelist(e1)%v(1)
      n2 = edgelist(e1)%v(2)
      edgelist(e1)%x(1) = x(n1)
      edgelist(e1)%x(2) = x(n2)
      edgelist(e1)%y(1) = y(n1)
      edgelist(e1)%y(2) = y(n2)
    end do
  End Subroutine Get_Edgelist

  !=====================================================================
  ! Order the edges edge 2 edge 
  !=====================================================================
  Subroutine Order_Edges 
    Integer :: i
    do i=1,N_bsegs
      call order_edge(edgelist(i))
    end do
  End Subroutine Order_Edges
    
  !=====================================================================
  ! Determine Left and Right Neighbors 
  !=====================================================================
  Subroutine Set_Seg_Neighbors 
    Call set_neighbors(N_bsegs,edgelist)
  End Subroutine Set_Seg_Neighbors

  !=====================================================================
  ! Generate segments of closed loops of edges 
  !=====================================================================
  Subroutine Get_Seglist
    Integer :: i,n_open,efrst,elast,ii,j,ee
    Type(segment_type), allocatable :: ltmp(:)
    Integer, allocatable :: seg(:)


    !build segments by grabbing connected edges => ltmp
    allocate(ltmp(N_bsegs))
    do i=1,N_bsegs
      if(edgelist(i)%lney == 0)then
        N_segs = N_segs + 1
        !setup seg
        allocate(seg(N_bsegs)); seg = 0
        ii = 1
        seg(ii) = i
        do while(edgelist(seg(ii))%rney /= 0)
          ii = ii + 1
          seg(ii) = edgelist(seg(ii-1))%rney
        end do
        !transfer 
        allocate(ltmp(N_segs)%e(ii))
        ltmp(N_segs)%n_edges = ii
        ltmp(N_segs)%e(1:ii) = seg(1:ii)
  
        deallocate(seg)
      endif
    end do

    !transfer to true segment array
    allocate( seglist(N_segs) )
    do i=1,N_segs
      seglist(i) = ltmp(i)
    end do
    deallocate(ltmp)

    !see if segment is closed (begin point /= end point)
    !set segment closed logical to T
    do i=1,N_segs 
      efrst = seglist(i)%e(1)
      elast = seglist(i)%e(seglist(i)%N_edges)
      seglist(i)%beg_point = edgelist(efrst)%v(1)
      seglist(i)%end_point = edgelist(elast)%v(2)
      if(seglist(i)%beg_point == seglist(i)%end_point) seglist(i)%closed = .true.
    end do
    
    !count open segments
    n_open = 0
    do i=1,N_segs
      if(.not. seglist(i)%closed)then
        write(ipt,*)'open segment! :',i  
        write(ipt,*)'halting' 
        stop
        n_open = n_open + 1
      endif
    end do

    !determine which segment includes the open boundary
    do i=1,N_segs
      seglist(i)%stype = 0 
      do j=1,seglist(i)%n_edges 
        ee = seglist(i)%e(j)
        if(edgelist(ee)%etype ==2)then
          seglist(i)%stype = 1
        endif
      end do
    end do

    
  
    !report segment count
    write(ipt,*) '#Bndry Segs  : ',N_segs           
    do i=1,N_segs
      write(ipt,*)'seg info     : ',i,seglist(i)%n_edges,seglist(i)%stype
    end do


  End Subroutine Get_Seglist 

  !=====================================================================
  ! Generate the GNOME boundary array 
  !=====================================================================
  Subroutine Gen_GNOME_bndry 

    Integer :: cnt,i,j,ee,n1,n2

    !allocate the array
    allocate(GNOME_bndry(4,N_bsegs)) ; GNOME_bndry = 0

    cnt = 0

    !dump segment with open boundary first
    do i=1,N_segs
      if(seglist(i)%stype == 1)then
        do j=1,seglist(i)%n_edges 
          cnt = cnt + 1
          ee = seglist(i)%e(j)
          GNOME_bndry(1,cnt) = edgelist(ee)%v(1)
          GNOME_bndry(2,cnt) = edgelist(ee)%v(2) 
          GNOME_bndry(3,cnt) = i- 1 
          GNOME_bndry(4,cnt) = edgelist(ee)%etype -1
        end do
      endif
    end do

    !dump islands next 
    do i=1,N_segs
      if(seglist(i)%stype == 0)then
        do j=1,seglist(i)%n_edges
          cnt = cnt + 1
          ee = seglist(i)%e(j)
          GNOME_bndry(1,cnt) = edgelist(ee)%v(1)
          GNOME_bndry(2,cnt) = edgelist(ee)%v(2)
          GNOME_bndry(3,cnt) = i- 1
          GNOME_bndry(4,cnt) = edgelist(ee)%etype -1
        end do
      endif
    end do

    !dump to check
    open(unit=35,file='check.dat',form='formatted')
    do i=1,N_bsegs
      n1 = GNOME_bndry(1,i)
      n2 = GNOME_bndry(2,i)
      write(34,*)GNOME_bndry(1:4,i)
      write(35,'(4F14.7,2I5)')lon(n1),lon(n2),lat(n1),lat(n2),GNOME_bndry(3,i),GNOME_bndry(4,i)
    end do
    close(35)
    


  End Subroutine Gen_GNOME_bndry
   

  !=====================================================================
  ! open the output file type and setup dimensions/variables 
  !=====================================================================
  Subroutine Setup_Output
    use netcdf

    !create new file and set name
    call cfcheck( nf90_create(trim(outputfile),nf90_clobber,ofid) ) 

    !global attributes
    call cfcheck( nf90_put_att(ofid,nf90_global,"file_type", "FEM"))
    call cfcheck( nf90_put_att(ofid,nf90_global,"Conventions", "COARDS"))
    call cfcheck( nf90_put_att(ofid,nf90_global,"grid_type", "Triangular"))
    !call cfcheck( nf90_put_att(ofid,nf90_global,"model_type", "TidalCurrentCycle"))
    call cfcheck( nf90_put_att(ofid,nf90_global,"history","fvcom2gnome generated this file")) 

    !define the dimensions 
    call cfcheck(nf90_def_dim(ofid,"nbi",4, nbi_did) )
    call cfcheck(nf90_def_dim(ofid,"nbnd",N_bsegs, nbnd_did) )
    call cfcheck(nf90_def_dim(ofid,"time",nf90_unlimited, time_did) )
    call cfcheck(nf90_def_dim(ofid,"node",N_verts, node_did) )
    call cfcheck(nf90_def_dim(ofid,"nele",1, nele_did) )

    !variables
    call cfcheck( nf90_def_var(ofid,"bnd",nf90_int,(/nbi_did,nbnd_did/),bnd_vid) )
    call cfcheck( nf90_put_att(ofid, bnd_vid,"long_name","Boundary Segment Node List") )
    call cfcheck( nf90_put_att(ofid, bnd_vid,"units","nondimensional") )

    call cfcheck( nf90_def_var(ofid,"time",nf90_float,(/time_did/), time_vid) )
    call cfcheck( nf90_put_att(ofid, time_vid,"long_name","Time") )
    call cfcheck( nf90_put_att(ofid, time_vid,"units",tstring) )
    if(with_basedate)then
      call cfcheck( nf90_put_att(ofid, time_vid,"base_date",trim(basedate)) )
    endif

    call cfcheck( nf90_def_var(ofid,"lon",nf90_float,(/node_did/), lon_vid) )
    call cfcheck( nf90_put_att(ofid, lon_vid,"long_name","Longitude") )
    call cfcheck( nf90_put_att(ofid, lon_vid,"units","degrees_east") )

    call cfcheck( nf90_def_var(ofid,"lat",nf90_float,(/node_did/), lat_vid) )
    call cfcheck( nf90_put_att(ofid, lat_vid,"long_name","Latitude") )
    call cfcheck( nf90_put_att(ofid, lat_vid,"units","degrees_north") )

    allocate(un(N_verts)) ; un = zero
    call cfcheck( nf90_def_var(ofid,"u",nf90_float,(/node_did,time_did/), u_vid) )
    call cfcheck( nf90_put_att(ofid, u_vid,"long_name","Eastward Water Velocity") )
    call cfcheck( nf90_put_att(ofid, u_vid,"units","m/s") )
    call cfcheck( nf90_put_att(ofid, u_vid,"missing_value",-99999.) )
    call cfcheck( nf90_put_att(ofid, u_vid,"dry_value",999.) )

    allocate(vn(N_verts)) ; vn = zero
    call cfcheck( nf90_def_var(ofid,"v",nf90_float,(/node_did,time_did/), v_vid) )
    call cfcheck( nf90_put_att(ofid, v_vid,"long_name","Northward Water Velocity") )
    call cfcheck( nf90_put_att(ofid, v_vid,"units","m/s") )
    call cfcheck( nf90_put_att(ofid, v_vid,"missing_value",-99999.) )
    call cfcheck( nf90_put_att(ofid, v_vid,"dry_value",999.) )

    allocate(z(N_verts)) ; z = zero
    call cfcheck( nf90_def_var(ofid,"z",nf90_float,(/node_did,time_did/), z_vid) )
    call cfcheck( nf90_put_att(ofid, z_vid,"long_name","Tidal Elevations")) 
    call cfcheck( nf90_put_att(ofid, z_vid,"units","m") )
    call cfcheck( nf90_put_att(ofid, z_vid,"missing_value",-99999.) )
    call cfcheck( nf90_enddef(ofid) )

    
    !dump statics
    call cfcheck(nf90_put_var(ofid, bnd_vid,GNOME_bndry))
    call cfcheck(nf90_put_var(ofid, lon_vid,lon))
    call cfcheck(nf90_put_var(ofid, lat_vid,lat))

    !close up
    call cfcheck( nf90_close(ofid) )
    

  End Subroutine Setup_Output
    

  !=====================================================================
  ! open the output file type and setup dimensions/variables 
  !=====================================================================
  Subroutine Setup_Output_NCTOOLS
    type(ncfile), pointer :: ncf
    type(ncatt) , pointer :: att
    type(ncvar) , pointer :: var
    type(ncdim) , pointer :: dim

    !create new file and set name
    ncf => new_file()
    ncf%fname = trim(outputfile)
    ncf%ftime=>new_ftime()
    ncf%ftime%next_stkcnt = 0
         
    !make dimensions
    !DIM_nele_out => NC_MAKE_DIM(name='nele', len=1)
    DIM_node_out => NC_MAKE_DIM(name='node', len=N_verts)
    DIM_nbnd_out => NC_MAKE_DIM(name='nbnd', len=N_bsegs)
    DIM_nbi_out  => NC_MAKE_DIM(name='nbi',  len=4)
    DIM_time_out => NC_MAKE_DIM(name='time', len=NF90_UNLIMITED)

    !add global attributes    
    att => nc_make_att(name='file_type',values="FEM") 
    ncf => add(ncf,att)
    att => nc_make_att(name='Conventions',values="COARDS") 
    ncf => add(ncf,att)
    att => nc_make_att(name='grid_type',values="Triangular") 
    ncf => add(ncf,att)
    !att => nc_make_att(name='model_type',values="TidalCurrentCycle") 
    !ncf => add(ncf,att)
    att => nc_make_att(name='history',values="fvcom2gnome generated this file") 
    ncf => add(ncf,att)
         
    !make dimensions
    nc_outfile => ncf

    ! boundary array 
    var => nc_make_avar(name='bnd',VALUES=GNOME_bndry,dim1=DIM_nbi_out,dim2=DIM_nbnd_out)
    att => nc_make_att(name='long_name',values='Boundary Segment Node List')
    var => add(var,att)
    att => nc_make_att(name='units',values='nondimensional')
    var => add(var,att)
    ncf => add(ncf,var)


    ! time 
    var => nc_make_avar(name='time',VALUES=gtime,dim1=DIM_time_out)
    att => nc_make_att(name='long_name',values='Time')
    var => add(var,att)  
    att => nc_make_att(name='units',values=trim(tstring))
    var => add(var,att)  
    ncf => add(ncf,var)

    ! Longitude
    var => nc_make_avar(name='lon',VALUES=lon,dim1=DIM_node_out)
    att => nc_make_att(name='long_name',values='Longitude')
    var => add(var,att)  
    att => nc_make_att(name='units',values='degrees_east')
    var => add(var,att)  
    ncf => add(ncf,var)

    ! Latitude 
    var => nc_make_avar(name='lat',VALUES=lat,dim1=DIM_node_out)
    att => nc_make_att(name='long_name',values='Latitude')
    var => add(var,att)  
    att => nc_make_att(name='units',values='degrees_north')
    var => add(var,att)  
    ncf => add(ncf,var)

    ! U values on the vertices
    allocate(un(N_verts)) ; un = zero
    var  => nc_make_avar(name='u',&
         & values=un, dim1= dim_node_out ,dim2 = dim_TIME_out)
    att  => nc_make_att(name='long_name',values='Eastward Water Velocity')
    var  => add(var,att)
    att  => nc_make_att(name='units',values='m/s')
    var  => add(var,att)
    att  => nc_make_att(name='missing_value',values=-99999.)
    var  => add(var,att)
    att  => nc_make_att(name='dry_value',values=999.)
    var  => add(var,att)
    ncf  => add(ncf,var)

    ! V values on the vertices
    allocate(vn(N_verts)) ; vn = zero
    var  => nc_make_avar(name='v',&
         & values=vn, dim1= dim_node_out, dim2 = dim_TIME_out)
    att  => nc_make_att(name='long_name',values='Northward Water Velocity')
    var  => add(var,att)
    att  => nc_make_att(name='units',values='m/s')
    var  => add(var,att)
    att  => nc_make_att(name='missing_value',values=-99999.)
    var  => add(var,att)
    att  => nc_make_att(name='dry_value',values=999.)
    var  => add(var,att)
    ncf  => add(ncf,var)

    ! zeta values on the vertices
    allocate(z(N_verts)) ; z = zero
    var  => nc_make_avar(name='z',&
         & values=z, dim1= dim_node_out, dim2 = dim_TIME_out)
    att  => nc_make_att(name='long_name',values='Tidal Elevations')
    var  => add(var,att)
    att  => nc_make_att(name='units',values='m')
    var  => add(var,att)
    att  => nc_make_att(name='missing_value',values=-99999.)
    var  => add(var,att)
    ncf  => add(ncf,var)

    !a junk array to force number of elements to dump 
    !allocate(junk(1)) ; junk = zero
    !var  => nc_make_avar(name='junk',&
    !     & values=junk, dim1= dim_nele_out ,dim2 = dim_TIME_out)
    !ncf  => add(ncf,var)


    !write the file
    call nc_write_file(ncf) 

    !connect the file pointer
    nc_outfile => ncf
    

  End Subroutine Setup_Output_NCTOOLS
    
  !=====================================================================
  ! Loop over frames, read velocity, interp to nodes and dump 
  !=====================================================================
  Subroutine Extract 
    Integer :: i,ii
    Type(ncvar), pointer :: var
    Logical :: found

    write(ipt,*) '============ extraction ======================='
    ii = 1
    do i=fbeg,fend,fint
      write(ipt,*)'extracting frame: ',i,ii

      !-------------------------------------------------------
      !read FVCOM data from frame i
      !-------------------------------------------------------
      !read time 
      var => find_var(nc_infile,'time',FOUND)
      if(.not.FOUND) call fatal_error("no variable 'time' in inputfile")
      call nc_connect_avar(var, gtime)
      call nc_read_var(var,i)
      nullify(var)

      !displace time to start from 00:00:00
      gtime = gtime - T0

      !surface elevation  
      var => find_var(nc_infile,'zeta',FOUND)
      if(.not.FOUND) call fatal_error("no variable 'zeta' in inputfile")
      call nc_connect_avar(var, z)
      call nc_read_var(var,i)
      nullify(var)

      ! ------- 3D option --------------
      if(fvcom_dim == 3)then
        !read u and interpolate to nodes
        var => find_var(nc_infile,'u',FOUND)
        if(.not.FOUND) call fatal_error("no variable 'u' in inputfile")
        call nc_connect_avar(var, u)
        call nc_read_var(var,i)
        call elem_2_vtx(u(:,1),un)
        nullify(var)

        !read v and interpolate to nodes
        var => find_var(nc_infile,'v',FOUND)
        if(.not.FOUND) call fatal_error("no variable 'v' in inputfile")
        call nc_connect_avar(var, v)
        call nc_read_var(var,i)
        call elem_2_vtx(v(:,1),vn)
        nullify(var)
      else   !2D
        !read ua and interpolate to nodes
        var => find_var(nc_infile,'ua',FOUND)
        if(.not.FOUND) call fatal_error("no variable 'ua' in inputfile")
        call nc_connect_avar(var, ua)
        call nc_read_var(var,i)
        call elem_2_vtx(ua(:),un)
        nullify(var)

        !read va and interpolate to nodes
        var => find_var(nc_infile,'va',FOUND)
        if(.not.FOUND) call fatal_error("no variable 'va' in inputfile")
        call nc_connect_avar(var, va)
        call nc_read_var(var,i)
        call elem_2_vtx(va(:),vn)
        nullify(var)
      endif
      

      !-------------------------------------------------------
      !interpolate surface velocity to nodes (u,v) => (un,vn)
      !-------------------------------------------------------

      !-------------------------------------------------------
      !dump to GNOME file  at frame ii
      !-------------------------------------------------------
!     NC tools output
      !nc_outfile%ftime%next_stkcnt = nc_outfile%ftime%next_stkcnt + 1
      !call nc_write_file(nc_outfile)
!     NC tools output

!     OLD SCHOOL NF90 OUTPUT
       call cfcheck( nf90_open(trim(outputfile),nf90_write,ofid) )
       call cfcheck( nf90_put_var(ofid, time_vid,gtime ,START=(/ii/)))
       call cfcheck( nf90_put_var(ofid, u_vid,un ,START=(/1,ii/)))
       call cfcheck( nf90_put_var(ofid, v_vid,vn ,START=(/1,ii/)))
       call cfcheck( nf90_put_var(ofid, z_vid,z ,START=(/1,ii/)))
       call cfcheck (nf90_close(ofid))
!     OLD SCHOOL NF90 OUTPUT

       ii = ii + 1


    end do

  End Subroutine Extract
  
  !=====================================================================
  ! Cleanup 
  !=====================================================================
  Subroutine Cleanup 
    write(ipt,*) '============ done and cleaning up ============='
    !add nullify/deallocates here
  End Subroutine Cleanup

  !=====================================================================
  ! Calculate the area of the cell, make sure it is positive 
  !=====================================================================
  Function cell_area(xin,yin) 
    real(sp), intent(in) :: xin(3)
    real(sp), intent(in) :: yin(3)
    real(sp) :: cell_area

    !note negative sign in cross product since FVCOM connectivity is clockwise
    cell_area = -ahalf*( (xin(2)-xin(1))*(yin(3)-yin(1)) - &
                (xin(3)-xin(1))*(yin(2)-yin(1)) )
    if(cell_area < zero)then
      call fatal_error("area is negative, mesh connectivity may be backwards")
    endif

  End Function cell_area 

  !=====================================================================
  ! Area-Weighted Interpolation of velocity from elements to nodes
  !=====================================================================
  Subroutine Elem_2_Vtx(fc,fv) 
    Real(sp), intent(in ) :: fc(1:N_elems)
    Real(sp), intent(out) :: fv(1:N_verts)
    Integer :: i,j,cell

    Do i=1,N_verts
      fv(i) = sum( art(nbve(i,1:ntve(i)))*fc(nbve(i,1:ntve(i)))  )/art1(i) 
    End Do

  End Subroutine Elem_2_Vtx 


  !-----------------------------------------------
  ! runtime errors  - netcdf
  !-----------------------------------------------
  subroutine cfcheck(status)
      use netcdf
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop
      end if
  end subroutine cfcheck


End Module mod_fvcom2gnome
