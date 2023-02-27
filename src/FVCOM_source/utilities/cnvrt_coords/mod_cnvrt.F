module mod_cnvrt
  use lims
  use mod_utils
  use mod_par
  use control
  use mod_prec
  use mod_ncll
  use mod_nctools
  implicit none

  INTEGER, PARAMETER :: makesphere=-1
  INTEGER, PARAMETER :: makecartesian=1


  CHARACTER(len=80) :: cref
  CHARACTER(LEN=120) :: filename
  integer :: ctype
  integer :: nNode



  TYPE(NCFILE), POINTER ::NCF
  TYPE(NCVAR), POINTER :: VX
  TYPE(NCVAR), POINTER :: VY
  TYPE(NCVAR), POINTER :: Vlat
  TYPE(NCVAR), POINTER :: Vlon

  TYPE(NCATT), POINTER :: ATT
  TYPE(NCDIM), POINTER :: DIM



  REAL(SP), ALLOCATABLE :: X(:),Y(:),LAT(:),LON(:)

contains



  SUBROUTINE GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)
    use mod_sng


    character(len=*), INTENT(IN)::CVS_Id  ! [sng] CVS Identification
    character(len=*), INTENT(IN)::CVS_Date ! [sng] Date string
    character(len=*), INTENT(IN)::CVS_Name ! [sng] File name string
    character(len=*), INTENT(IN)::CVS_Revision ! [sng] File revision string

    character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
    ! Command-line parsing
    character(80)::arg_val ! [sng] command-line argument value
    character(200)::cmd_ln ! [sng] command-line
    character(80)::opt_sng ! [sng] Option string
    character(2)::dsh_key ! [sng] command-line dash and switch
    character(200)::prg_ID ! [sng] Program ID

    integer::arg_idx ! [idx] Counting index
    integer::arg_nbr ! [nbr] Number of command-line arguments
    integer::opt_lng ! [nbr] Length of option

    ! Main code
    call ftn_strini(cmd_ln) ! [sng] sng(1:len)=NUL

    call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
    call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID

    arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments

    if (arg_nbr .LE. 0 ) then
       if(MSR) WRITE(IPT,*) "You must specify an arugument:"
       if(MSR) Call MYHelpTxt
       call PSHUTDOWN
    end if

    arg_idx=1 ! [idx] Counting index
    do while (arg_idx <= arg_nbr)
       call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg, increment arg_idx
       dsh_key=arg_val(1:2) ! [sng] First two characters of option
       if (dsh_key == "--") then
          opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
          if (opt_lng <= 0) then
             if(MSR) write(IPT,*) "Long option has no name"
             call PSHUTDOWN
          end if

          opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
          if (dbg_lvl >= dbg_io) then
             if(MSR) write (6,"(5a,i3)") prg_nm(1:ftn_strlen(prg_nm)), &
                  ": DEBUG Double hyphen indicates multi-character option: ", &
                  "opt_sng = ",opt_sng(1:ftn_strlen(opt_sng)),", opt_lng = ",opt_lng
          end if
          if (opt_sng == "dbg" .or. opt_sng == "dbg_lvl" ) then
             call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level

             !          else if (opt_sng == "dbg_par" .or.opt_sng == "Dbg_Par"&
             !               & .or.opt_sng == "DBG_PAR") then

             !             dbg_par = .true.

          else if (opt_sng == "FILENAME" .or.opt_sng == "filename"&
               & .or.opt_sng == "Filename") then

             call ftn_arg_get(arg_idx,arg_val,filename) ! [sng] Input file
             filename=filename(1:ftn_strlen(filename))
             ! Convert back to a fortran string!

          else if (opt_sng == "REFERENCE" .or.opt_sng == "reference"&
               & .or.opt_sng == "Reference") then

             call ftn_arg_get(arg_idx,arg_val,cref) ! [sng] Input file
             cref=cref(1:ftn_strlen(cref))
             ! Convert back to a fortran string!

          else if (opt_sng == "SPHR2CART" .or.opt_sng == "sphr2cart"&
               & .or.opt_sng == "Sphr2cart") then

             ctype = makecartesian

          else if (opt_sng == "CART2SPHR" .or.opt_sng == "cart2sphr"&
               & .or.opt_sng == "Cart2sphr") then

             ctype =makesphere

          else if (opt_sng == "help" .or.opt_sng == "HELP" .or. opt_sng&
               & == "Help") then

             if(MSR) call MYHelpTxt
             call PSHUTDOWN

          else ! Option not recognized
             arg_idx=arg_idx-1 ! [idx] Counting index
             if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
          endif ! endif option is recognized
          ! Jump to top of while loop
          cycle ! C, F77, and F90 use "continue", "goto", and "cycle"
       endif ! endif long option

       if (dsh_key == "-V" .or.dsh_key == "-v" ) then

          if(MSR) write(IPT,*) prg_id
          call PSHUTDOWN

       else if (dsh_key == "-H" .or.dsh_key == "-h" ) then

          if(MSR) Call MYHelpTxt
          Call PSHUTDOWN

       else ! Option not recognized
          arg_idx=arg_idx-1 ! [idx] Counting index
          if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
       endif ! endif arg_val


    end do ! end while (arg_idx <= arg_nbr)

    CALL dbg_init(IPT_BASE,.false.)

  END SUBROUTINE GET_COMMANDLINE

  SUBROUTINE MYHELPTXT
    IMPLICIT NONE


    write(IPT,*) "! ARGUMENTS FOR cnvrt_coords:"
    write(IPT,*) "! --sphr2cart"
    write(IPT,*) "!    or"
    write(IPT,*) "! --cart2sphr"
    write(IPT,*) "! "
    write(IPT,*) "! --reference=XXXXX : user specified projection type"
    write(IPT,*) "! "
    write(IPT,*) "! --filename=XXXX : user specified file name"

    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT


  SUBROUTINE CONVERT
    IMPLICIT NONE
    LOGICAL FOUND
    integer status

    NCF => NEW_FILE()
    NCF%FNAME=trim(filename)
    NCF%WRITABLE=.true.
    ! OPEN THE FILE AND LOAD METADATA       

    Call NC_OPEN(NCF)
    CALL NC_LOAD(NCF)

    if (dbg_set(dbg_io)) then
       CALL print_file(NCF)
    END if


    DIM => FIND_DIM(NCF,'node',FOUND)
    IF (.NOT. FOUND) CALL FATAL_ERROR("FINLE DOES NOT HAVE DIMENSION 'node' ?")

    nNode= Dim%dim

    allocate(X(nNode))
    allocate(Y(nNode))
    allocate(LAT(nNode))
    allocate(LON(nNode))


    VX => find_var(NCF,'x',FOUND)
    IF (.NOT. FOUND) CALL FATAL_ERROR("FINLE DOES NOT HAVE VARIABLE 'x' ?")
    CALL NC_CONNECT_AVAR(VX,X)
        
    VY => find_var(NCF,'y',FOUND)
    IF (.NOT. FOUND) CALL FATAL_ERROR("FINLE DOES NOT HAVE VARIABLE 'y' ?")
    CALL NC_CONNECT_AVAR(VY,Y)
        
    VLAT => find_var(NCF,'lat',FOUND)
    IF (.NOT. FOUND) CALL FATAL_ERROR("FINLE DOES NOT HAVE VARIABLE 'lat' ?")
    CALL NC_CONNECT_AVAR(VLAT,lat)
        
    VLON => find_var(NCF,'lon',FOUND)
    IF (.NOT. FOUND) CALL FATAL_ERROR("FINLE DOES NOT HAVE VARIABLE 'lon' ?")
    CALL NC_CONNECT_AVAR(VLON,lon)
    

    select case(ctype)
    case(makesphere)
       
       ATT => FIND_ATT(NCF,'CoordinateSystem',FOUND)
       IF (FOUND) THEN
          If(ATT%CHR(1) == "GeoReferenced") THEN
             CALL WARNING("This file says it is already 'GeoReferenced' ?",&
                  & "Converting from Cartesian to LAT/LON and overwrite:y")
          ELSE
             ATT%CHR(1) = "GeoReferenced"
          END If

       ELSE
          ATT => NC_MAKE_ATT(name='CoordinateSystem',values="GeoReferenced" )
       END IF

       CALL NC_READ_VAR(VX)
       CALL NC_READ_VAR(VY)

       CALL METERS2DEGREES(X,Y,cref,LON,LAT,nNode)
       
       call NC_WRITE_VAR(VLAT,DUMP=.TRUE.,COLLECTDATA=.FALSE.,COLLECTOR=MYID)
       call NC_WRITE_VAR(VLON,DUMP=.TRUE.,COLLECTDATA=.FALSE.,COLLECTOR=MYID)


       status = nf90_redef(NCF%ncid)
       CALL write_att_type(NCF%ncid,NF90_GLOBAL,ATT)
       status = nf90_enddef(NCF%ncid)


    case(makecartesian)
       
       ATT => FIND_ATT(NCF,'CoordinateSystem',FOUND)
       IF (FOUND) THEN
          If(ATT%CHR(1) == "GeoReferenced") THEN
             ATT%CHR(1) = "Cartesian"
          ELSE
             CALL WARNING("This file says it is already 'Cartesian' ?",&
                  & "Converting from LAT/LON to Cartesian and overwrite:")

          END If

       ELSE
          ATT => NC_MAKE_ATT(name='CoordinateSystem',values="GeoReferenced" )
       END IF

       CALL NC_READ_VAR(VLAT)
       CALL NC_READ_VAR(VLON)

       CALL DEGREES2METERS(LON,LAT,cref,X,Y,nNode)

       call NC_WRITE_VAR(VX,DUMP=.TRUE.,COLLECTDATA=.FALSE.,COLLECTOR=MYID)
       call NC_WRITE_VAR(VY,DUMP=.TRUE.,COLLECTDATA=.FALSE.,COLLECTOR=MYID)


       status = nf90_redef(NCF%ncid)
       CALL write_att_type(NCF%ncid,NF90_GLOBAL,ATT)
       status = nf90_enddef(NCF%ncid)

    case default
       call fatal_error("unknown conversion type ???")
    END select

    CALL NC_CLOSE(NCF)


  END SUBROUTINE CONVERT

end module mod_cnvrt
