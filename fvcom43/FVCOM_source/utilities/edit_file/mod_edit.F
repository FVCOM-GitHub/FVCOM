module mod_edit
  use lims
  use mod_utils
  use mod_par
  use control
  use mod_prec
  use mod_ncll
  use mod_nctools
  use all_vars
  implicit none


  CHARACTER(LEN=120) :: filename


  TYPE(NCFILE), POINTER :: NCF
  TYPE(NCVAR), POINTER :: VAR

  TYPE(NCATT), POINTER :: ATT
  TYPE(NCDIM), POINTER :: DIM



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
    write(IPT,*) "! "
    write(IPT,*) "! --filename=XXXX : user specified file name"
    write(IPT,*) "! "
    write(IPT,*) "! Beware - this is a case specific program"
    write(IPT,*) "! It should be edited for each use..."
    
    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT


  SUBROUTINE EDIT
    IMPLICIT NONE
    LOGICAL FOUND
    integer status, I

    INTEGER, ALLOCATABLE :: IBCN(:)

    NCF => NEW_FILE()
    NCF%FNAME=trim(filename)
    NCF%WRITABLE=.true.
    ! OPEN THE FILE AND LOAD METADATA       

    Call NC_OPEN(NCF)
    CALL NC_LOAD(NCF)

    if (dbg_set(dbg_io)) then
       CALL print_file(NCF)
    END if

!!$ THIS EXAMPLE WOULD SET THE OPEN BOUNDARY NODE NUMBER IN A RESTART FILE
!!$    DIM => FIND_DIM(NCF,'nobc',FOUND)
!!$    IF (.NOT. FOUND) CALL FATAL_ERROR("FINLE DOES NOT HAVE DIMENSION 'nobc' ?")
!!$
!!$    VAR => find_var(NCF,'obc_nodes',FOUND)
!!$    IF (.NOT. FOUND) CALL FATAL_ERROR("FINLE DOES NOT HAVE VARIABLE 'obc_nodes' ?")
!!$
!!$    ALLOCATE(IBCN(DIM%DIM))
!!$    DO I = 1, DIM%DIM
!!$       IBCN(I) = I
!!$    END DO
!!$
!!$    CALL NC_CONNECT_AVAR(VAR,IBCN)
!!$
!!$    call NC_WRITE_VAR(VAR,.TRUE.,.FALSE.,0)


    CALL NC_CLOSE(NCF)


  END SUBROUTINE EDIT


end module mod_edit
