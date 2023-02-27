!/===========================================================================/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/
! Purpose: Library of Fortran string manipulation routines
!
! Copyright (C) 1997--2005 Charlie Zender
! License Summary: X11-style free, non-restrictive, non-copyleft. 
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, and/or sell copies of the Software, and to permit persons
! to whom the Software is furnished to do so, provided that the above
! copyright notice(s) and this permission notice appear in all copies of
! the Software and that both the above copyright notice(s) and this
! permission notice appear in supporting documentation.
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT
! OF THIRD PARTY RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
! HOLDERS INCLUDED IN THIS NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL
! INDIRECT OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING
! FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
! NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
! WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
! Except as contained in this notice, the name of a copyright holder
! shall not be used in advertising or otherwise to promote the sale, use
! or other dealings in this Software without prior written authorization
! of the copyright holder.

! The original author of this software, Charlie Zender, wants to receive
! your suggestions, thanks, bug-reports, and patches to improve it.
! Charlie Zender <zender at uci dot edu>
! Department of Earth System Science
! University of California, Irvine
! Irvine, CA 92697-3100

! These routines contain the functionality required to parse command-line arguments
! http://www.winteracter.com/f2kcli/index.htm provides f2kcli module 
! required by some (non-UNIX) F9X compilers. F2K provides functions natively.
! A standalone test program, sng.F90, is used to demonstrate sng_mdl.F90

! Recursive I/O:
! Do not embed print statments in a function whose result is being printed

! NyL97 p. 318: POT, TIP, APT
! For formatted I/O, characters are truncated or blanks are added,
! depending on whether the field width is too small or too large.
! For input, truncation occurs on the left, and blank padding on the right;
! for output, truncation occurs on the right, and blak padding on the left

! Notes on PRC_DBL, auto-promotion, overloading:
! Assume PRC_DBL is defined when compiler is passed autopromotion flags ("-r8")
! This means all floats are promoted to double precision
! As a result, ftn_arg_get_flt() is identical to ftn_arg_get_dbl() in precision
! Thus compiler cannot disambiguate overloading these two functions
! Result is OK on all known compilers except Lahey lf95 which complains
! Solution is to define ftn_arg_get_flt() only when auto-promotion is not enabled
! Alternate solution is never use auto-promotion with lf95

! Usage:
!use sng_mdl ! [mdl] String manipulation

! Set F2K token on compilers known to support Fortran 2000/2003 intrinsics
! gfortran 4.1.2 does not recognize __GFORTRAN__ so use __GNUC__ instead
#if defined(__G95__) || defined(__GNUC__) || defined(F2KCLI)
# define F2K
#endif /* !g95 || gfortran */

module mod_sng ! [mdl] String manipulation
#ifdef F2KCLI
  ! If Fortran2000 compiler does not have command_argument_count() and 
  ! get_command_argument() intrinsics, then use Winteracter's f2kcli module
  use f2kcli ! [mdl] F2K Command-Line Interpreter (command_argument_count()...)
#endif /* !F2KCLI */
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private
  private::ftn_arg_get_dbl,ftn_arg_get_int,ftn_arg_get_lgc,ftn_arg_get_sng
  private::ftn_arg_get_flt
  
  ! Overloaded command-line argument retrieval functions
  interface ftn_arg_get
     module procedure ftn_arg_get_dbl,ftn_arg_get_int,ftn_arg_get_lgc,ftn_arg_get_sng
     
     module procedure ftn_arg_get_flt
  end interface ! ftn_arg_get
  
contains
  
#ifndef F2K
  ! Declare forward-compatible F2K functions
  integer function command_argument_count() ! [fnc] Fortran intrinsic (F2K)
    ! Purpose: Return number of command-line arguments (not counting program name)
    ! Standardized F2K replacement for F77/F9X iargc()
    interface
       integer function iargc() ! [fnc] Fortran intrinsic
       end function iargc
    end interface
    command_argument_count=iargc()
    return
  end function command_argument_count
  
  subroutine get_command_argument(arg_idx,arg_val)
    ! Purpose: Return specified command-line argument
    ! Standardized F2K replacement for F77/F9X iargc()
    integer,intent(in)::arg_idx ! I [idx] Argument counter 
    character(len=*),intent(out)::arg_val ! I/O [sng] String to copy into opt_val
    interface
       subroutine getarg(idx,val) ! [fnc] Fortran intrinsic
         integer,intent(in)::idx
         character(len=*),intent(out)::val
       end subroutine getarg
    end interface
    call getarg(arg_idx,arg_val)
    return
  end subroutine get_command_argument
#endif /* F2K */

  integer function ftn_strlen(sng) ! [nbr] Length of string
    ! Purpose: Return length of string preceding first null character
    ! Fortran intrinsic len(sng) is returned if string is not NUL-terminated
    ! Prototype:
    ! integer ftn_strlen ! [nbr] Length of string
    implicit none
    ! Input
    character(len=*),intent(in)::sng
    ! Local
    ! Main Code
    ! Normally, length of string is position of first insignificant character minus one
    ftn_strlen=ftn_strfic(sng)-1 ! [idx] First NUL or 8-bit character position
    ! String length may not exceed buffer size
    if (ftn_strlen > len(sng)) ftn_strlen=len(sng)
    ! String length may not be less than 0
    ! Strings initialized to '' should have length 0 not -1
    if (ftn_strlen < 0 .or. sng == '') ftn_strlen=0
    return
  end function ftn_strlen ! end ftn_strlen()
  
  integer function ftn_opt_lng_get(sng)
    ! Purpose: Return length of option string preceding space or equals character
    ! Any preceding dashes are NOT counted towards option length
    ! Fortran intrinsic len(sng) is returned if string is not NUL-terminated
    ! Prototype:
    ! integer ftn_opt_lng_get ! [nbr] Length of option
    ! Usage:
    ! Call ftn_opt_lng_get with the full option string, e.g.,
    ! ftn_opt_lng_get('--dbg_lvl=5') and ftn_opt_lng_get will return the
    ! length of the option string 
    ! opt_sng=arg_val(3:2+ftn_opt_lng_get(arg_val)) ! [sng] Option string
    ! if (opt_sng == 'dbg') then ...
    implicit none
    ! Input
    character(len=*),intent(in)::sng
    ! Local
    integer idx
    integer len_sng
    integer::dsh_nbr=0        ! [nbr] Number of dashes
    ! Main Code
    ! Length of option string is one less than position of first space or equal sign
    len_sng=len(sng)
    if (len_sng >= 1) then
       if (sng(1:1) == '-') dsh_nbr=1 ! [nbr] Number of dashes
    endif
    if (len_sng >= 2) then
       if (sng(1:2) == '--') dsh_nbr=2 ! [nbr] Number of dashes
    endif
    do idx=1,len_sng
       ! Look for first abnormal character
       if (sng(idx:idx) == ' ' .or. sng(idx:idx) == '=') goto 100
    end do                    ! end loop over characters
100 continue
    ftn_opt_lng_get=idx-1-dsh_nbr
    ! Strings initialized to '' should have length 0
    if (ftn_opt_lng_get < 0 .or. sng == '') ftn_opt_lng_get=0
    return
  end function ftn_opt_lng_get                       ! end ftn_opt_lng_get()
  
  subroutine ftn_getarg_wrp( & ! [sbr] Call getarg() and increment arg_idx
       arg_idx, & ! I/O [idx] Argument counter
       arg_val) ! I/O [sng] String to copy into opt_val
    ! Purpose: Wrapper for getarg() intrinsic that increments argument counter
    ! Argument index is incremented after getarg() is called
    ! Usage:
    ! call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg() and increment arg_idx
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Input
    ! Input/Output
    integer,intent(inout)::arg_idx ! I/O [idx] Argument counter 
    character(len=*),intent(inout)::arg_val ! I/O [sng] String to copy into opt_val
    ! Output
    ! Local
    ! Main Code
    call get_command_argument(arg_idx,arg_val)
    if (dbg_lvl >= dbg_io) write (6,'(2a,i2,2a)') prg_nm(1:ftn_strlsc(prg_nm)),   &
         ': DEBUG ftn_getarg_wrp() reports arg_idx = ',arg_idx,', arg_val = ',arg_val(1:len_trim(arg_val))
    ! Increment counter for next read
    arg_idx=arg_idx+1
    return
  end subroutine ftn_getarg_wrp
  
  subroutine ftn_getarg_err( & ! [sbr] Error handler for getarg()
       arg_idx, & ! I/O [idx] Argument counter
       arg_val) ! I/O [sng] String to copy into opt_val
    ! Purpose: Handle option-not-found errors for getopt()
    ! Routine assumes arg_idx points to arg_val
    ! Typically arg_idx is advanced before option recognition is attempted,
    ! thus arg_idx may need to be decremented before ftn_getarg_err() is called
    ! Usage:
    ! call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Input
    ! Input/Output
    integer,intent(inout)::arg_idx ! I/O [idx] Argument counter 
    character(len=*),intent(inout)::arg_val ! I/O [sng] String to copy into opt_val
    ! Output
    ! Local
    integer arg_nbr           ! [nbr] Number of command-line arguments
    ! integer exit_status       ! [enm] Program exit status (non-standard Fortran)
    ! Main Code
    arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments
    write (6,'(2a)') prg_nm(1:ftn_strlsc(prg_nm)),': ERROR Option not recognized'
    write (6,'(2a)') prg_nm(1:ftn_strlsc(prg_nm)),': HINT Option syntax is nearly POSIX/GNU compliant:'
    write (6,'(2a)') prg_nm(1:ftn_strlsc(prg_nm)), &
         ': Short option syntax is dash-key-space-value, e.g., -D 2'
    write (6,'(3a)') prg_nm(1:ftn_strlsc(prg_nm)), &
         ': Long option syntax is dash-dash-key-space-value, e.g., --dbg 2', &
         ', or dash-dash-key-equal-value, e.g., --dbg=2 (preferred)'
    ! Subtract 1 from arg_idx because rule is that arg_idx is incremented before option recognition is attempted
    write (6,'(2a,i2,2a)') prg_nm(1:ftn_strlsc(prg_nm)),   &
         ': DEBUG ftn_getarg_err() reports arg_idx = ',arg_idx,', arg_val = ',arg_val(1:len_trim(arg_val))
    if (arg_idx+1 <= arg_nbr) then
       arg_idx=arg_idx+1      ! [idx] Counting index
       call get_command_argument(arg_idx,arg_val)
       write (6,'(2a,i2,2a)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG Next argument value is arg_val(',arg_idx,') = ',arg_val
    else
       write (6,'(2a)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG No arguments follow this one'
    endif                     ! endif not last argument
    ! exit_status=-1            ! [enm] Program exit status (non-standard Fortran)
    ! call exit(exit_status)    ! [enm] Exit with current exit status (non-standard Fortran)
    stop 'Exit on error from ftn_getarg_err()'
    return
  end subroutine ftn_getarg_err
  
  ! Notes common to all arg_get() routines:
  ! arg_val is assumed to be present as a command-line line option 
  ! On input, first two characters of arg_val are assumed to be dashes,
  ! although this would not be too difficult to change.
  ! Routine accepts options specified in POSIX format: "--opt_sng val" or "--opt_sng=val" (preferred)
  ! Case 1: Option string and value are separated by space " " conjunction
  ! On input, arg_val is assumed to hold option string with dashes, e.g., "--opt_sng"
  ! On output, arg_val holds corresponding argument, e.g., "73"
  ! Routine uses getarg() to obtain arg_val, then increments arg_idx
  ! Case 2: Option string and value are separated by equals "=" conjunction
  ! On input and output, arg_val holds entire option and value, e.g., "--opt_sng=73"
  ! Routines do not call getarg() to change arg_val and do not increment arg_idx
  ! Routines set optional logical argument opt_flg to true if called
  
  subroutine ftn_arg_get_dbl( & ! [sbr] Process double-valued command-line argument
       arg_idx, & ! I/O [idx] Argument counter
       arg_val, & ! I/O [sng] Double to copy into opt_val
       opt_val, & ! O [frc] Variable to receive copy of arg_val
       opt_flg) ! O [frc] Variable set by command-line
    ! Purpose: Copy command-line double argument arg_val into variable opt_val
    ! Usage:
    ! call ftn_arg_get_dbl(arg_idx,arg_val,opt_val) ! [sbr] Process double-valued command-line argument
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='ftn_arg_get_dbl' ! [sng] Subroutine name
    ! Commons
    ! Input
    ! Input/Output
    integer,intent(inout)::arg_idx ! I/O [idx] Argument counter 
    character(len=*),intent(inout)::arg_val ! I/O [sng] Double to copy into opt_val
    ! Output
    real(selected_real_kind(p=12)),intent(out)::opt_val ! O [frc] Variable to receive copy of arg_val
    logical,optional,intent(out)::opt_flg ! O [flg] Variable set by command-line
    ! Local
    integer arg_lng           ! [nbr] Length of argument
    integer opt_lng           ! [nbr] Length of option
    integer arg_val_srt_idx   ! [idx] Starting position of argument value
    logical opt_cnj_is_spc ! [flg] Option conjunction is space character
    ! Main Code
    ! Print any diagnostics before current value of arg_val is overwritten
    arg_lng=len_trim(arg_val) ! [nbr] Length of argument
    opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
    if (dbg_lvl >= dbg_io) then
       write (6,'(2a,i2,3a,i2)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() reports arg_idx = ',arg_idx,  &
            ', Full option = ''',arg_val(1:arg_lng),''', Length = ',arg_lng
    endif                     ! endif dbg
    ! Determine whether format is --opt_sng=val or --opt_sng val
    ! Add two to cover preliminary dashes
    if (opt_lng+2 < arg_lng) then
       opt_cnj_is_spc=.false. ! [flg] Option conjunction is space character
       arg_val_srt_idx=3+opt_lng+1 ! [idx] Starting position of argument value
       if (dbg_lvl >= dbg_io) then
          write (6,'(8a)') prg_nm(1:ftn_strlsc(prg_nm)), &
               ': DEBUG '//sbr_nm//'() diassembles argument into Option = ''', &
               arg_val(3:2+opt_lng),''', conjunction is ''', &
               arg_val(3+opt_lng:3+opt_lng),''', argument = ''', &
               arg_val(arg_val_srt_idx:arg_lng),'''' ! Quote syntax is mystical
       endif                  ! endif dbg
    else
       opt_cnj_is_spc=.true.  ! [flg] Option conjunction is space character
       arg_val_srt_idx=1      ! [idx] Starting position of argument value
    endif                     ! endif
    if (opt_cnj_is_spc) then 
       ! Obtain value of option by getting next command-line argument
       ! New arguments will alter input values of arg_idx,arg_val
       call get_command_argument(arg_idx,arg_val)
       ! Increment counter for next read
       arg_idx=arg_idx+1
       ! Sanity checks
       arg_lng=len_trim(arg_val) ! [nbr] Length of argument
       if (arg_lng <= 0) stop 'ftn_???_arg_get() reports option lacks argument'
    endif                     ! endif opt_cnj_is_spc
    ! Read argument into double and return
    read (arg_val(arg_val_srt_idx:arg_lng),*) opt_val ! [frc] Variable to receive copy of arg_val
    if (dbg_lvl >= dbg_io) then
       write (6,'(2a,g13.6)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() assigned argument value = ',opt_val
    endif                     ! endif
    if(present(opt_flg)) opt_flg=.true. ! [flg] Variable set by command-line
    return
  end subroutine ftn_arg_get_dbl
  
  subroutine ftn_arg_get_flt( & ! [sbr] Process float-valued command-line argument
       arg_idx, & ! I/O [idx] Argument counter
       arg_val, & ! I/O [sng] Float to copy into opt_val
       opt_val, & ! O [frc] Variable to receive copy of arg_val
       opt_flg) ! O [frc] Variable set by command-line
    ! Purpose: Copy command-line float argument arg_val into variable opt_val
    ! Usage:
    ! call ftn_arg_get_flt(arg_idx,arg_val,opt_val) ! [sbr] Process float-valued command-line argument
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='ftn_arg_get_flt' ! [sng] Subroutine name
    ! Commons
    ! Input
    ! Input/Output
    integer,intent(inout)::arg_idx ! I/O [idx] Argument counter 
    character(len=*),intent(inout)::arg_val ! I/O [sng] Float to copy into opt_val
    ! Output
    real(selected_real_kind(p=6)),intent(out)::opt_val ! O [frc] Variable to receive copy of arg_val
    logical,optional,intent(out)::opt_flg ! O [flg] Variable set by command-line
    ! Local
    integer arg_lng           ! [nbr] Length of argument
    integer opt_lng           ! [nbr] Length of option
    integer arg_val_srt_idx   ! [idx] Starting position of argument value
    logical opt_cnj_is_spc ! [flg] Option conjunction is space character
    ! Main Code
    ! Print any diagnostics before current value of arg_val is overwritten
    arg_lng=len_trim(arg_val) ! [nbr] Length of argument
    opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
    if (dbg_lvl >= dbg_io) then
       write (6,'(2a,i2,3a,i2)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() reports arg_idx = ',arg_idx,  &
            ', Full option = ''',arg_val(1:arg_lng),''', Length = ',arg_lng
    endif                     ! endif dbg
    ! Determine whether format is --opt_sng=val or --opt_sng val
    ! Add two to cover preliminary dashes
    if (opt_lng+2 < arg_lng) then
       opt_cnj_is_spc=.false. ! [flg] Option conjunction is space character
       arg_val_srt_idx=3+opt_lng+1 ! [idx] Starting position of argument value
       if (dbg_lvl >= dbg_io) then
          write (6,'(8a)') prg_nm(1:ftn_strlsc(prg_nm)), &
               ': DEBUG '//sbr_nm//'() diassembles argument into Option = ''', &
               arg_val(3:2+opt_lng),''', conjunction is ''', &
               arg_val(3+opt_lng:3+opt_lng),''', argument = ''', &
               arg_val(arg_val_srt_idx:arg_lng),'''' ! Quote syntax is mystical
       endif                  ! endif dbg
    else
       opt_cnj_is_spc=.true.  ! [flg] Option conjunction is space character
       arg_val_srt_idx=1      ! [idx] Starting position of argument value
    endif                     ! endif
    if (opt_cnj_is_spc) then 
       ! Obtain value of option by getting next command-line argument
       ! New arguments will alter input values of arg_idx,arg_val
       call get_command_argument(arg_idx,arg_val)
       ! Increment counter for next read
       arg_idx=arg_idx+1
       ! Sanity checks
       arg_lng=len_trim(arg_val) ! [nbr] Length of argument
       if (arg_lng <= 0) stop 'ftn_???_arg_get() reports option lacks argument'
    endif                     ! endif opt_cnj_is_spc
    ! Read argument into float and return
    read (arg_val(arg_val_srt_idx:arg_lng),*) opt_val ! [frc] Variable to receive copy of arg_val
    if (dbg_lvl >= dbg_io) then
       write (6,'(2a,g13.6)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() assigned argument value = ',opt_val
    endif                     ! endif
    if(present(opt_flg)) opt_flg=.true. ! [flg] Variable set by command-line
    return
  end subroutine ftn_arg_get_flt
  
  subroutine ftn_arg_get_int( & ! [sbr] Process integer-valued command-line argument
       arg_idx, & ! I/O [idx] Argument counter
       arg_val, & ! I/O [sng] Integer to copy into opt_val
       opt_val, & ! O [nbr] Variable to receive copy of arg_val
       opt_flg) ! O [frc] Variable set by command-line
    ! Purpose: Copy command-line integer argument arg_val into variable opt_val
    ! Usage:
    ! call ftn_arg_get_int(arg_idx,arg_val,opt_val) ! [sbr] Process integer-valued command-line argument
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='ftn_arg_get_int' ! [sng] Subroutine name
    ! Commons
    ! Input
    ! Input/Output
    integer,intent(inout)::arg_idx ! I/O [idx] Argument counter 
    character(len=*),intent(inout)::arg_val ! I/O [sng] Integer to copy into opt_val
    ! 20030429: Change opt_val to inout since dbg_lvl itself may be an opt_val but is referenced on LHS
    integer,intent(inout)::opt_val ! O [nbr] Variable to receive copy of arg_val
    logical,optional,intent(out)::opt_flg ! O [flg] Variable set by command-line
    ! Output
    ! Local
    integer arg_lng           ! [nbr] Length of argument
    integer opt_lng           ! [nbr] Length of option
    integer arg_val_srt_idx   ! [idx] Starting position of argument value
    logical opt_cnj_is_spc ! [flg] Option conjunction is space character
    ! Main Code
    ! Print any diagnostics before current value of arg_val is overwritten
    arg_lng=len_trim(arg_val) ! [nbr] Length of argument
    opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
    if (dbg_lvl >= dbg_io) then
       write (6,'(2a,i2,3a,i2)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() reports arg_idx = ',arg_idx,  &
            ', Full option = ''',arg_val(1:arg_lng),''', Length = ',arg_lng
    endif                     ! endif dbg
    ! Determine whether format is --opt_sng=val or --opt_sng val
    ! Add two to cover preliminary dashes
    if (opt_lng+2 < arg_lng) then
       opt_cnj_is_spc=.false. ! [flg] Option conjunction is space character
       arg_val_srt_idx=3+opt_lng+1 ! [idx] Starting position of argument value
       if (dbg_lvl >= dbg_io) then
          write (6,'(8a)') prg_nm(1:ftn_strlsc(prg_nm)), &
               ': DEBUG '//sbr_nm//'() diassembles argument into Option = ''', &
               arg_val(3:2+opt_lng),''', conjunction is ''', &
               arg_val(3+opt_lng:3+opt_lng),''', argument = ''', &
               arg_val(arg_val_srt_idx:arg_lng),'''' ! Quote syntax is mystical
       endif                  ! endif dbg
    else
       opt_cnj_is_spc=.true.  ! [flg] Option conjunction is space character
       arg_val_srt_idx=1      ! [idx] Starting position of argument value
    endif                     ! endif
    if (opt_cnj_is_spc) then 
       ! Obtain value of option by getting next command-line argument
       ! New arguments will alter input values of arg_idx,arg_val
       call get_command_argument(arg_idx,arg_val)
       ! Increment counter for next read
       arg_idx=arg_idx+1
       ! Sanity checks
       arg_lng=len_trim(arg_val) ! [nbr] Length of argument
       if (arg_lng <= 0) stop 'ftn_arg_get_int() reports option lacks argument'
    endif                     ! endif opt_cnj_is_spc
    ! Read argument into integer and return
    read (arg_val(arg_val_srt_idx:arg_lng),*) opt_val ! [nbr] Variable to receive copy of arg_val
    if (dbg_lvl >= dbg_io) then
       write (6,'(2a,g13.6)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() assigned argument value = ',opt_val
    endif                     ! endif
    if(present(opt_flg)) opt_flg=.true. ! [flg] Variable set by command-line
    return
  end subroutine ftn_arg_get_int
  
  subroutine ftn_arg_get_lgc( & ! [sbr] Process logical-valued command-line argument
       arg_idx, & ! I/O [idx] Argument counter
       arg_val, & ! I/O [sng] Logical to copy into opt_val
       opt_val, & ! O [flg] Variable to receive copy of arg_val
       opt_flg) ! O [frc] Variable set by command-line
    ! Purpose: Copy command-line logical argument arg_val into variable opt_val
    ! Usage:
    ! call ftn_arg_get_lgc(arg_idx,arg_val,opt_val) ! [sbr] Process logical-valued command-line argument
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='ftn_arg_get_lgc' ! [sng] Subroutine name
    ! Commons
    ! Input
    ! Input/Output
    integer,intent(inout)::arg_idx ! I/O [idx] Argument counter 
    character(len=*),intent(inout)::arg_val ! I/O [sng] Logical to copy into opt_val
    ! Output
    logical,intent(out)::opt_val ! O [flg] Variable to receive copy of arg_val
    logical,optional,intent(out)::opt_flg ! O [flg] Variable set by command-line
    ! Local
    integer arg_lng           ! [nbr] Length of argument
    integer opt_lng           ! [nbr] Length of option
    integer arg_val_srt_idx   ! [idx] Starting position of argument value
    logical opt_cnj_is_spc ! [flg] Option conjunction is space character
    ! Main Code
    ! Print any diagnostics before current value of arg_val is overwritten
    arg_lng=len_trim(arg_val) ! [nbr] Length of argument
    opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
    if (dbg_lvl >= dbg_io) then
       write (6,'(2a,i2,3a,i2)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() reports arg_idx = ',arg_idx,  &
            ', Full option = ''',arg_val(1:arg_lng),''', Length = ',arg_lng
    endif                     ! endif dbg
    ! Determine whether format is --opt_sng=val or --opt_sng val
    ! Add two to cover preliminary dashes
    if (opt_lng+2 < arg_lng) then
       opt_cnj_is_spc=.false. ! [flg] Option conjunction is space character
       arg_val_srt_idx=3+opt_lng+1 ! [idx] Starting position of argument value
       if (dbg_lvl >= dbg_io) then
          write (6,'(8a)') prg_nm(1:ftn_strlsc(prg_nm)), &
               ': DEBUG '//sbr_nm//'() diassembles argument into Option = ''', &
               arg_val(3:2+opt_lng),''', conjunction is ''', &
               arg_val(3+opt_lng:3+opt_lng),''', argument = ''', &
               arg_val(arg_val_srt_idx:arg_lng),'''' ! Quote syntax is mystical
       endif                  ! endif dbg
    else
       opt_cnj_is_spc=.true.  ! [flg] Option conjunction is space character
       arg_val_srt_idx=1      ! [idx] Starting position of argument value
    endif                     ! endif
    if (opt_cnj_is_spc) then 
       ! Obtain value of option by getting next command-line argument
       ! New arguments will alter input values of arg_idx,arg_val
       call get_command_argument(arg_idx,arg_val)
       ! Increment counter for next read
       arg_idx=arg_idx+1
       ! Sanity checks
       arg_lng=len_trim(arg_val) ! [nbr] Length of argument
       if (arg_lng <= 0) stop 'ftn_arg_get_lgc() reports option lacks argument'
    endif                     ! endif opt_cnj_is_spc
    ! Read argument into integer and return
    read (arg_val(arg_val_srt_idx:arg_lng),*) opt_val ! [flg] Variable to receive copy of arg_val
    if (dbg_lvl >= dbg_io) then
       write (6,'(2a,g13.6)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() assigned argument value = ',opt_val
    endif                     ! endif
    if(present(opt_flg)) opt_flg=.true. ! [flg] Variable set by command-line
    return
  end subroutine ftn_arg_get_lgc
  
  subroutine ftn_arg_get_sng( & ! [sbr] Process string-valued command-line argument
       arg_idx, & ! I/O [idx] Argument counter
       arg_val, & ! I/O [sng] String to copy into opt_val
       opt_val, & ! O [sng] Variable to receive copy of arg_val
       opt_flg) ! O [frc] Variable set by command-line
    ! Purpose: Copy command-line string argument arg_val into variable opt_val
    ! opt_val is unchanged if arg_val cannot be copied
    ! Remainder of opt_val is NUL-initialized
    ! Usage:
    ! call ftn_arg_get_sng(arg_idx,arg_val,opt_val) ! [sbr] Process string-valued command-line argument
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='ftn_arg_get_sng' ! [sng] Subroutine name
    ! Commons
    ! Input
    ! Input/Output
    integer,intent(inout)::arg_idx ! I/O [idx] Argument counter 
    character(len=*),intent(inout)::arg_val ! I/O [sng] String to copy into opt_val
    ! Output
    character(len=*),intent(out)::opt_val ! O [sng] Variable to receive copy of arg_val
    logical,optional,intent(out)::opt_flg ! O [flg] Variable set by command-line
    ! Local
    integer idx
    integer len_var           ! [nbr] Length of string variable
    integer arg_lng           ! [nbr] Length of argument
    integer opt_lng           ! [nbr] Length of option
    integer arg_val_srt_idx   ! [idx] Starting position of argument value
    logical opt_cnj_is_spc ! [flg] Option conjunction is space character
    integer nll_srt_idx ! [idx] First NUL position in opt_val
    ! Main Code
    ! Print any diagnostics before current value of arg_val is overwritten
    len_var=len(opt_val)      ! [nbr] Length of string variable
    arg_lng=len_trim(arg_val) ! [nbr] Length of argument
    opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
    if (dbg_lvl >= dbg_io) then
       write (6,'(2a,i2,3a,i2)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() reports arg_idx = ',arg_idx,  &
            ', Option = ''',arg_val(1:arg_lng),''', Length = ',arg_lng
    endif                     ! endif
    ! Determine whether format is --opt_sng=val or --opt_sng val
    ! Add two to cover preliminary dashes
    if (opt_lng+2 < arg_lng) then
       opt_cnj_is_spc=.false. ! [flg] Option conjunction is space character
       arg_val_srt_idx=3+opt_lng+1 ! [idx] Starting position of argument value
       if (dbg_lvl >= dbg_io) then
          write (6,'(8a)') prg_nm(1:ftn_strlsc(prg_nm)), &
               ': DEBUG '//sbr_nm//'() diassembles argument into Option = ''', &
               arg_val(3:2+opt_lng),''', conjunction is ''', &
               arg_val(3+opt_lng:3+opt_lng),''', argument = ''', &
               arg_val(arg_val_srt_idx:arg_lng),'''' ! Quote syntax is mystical
       endif                  ! endif dbg
    else
       opt_cnj_is_spc=.true.  ! [flg] Option conjunction is space character
       arg_val_srt_idx=1      ! [idx] Starting position of argument value
    endif                     ! endif
    if (opt_cnj_is_spc) then 
       ! Obtain value of option by getting next command-line argument
       ! New arguments will alter input values of arg_idx,arg_val
       call get_command_argument(arg_idx,arg_val)
       ! Increment counter for next read
       arg_idx=arg_idx+1
       ! Sanity checks
       arg_lng=len_trim(arg_val) ! [nbr] Length of argument
       if (arg_lng <= 0) stop 'ftn_???_arg_get() reports option lacks argument'
       if (arg_lng > len_var) then
          write (6,'(2a,i5,a,i5)') prg_nm(1:ftn_strlsc(prg_nm)), &
               ': ERROR '//sbr_nm//'() reports argument length = ',arg_lng, &
               ' too long to fit into variable length ',len_var
          stop
       endif ! endif arg_lng > len_var
    endif                     ! endif opt_cnj_is_spc
    ! Read argument into string and return
    read (arg_val(arg_val_srt_idx:arg_lng),'(a)') opt_val ! [sng] Variable to receive copy of arg_val
    ! NUL-initialize remainder of opt_val
    nll_srt_idx=arg_lng-arg_val_srt_idx+2 ! [idx] First NUL position in opt_val
    do idx=nll_srt_idx,len_var
       opt_val(idx:idx)=char(0) ! [sng] Variable to receive copy of arg_val
    end do                    ! end loop over characters
    if (dbg_lvl >= dbg_io) then
       write (6,'(4a)') prg_nm(1:ftn_strlsc(prg_nm)), &
            ': DEBUG '//sbr_nm//'() assigned argument value = ',opt_val(1:len_trim(opt_val)),'<--End of opt_val string'
       call ftn_strprn(opt_val)   ! [fnc] Print character values of string
    endif ! endif
    if(present(opt_flg)) opt_flg=.true. ! [flg] Variable set by command-line
    return
  end subroutine ftn_arg_get_sng
  
  integer function ftn_strstr(sng1,sng2) ! [idx] Location of sng2 in sng1
    ! Purpose: Return the location of second string in the first
    ! Fortran intrinsic len(sng) is returned if string is not NUL-terminated
    ! Return -1 if sng1 does not contain sng1
    ! Prototype:
    ! integer ftn_strstr ! [idx] Location of sng1 in sng2
    implicit none
    ! Input
    character(len=*),intent(in)::sng1
    character(len=*),intent(in)::sng2
    ! Local
    integer lng1
    integer lng2
    integer idx1
    integer idx2
    ! Initialize
    idx2=0 ! CEWI for lf95
    ! Main Code
    lng1=ftn_strlen(sng1) ! [nbr] Length of string
    lng2=ftn_strlen(sng2) ! [nbr] Length of string
    do idx1=1,lng1
       if (sng1(idx1:idx1) == sng2(1:1)) then
          do idx2=2,lng2
             if (sng1(idx1+idx2-1:idx1+idx2-1) /= sng2(idx2:idx2)) goto 200
          end do              ! end loop over characters
          goto 100
       endif                  ! endif
200    continue
    end do                    ! end loop over characters
100 continue
    if (idx2 == lng2+1) then
       ftn_strstr=idx1
    else
       ftn_strstr=-1
    endif ! endif
    return
  end function ftn_strstr
  
  subroutine ftn_strini(sng) ! [sng] sng(1:len)=NUL
    ! Purpose: Initialize all elements of a character array to char(0)
    ! Usage:
    ! call ftn_strini(sng) ! [sng] sng(1:len)=NUL
    implicit none
    ! Input/Output
    character(len=*),intent(out)::sng ! I/O [sng] String to initialize
    ! Local
    integer idx
    integer lng
    ! Main Code
    lng=len(sng)
    do idx=1,lng
       sng(idx:idx)=char(0)
    end do                    ! end loop over characters
    ! write (6,*) 'ftn_strini(): Initialized string of length',lng
    return
  end subroutine ftn_strini
  
  subroutine ftn_strnul(sng) ! [sbr] NUL-initialize all characters after LSC
    ! Purpose: Change space characters to NUL characters in a string
    ! Spaces after last significant character (LSC) in string are changed to NUL
    ! Spaces before LSC are not affected
    ! Usage:
    ! call ftn_strnul(sng) ! [sbr] NUL-initialize all characters after LSC
    implicit none
    ! Input/Output
    character(len=*),intent(inout)::sng ! I/O [sng] String to NUL-initialize
    ! Local
    integer idx
    integer lng
    integer lst_sgn_chr_idx
    ! Main Code
    lng=ftn_strlen(sng) ! [nbr] Length of string
    lst_sgn_chr_idx=ftn_strlsc(sng)
    do idx=lst_sgn_chr_idx+1,lng
       if(sng(idx:idx) == ' ') sng(idx:idx)=char(0)
    end do                    ! end loop over characters
    return
  end subroutine ftn_strnul
  
  integer function ftn_strlsc(sng)
    ! Purpose: Return position of last significant character (LSC) in string
    ! LSC is last character that is not a space or a NUL character
    ! Routine is useful in determining effective lengths of strings which
    ! may have been padded by compiler or shell.
    ! Automatic padding is always done with spaces or NULs, so this routine
    ! tells us how long the non-padded portion of the string is.
    ! Routine will mis-handle user-defined strings that end in spaces
    ! LSC is always <= len(sng)
    ! Usage:
    ! psn=ftn_strlsc(sng) ! [idx] Last significant character position
    ! Prototype:
    ! integer ftn_strlsc ! [idx] Last significant character position
    implicit none
    ! Input
    character(len=*),intent(in)::sng
    ! Local
    integer idx
    integer lng
    ! Main Code
    lng=len(sng)
    do idx=lng,1,-1
       ! Look for last normal character
       if (sng(idx:idx) /= ' '.and.sng(idx:idx) /= char(0)) goto 100
    end do                    ! end loop over characters
100 continue
    ftn_strlsc=idx
    return
  end function ftn_strlsc
  
  integer function ftn_strfic(sng) ! [idx] First NUL or 8-bit character position
    ! Purpose: Return position of first insignificant character (FIC) in string
    ! FIC is first NUL character in string, or first 8-bit character, i.e., iachar > 127
    ! FIC is len(sng)+1 if string has no NUL character and no iachar > 127
    ! Usage:
    ! psn=ftn_strfic(sng) ! [idx] First NUL or 8-bit character position
    ! Prototype:
    ! integer ftn_strfic ! [idx] First NUL or 8-bit character position
    implicit none
    ! Input
    character(len=*),intent(in)::sng ! I [sng] String to examine
    ! Local
    integer idx
    integer lng
    ! Main Code
    lng=len(sng)
    do idx=1,lng
       ! Look for first abnormal character
       if (sng(idx:idx) == char(0) .or. iachar(sng(idx:idx)) > 127) goto 100
    end do                    ! end loop over characters
100 continue
    ftn_strfic=idx ! [idx] First NUL or 8-bit character position
    return
  end function ftn_strfic
  
  subroutine ftn_strprn(sng) ! [fnc] Print character values of string
    ! Purpose: Print character values of string
    ! Usage:
    ! call ftn_strprn(sng) ! [fnc] Print character values of string
    implicit none
    ! Input
    character(len=*),intent(in)::sng
    ! Local
    character(1) chr_crr ! [sng] Current character
    integer idx
    integer lng
    integer iachr_crr ! [enm] Integer representing current character
    ! Main Code
    lng=len(sng)
    do idx=1,lng
       ! Print each normal character
       chr_crr=sng(idx:idx) ! [sng] Current character
       iachr_crr=iachar(chr_crr) ! [enm] Integer representing current character
       if (idx < lng) then
          if (iachr_crr == 0) then
             write (6,'(a3,a3,i3,a2)',advance="no") 'NUL',' = ',iachr_crr,', '
          else if (iachr_crr == 32) then
             write (6,'(a5,a3,i3,a2)',advance="no") 'SPACE',' = ',iachr_crr,', '
          else
             write (6,'(a1,a3,i3,a2)',advance="no") chr_crr,' = ',iachr_crr,', '
          endif               ! endif
       else
          if (iachr_crr == 0) then
             write (6,'(a3,a3,i3)') 'NUL',' = ',iachr_crr,', '
          else if (iachr_crr == 32) then
             write (6,'(a5,a3,i3)') 'SPACE',' = ',iachr_crr,', '
          else
             write (6,'(a1,a3,i3)') chr_crr,' = ',iachr_crr,', '
          endif               ! endif
       endif                  ! endif
    end do                    ! end loop over characters
    return
  end subroutine ftn_strprn
  
  subroutine ftn_strcpy(sng1,sng2)
    ! Purpose: Copy sng2 into sng1
    ! Space remaining at end of sng1 is NUL-initialized
    ! This function works well at copying fixed strings into variables, e.g.,
    ! Usage:
    ! call ftn_strcpy(sng,'Hello World') ! [fnc] Copy sng2 into sng1, NUL-terminate unused space
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Input
    character(len=*),intent(in)::sng2
    ! Input/Output
    character(len=*),intent(inout)::sng1
    ! Local
    integer len1
    integer len2
    integer lng1
    integer lng2
    ! Main Code
    len1=len(sng1)
    len2=len(sng2)
    ! Initialize sng1 before diagnostic ftn_strlen(sng1) is evaluated
    ! This handles situation when sng1 is large enough, but not-initialized in calling routine
    ! Otherwise some compilers, e.g., lf95 in debug mode, may complain in ftn_strlen(sng1)->ftn_strfic()
    ! when they are unable to determine string length
    ! Initialize sng1 before copying sng2 into sng1 WLOG since ftn_strcpy() behavior
    ! is intended to mimic C library strcpy() which NUL-terminates the copy, and thus
    ! all material after the NUL may be assumed to be undefined, i.e., set to NUL.
    ! This is not equivalent to calling ftn_strnul(sng1) after the copy since when sng1 is undefined
    ! There could be corner cases where contents of undefined sng1 character immediately
    ! following copy of sng2 is not NUL or 8-bit. Results are unpredicatable in that case.
    ! Hence initialize sng1 to NUL before using ftn_strlen(sng1) and before copying sng2
    call ftn_strini(sng1) ! [sng] sng(1:len)=NUL
    lng1=ftn_strlen(sng1) ! [nbr] Length of string
    lng2=ftn_strlen(sng2) ! [nbr] Length of string
    if (len1 < lng2) then
       write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR len1 < lng2 in ftn_strcpy()'
       write (6,'(a,i3,a,i3,2a)') 'len1 = ',len1,', lng1 = ',lng1,', sng1 = ',sng1
       write (6,'(a,i3,a,i3,2a)') 'len2 = ',len2,', lng2 = ',lng2,', sng2 = ',sng2
       stop 'EXIT_FAILURE from ftn_strcpy()'
    endif                     ! endif
    sng1(1:lng2)=sng2(1:lng2)
    return
  end subroutine ftn_strcpy
  
  integer function ftn_strcmp(sng1,sng2) ! [fnc] Compare sng1 to sng2: -1,0,1 iff sng1 <,=,> sng2
    ! Purpose: Compare sng2 to sng1
    ! Returns an integer less than, equal to, or greater than
    ! zero  if  sng1  is  found, respectively, to be less than, to
    ! match, or be greater than sng2.
    ! Usage:
    ! if(ftn_strcmp(sng1,"Hello World")==0) then ! [fnc] Compare sng1 to sng2: -1,0,1 iff sng1 <,=,> sng2
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Input
    character(len=*),intent(in)::sng1
    character(len=*),intent(in)::sng2
    ! Input/Output
    ! Local
    integer idx               ! [idx] Counting idx
    integer len1
    integer len2
    integer lng1
    integer lng2
    integer lng_min
    ! Main Code
    len1=len(sng1)
    len2=len(sng2)
    lng1=ftn_strlen(sng1) ! [nbr] Length of string
    lng2=ftn_strlen(sng2) ! [nbr] Length of string
    if (.false.) then 
       ! If this debugging block is turned on then attempting to print,ftn_strcmp()
       ! will cause a recursive I/O error at runtime
       write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),': DEBUG Diagnostics from ftn_strcmp()'
       write (6,'(2(a,i3),2a)') 'len1 = ',len1,', lng1 = ',lng1,', sng1 = ',sng1
       write (6,'(2(a,i3),2a)') 'len2 = ',len2,', lng2 = ',lng2,', sng2 = ',sng2
       call ftn_strprn(sng1) ! [fnc] Print character values of string
       call ftn_strprn(sng2) ! [fnc] Print character values of string
    endif                     ! endif
    lng_min=min(lng1,lng2)
    do idx=1,lng_min
       ! Prevent memory overruns
       ! Characters are equal, proceed to next character
       if (iachar(sng1(idx:idx)) == iachar(sng2(idx:idx))) cycle ! Goto next iteration
       ! Characters are unequal, determine relationship
       if (iachar(sng1(idx:idx)) < iachar(sng2(idx:idx))) then
          ftn_strcmp=-1    ! [fnc] Compare sng1 to sng2: -1,0,1 iff sng1 <,=,> sng2
       else                ! (iachar(sng1(idx:idx)) > iachar(sng2(idx:idx)))
          ftn_strcmp=1     ! [fnc] Compare sng1 to sng2: -1,0,1 iff sng1 <,=,> sng2
       endif               ! endif
       
       ! Characters were unequal, exit with appropriate value
       goto 100            ! Exit Branch
    end do
    
    if (lng1 == lng2) then
       ! All significant values of both strings were equal in all positions
       ! Thus the strings are equal by our criterion
       ftn_strcmp=0           ! [fnc] Compare sng1 to sng2: -1,0,1 iff sng1 <,=,> sng2
    else
       ! All significant values of both strings were equal in all mutually legal positions
       ! However, one string has more significant characters than the other
       ! Thus return value depends on which string is longer
       if (lng1 > lng2) then
          ftn_strcmp=1        ! [fnc] Compare sng1 to sng2: -1,0,1 iff sng1 <,=,> sng2
       else 
          ftn_strcmp=-1       ! [fnc] Compare sng1 to sng2: -1,0,1 iff sng1 <,=,> sng2
       endif                  ! endif
    endif ! endif lng1 == lng2
    
100 continue                  ! Exit branch for unequal values inside loop
    return
  end function ftn_strcmp
  
  subroutine ftn_strcpylsc(sng1,sng2)
    ! Purpose: Copy sng2 into sng1
    ! Space remaining at end of sng1 is NUL-initialized
    ! This function works well at copying fixed strings into variables, e.g.,
    ! Usage:
    ! call ftn_strcpylsc(sng,'Hello World') ! [fnc] Copy up to LSC of sng2 into sng1, NUL-initialize unused space
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Input
    character(len=*),intent(in)::sng2
    ! Input/Output
    character(len=*),intent(inout)::sng1
    ! Local
    integer idx               ! [idx] Counting idx
    integer len1
    integer len2
    integer lsc1
    integer lsc2
    ! Main Code
    len1=len(sng1)
    len2=len(sng2)
    lsc1=ftn_strlsc(sng1)
    lsc2=ftn_strlsc(sng2)
    if (len1 < lsc2) then
       write (6,'(a,a)') prg_nm(1:ftn_strlsc(prg_nm)),': ERROR len1 < lsc2 in ftn_strcpylsc()'
       write (6,'(a,i3,a,i3,2a)') 'len1 = ',len1,', lsc1 = ',lsc1,', sng1 = ',sng1
       write (6,'(a,i3,a,i3,2a)') 'len2 = ',len2,', lsc2 = ',lsc2,', sng2 = ',sng2
       stop 'EXIT_FAILURE from ftn_strcpylsc()'
    endif                     ! endif
    sng1(1:lsc2)=sng2(1:lsc2)
    ! NUL-initialize rest of string
    do idx=lsc2+1,len1
       sng1(idx:idx)=char(0)
    end do                    ! end loop over characters
    return
  end subroutine ftn_strcpylsc
  
  subroutine ftn_strcat( &   ! [fnc] sng1 := sng1 // sng2
       sng1, &                ! I/O [sng] String to affix second string to
       sng2)                ! I [sng] String to affix to first string
    ! Purpose: sng1 := sng1 // sng2
    ! The case where sng1=sng2 is handled correctly
    ! Usage:
    ! call ftn_strcat(sng1,sng2) ! [fnc] sng1 := sng1 // sng2
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Input
    character(len=*),intent(in)::sng2 ! I [sng] String to affix to first string
    ! Input/Output
    character(len=*),intent(inout)::sng1 ! I/O [sng] String to affix second string to
    ! Local
    integer len1
    integer len2
    integer lng1
    integer lng2
    ! Main Code
    len1=len(sng1)
    len2=len(sng2)
    lng1=ftn_strlen(sng1) ! [nbr] Length of string
    lng2=ftn_strlen(sng2) ! [nbr] Length of string
    if (lng1+lng2 >= len1) then
       write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR lng1+lng2 >= len1 in ftn_strcat()'
       write (6,'(a,i3,a,i3,2a)') 'len1 = ',len1,', lng1 = ',lng1,', sng1 = ',sng1
       write (6,'(a,i3,a,i3,2a)') 'len2 = ',len2,', lng2 = ',lng2,', sng2 = ',sng2
       stop 'EXIT_FAILURE from ftn_strcat()'
    endif                     ! endif
    sng1(lng1+1:lng1+lng2)=sng2(1:lng2) ! I/O [sng] String to affix second string to
    return
  end subroutine ftn_strcat
  
  subroutine ftn_strpfx(sng1,sng2) ! [sbr] sng2 := sng1 // sng2
    ! Purpose: sng2 := sng1 // sng2
    ! Differences with ftn_strcat(): 
    ! 1. Result is stored in second string not first string
    ! 2. Strings are trimmed before concatenation
    ! The case where sng1=sng2 is handled correctly
    ! Usage:
    ! call ftn_strpfx(sng1,sng2) ! [sbr] sng2 := sng1 // sng2
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Input
    character(len=*),intent(in)::sng1
    ! Output
    character(len=*),intent(inout)::sng2
    ! Local
    integer len1
    integer len2
    integer len_trim1
    integer len_trim2
    integer lng1
    integer lng2
    integer lsc1
    integer lsc2
    ! Main Code
    len1=len(sng1)
    len2=len(sng2)
    lng1=ftn_strlen(sng1) ! [nbr] Length of string
    lng2=ftn_strlen(sng2) ! [nbr] Length of string
    lsc1=ftn_strlsc(sng1)
    lsc2=ftn_strlsc(sng2)
    len_trim1=len_trim(sng1)
    len_trim2=len_trim(sng2)
    if (lsc1+lsc2 >= len2) then
       write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR lsc1+lsc2 >= len2 in ftn_strpfx()'
       write (6,'(4(a,i3),2a)') 'len1 = ',len1,', len_trim1 = ',len_trim1,', lng1 = ',lng1,', sng1 = ',sng1
       write (6,'(4(a,i3),2a)') 'len2 = ',len2,', len_trim2 = ',len_trim2,', lng2 = ',lng2,', sng2 = ',sng2
       call ftn_strprn(sng1)
       call ftn_strprn(sng2)
       stop 'EXIT_FAILURE from ftn_strpfx()'
    endif                     ! endif
    ! Order is important
    sng2(lsc1+1:lsc1+lsc2)=sng2(1:lsc2)
    sng2(1:lsc1)=sng1(1:lsc1)
    return
  end subroutine ftn_strpfx
  
  subroutine ftn_drcpfx( & ! [sbr] fl_nm := drc/fl_nm 
       drc, & ! I [sng] Directory to prepend
       fl_nm) ! I/O [sng] File name
    ! Purpose: fl_nm := drc/fl_nm, more or less
    ! Differences with ftn_strpfx(): 
    ! 1. Result is stored in second string not first string
    ! 2. Strings are trimmed before concatenation
    ! 3. Filenames containing slashes are unaltered
    ! Usage:
    ! call ftn_drcpfx(drc,fl_nm) ! [sbr] fl_nm := drc/fl_nm
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Input
    character(len=*),intent(in)::drc ! I [sng] Directory to prepend
    ! Input/Output
    character(len=*),intent(inout)::fl_nm ! I/O [sbr] File name
    ! Local
    integer len_drc
    integer len_fl
    integer lng_drc
    integer lng_fl
    ! Main Code
    len_drc=len(drc)
    len_fl=len(fl_nm)
    lng_drc=ftn_strlen(drc) ! [nbr] Length of string
    lng_fl=ftn_strlen(fl_nm) ! [nbr] Length of string
    if (lng_drc+lng_fl >= len_fl) then
       write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR lng_drc+lng_fl >= len_fl in ftn_drcpfx()'
       write (6,'(2(a,i3),2a)') 'len_drc = ',len_drc,', lng_drc = ',lng_drc,', drc = ',drc
       write (6,'(2(a,i3),2a)') 'len_fl = ',len_fl,', lng_fl = ',lng_fl,', fl_nm = ',fl_nm
       call ftn_strprn(drc)
       call ftn_strprn(fl_nm)
       stop 'EXIT_FAILURE from ftn_drcpfx()'
    endif                     ! endif
    
    if ( &                     ! Conditions for prepending fl_nm with drc:
         lng_drc > 0 .and. &   ! Directory must be non-NUL 
         ftn_strstr(fl_nm,'/') == -1 & ! fl_nm does not look like directory itself
         ) then
       ! If directory does not have trailing slash...
       if (drc(lng_drc:lng_drc) /= '/') then ! Directory does not have trailing slash
          ! ...and if room exists for one...
          if (lng_drc+lng_fl+1 < len_fl) then
             ! ...then add one
             ! Order is important
             fl_nm(lng_drc+2:lng_drc+2+lng_fl)=fl_nm(1:lng_fl)
             fl_nm(1:lng_drc)=drc(1:lng_drc)
             fl_nm(lng_drc+1:lng_drc+1)='/'
          endif               ! endif adding trailing slash
       else                   ! Directory does have trailing slash
          fl_nm(lng_drc+1:lng_drc+1+lng_fl)=fl_nm(1:lng_fl)
          fl_nm(1:lng_drc)=drc(1:lng_drc)
       endif                  ! Directory does have trailing slash
    endif                     ! endif prepending drc
    call ftn_strnul(fl_nm) ! [sbr] NUL-initialize all characters after LSC
    
    return
  end subroutine ftn_drcpfx
  
  subroutine ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
    ! Purpose: Create identity string from CVS input
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Parameter
    integer,parameter::CVS_kk=0
    integer,parameter::CVS_kv=1
    integer,parameter::CVS_kkv=2
    ! Input
    character(len=*),intent(in)::CVS_Id
    character(len=*),intent(in)::CVS_Revision
    character(len=*),intent(in)::CVS_Date
    ! Output
    character(len=*),intent(out)::prg_ID
    ! Local
    integer CVS_typ
    integer Date_ptr
    integer slash_ptr
    integer prg_ptr
    integer vrs_ptr
    ! Main Code
    ! Decide how CVS keywords were expanded
    Date_ptr=ftn_strstr(CVS_Date,'Date')
    if (Date_ptr < 0) then
       ! String 'Date' was not found---CVS expansion must be -kv
       CVS_typ=CVS_kv
    else
       ! String 'Date' was found---CVS expansion is -kk or -kkv
       slash_ptr=ftn_strstr(CVS_Date,'/')
       if (slash_ptr > 0) then 
          CVS_typ=CVS_kkv 
       else 
          CVS_typ=CVS_kk
       endif                  ! endif
    endif                     ! endif
    
    if (.false.) then
       write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': DEBUG ftn_prg_ID_mk() reports CVS strings:'
       write (6,'(a,a)') 'CVS_Id: ',CVS_Id
       write (6,'(a,a)') 'CVS_Revision: ',CVS_Revision
       write (6,'(a,a)') 'CVS_Date: ',CVS_Date
       write (6,'(a,i1)') 'CVS_typ: ',CVS_typ
       write (6,'(a,i2)') 'Date_ptr: ',Date_ptr
       write (6,'(a,i2)') 'slash_ptr: ',slash_ptr
    endif
    
    ! CVS_Revision and CVS_Date are parameters---do not attempt to write them
    !  call ftn_strnul(CVS_Revision) ! [sbr] NUL-initialize all characters after LSC
    !  call ftn_strnul(CVS_Date) ! [sbr] NUL-initialize all characters after LSC
    call ftn_strini(prg_ID) ! [sng] sng(1:len)=NUL
    if (CVS_typ == CVS_kk) then
       prg_ID='Source file unknown Version unknown Date unknown' // char(0)
    else 
       prg_ptr=ftn_strstr(CVS_Id,',v') ! ',v' is right after program name
       if (prg_ptr < 0) stop 'ERROR: ftn_prg_ID_mk() unable to find source code name'
       if (CVS_typ == CVS_kv) then
          prg_ID=CVS_Id(1:prg_ptr-1) //  &
               ' version ' // CVS_Revision(1:ftn_strlen(CVS_Revision)) //  &
               ' dated ' // CVS_Date(1:ftn_strlen(CVS_Date)) // ' GMT' // char(0)
       else if (CVS_typ == CVS_kkv) then
          vrs_ptr=ftn_strstr(CVS_Revision,'ion:') ! 'ion:' should appear two characters before version
          if (vrs_ptr < 0) stop 'ERROR: ftn_prg_ID_mk() unable to find revision number'
          prg_ID=CVS_Id(6:prg_ptr-1) //  &
               ' version ' // CVS_Revision(vrs_ptr+5:ftn_strlen(CVS_Revision)-2) //  &
               ' dated ' // CVS_Date(8:26) // ' GMT' // char(0)
       else 
          write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR ftn_prg_ID_mk() reports unknown CVS_typ'
       endif                  ! endif
    endif                     ! endif
    return
  end subroutine ftn_prg_ID_mk
  
  subroutine ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
    ! Purpose: Return a copy of command-line and initialize program name
    ! See http://www.winteracter.com/f2kcli/index.htm for alternatives
    ! Usage:
    ! call ftn_cmd_ln_sng(cmd_ln) ! [sng] command-line
    use mod_utils ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Commons
    ! Output
    character(len=*),intent(out)::cmd_ln ! O [sng] command-line
    ! Local
    character(len=80) arg_val         ! [sng] command-line argument value
    integer arg_idx           ! [idx] Counting index
    integer arg_nbr           ! [nbr] Number of command-line arguments
    integer cmd_ln_len        ! [nbr] Length of command-line
    integer cmd_ln_lng        ! [nbr] Length of command-line
    ! Main Code
    ! Initialize defaults
    call ftn_strini(cmd_ln) ! [sng] sng(1:len)=NUL
    call ftn_strini(prg_nm) ! [sng] sng(1:len)=NUL
    cmd_ln_len=len(cmd_ln)  ! [nbr] Length of command-line
    cmd_ln_lng=0            ! [nbr] Length of command-line CEWI
    
    arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments
    if (dbg_lvl >= dbg_vec) then
       write (6,'(a,i2)') 'ftn_cmd_ln_sng(): arg_nbr = ',arg_nbr
    endif                     ! endif dbg
    
    ! Loop over arguments
    do arg_idx=0,arg_nbr      ! NB: Loop starts with 0
       ! Argument 0 is program name
       ! if (arg_idx > 0) call ftn_strcat(cmd_ln,' ')
       ! Insert space between arguments
       if (arg_idx > 0) cmd_ln(cmd_ln_lng+1:cmd_ln_lng+1)=' '
       call ftn_strini(arg_val) ! [sng] sng(1:len)=NUL
       call get_command_argument(arg_idx,arg_val)
       call ftn_strnul(arg_val) ! [sbr] NUL-initialize all characters after LSC
       call ftn_strcat(cmd_ln,arg_val)
       if (arg_idx == 0) call ftn_strcpy(prg_nm,arg_val) ! [sng] Program name
       cmd_ln_lng=ftn_strlen(cmd_ln) ! [nbr] Length of string
       if (cmd_ln_lng > cmd_ln_len) stop 'cmd_ln_lng > cmd_ln_len in ftn_cmd_ln_sng()'
       if (dbg_lvl >= dbg_vrb) then
          write (6,'(a,i3)') 'arg_idx = ',arg_idx
          write (6,'(a,a)') 'arg_val = ',arg_val
          write (6,'(a,i3)') 'len(arg_val) = ',len(arg_val)
          write (6,'(a,i3)') 'ftn_strlen(arg_val) = ',ftn_strlen(arg_val)
          write (6,'(a,i3)') 'ftn_strlen(cmd_ln) = ',ftn_strlen(cmd_ln)
          write (6,'(a,a)') 'cmd_ln = ',cmd_ln
       endif                  ! endif dbg
    end do                     ! end loop over arg
    
    if (dbg_lvl >= dbg_io) write (6,'(a,a)') 'cmd_ln = ',cmd_ln(1:ftn_strlen(cmd_ln))
    return
  end subroutine ftn_cmd_ln_sng
  
  character(len=10) function ftn_date2sng(idate)
    ! Purpose: Convert integer date in YYYYMMDD to a character array in the form
    ! 'YYYY-MM-DD'
    ! Author: Brian Eaton charutl.F:idate2char()
    ! Recoded: 19990901 Charlie Zender
    ! Usage:
    ! character(len=10) ftn_date2sng   ! [sng] Convert YYYYMMDD integer to YYYY-MM-DD string
    implicit none
    ! Input
    integer,intent(in)::idate
    ! Local
    integer yy 
    integer mm
    integer dd
    character dash
    ! Main code
    dash='-'
    ! Extract year, month, and day from date
    yy=idate/10000
    mm=mod(abs(idate),10000)/100
    dd=mod(abs(idate),100)
    write(ftn_date2sng,'(i4.4,a1,i2.2,a1,i2.2)') yy,dash,mm,dash,dd
    return
  end function ftn_date2sng
  
  character(len=8) function ftn_sec2sng(isec)
    ! Purpose: Convert integer seconds to a character array in the form 'HH:MM:SS'
    ! Author: Brian Eaton charutl.F:isec2char()
    ! Recoded: 19990901 Charlie Zender
    ! character(len=8) ftn_sec2sng   ! [sng] Convert integer to 'HH:MM:SS' string
    implicit none
    ! Input
    integer,intent(in)::isec
    ! Local
    integer hh
    integer mm
    integer ss
    character colon
    ! Main Code
    colon=':'
    ! Extract hour, minutes, seconds
    hh=isec/3600
    mm=mod(isec,3600)/60
    ss=mod(isec,60)
    write(ftn_sec2sng,'(i2.2,a1,i2.2,a1,i2.2)') hh,colon,mm,colon,ss
    return
  end function ftn_sec2sng
  
end module mod_sng ! [mdl] String manipulation
