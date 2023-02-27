module mod_sms2grid
  use all_vars
  use mod_utils
  use mod_input
  use mod_prec
  implicit none


  CHARACTER(LEN=80) ::FNAME


  CHARACTER(LEN=80) ::GRID_FILE_IN
  CHARACTER(LEN=80) ::DEPTH_FILE_IN
  CHARACTER(LEN=80) ::GRID_FILE_OUT
  CHARACTER(LEN=80) ::DEPTH_FILE_OUT
  CHARACTER(LEN=80) ::CORIOLIS_FILE_OUT

  CHARACTER(LEN=80) :: GRID_INPUT_UNITS
  CHARACTER(LEN=80) :: DEPTH_INPUT_UNITS
  CHARACTER(LEN=80) :: GRID_OUTPUT_UNITS

  REAL(SP) :: DEPTH_MINIMUM


  NAMELIST /NML_SMS2GRID/                      &
       & INPUT_DIR,                            &
       & OUTPUT_DIR,                           &
       & GRID_FILE_IN,                         &
       & DEPTH_FILE_IN,                        &
       & DEPTH_INPUT_UNITS,                    &
       & DEPTH_MINIMUM,                        &
       & GRID_FILE_OUT,                        &
       & DEPTH_FILE_OUT,                       &
       & CORIOLIS_FILE_OUT,                    &
       & GRID_INPUT_UNITS,                     &
       & GRID_OUTPUT_UNITS,                    &
       & PROJECTION_REFERENCE



  INTEGER, PARAMETER :: GRDINUNIT = 101
  INTEGER, PARAMETER :: DPTHINUNIT = 102

  INTEGER, PARAMETER :: GRDOUTUNIT = 103
  INTEGER, PARAMETER :: DPTHOUTUNIT = 104
  INTEGER, PARAMETER :: COROUTUNIT = 105

  INTEGER, PARAMETER :: NODESTRINGUNIT = 106
  LOGICAL ::  NSOPEN=.false.
  INTEGER :: NSCNT,NNS
  CHARACTER(LEN=80) :: NODESTRINGNAME


  LOGICAL :: DEPTH_ON, COR_ON, PROJ_ON, DEPTH_IN_GRID


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

          else if (opt_sng == "Fileame" .or.opt_sng == "filename"&
               & .or.opt_sng == "FILENAME") then

             call ftn_arg_get(arg_idx,arg_val,FName) ! [sng] Input file
             FName=FName(1:ftn_strlen(FName))
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

    WRITE(IPT,*) "THIS PROGRAM CONVERTS SMS FILES INTO FVCOM RUN FILES!"
    WRITE(IPT,*) "!    You can convert the grid file"
    WRITE(IPT,*) "!    Extract SMS nodestrings for OBC files"
    WRITE(IPT,*) "!    Create a coriolis file"
    WRITE(IPT,*) "!    Convert a depth file"
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    Any of the above can change coordinates"
    WRITE(IPT,*) "!    from meters => degrees or"
    WRITE(IPT,*) "!    degrees => meters"
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "! OPTIONS:"
    WRITE(IPT,*) "! --filename=XXX"
    WRITE(IPT,*) "!   The namelist runfile for the program! "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   Namelist OPTIONS: "
    WRITE(IPT,*) "!    INPUT_DIR"
    WRITE(IPT,*) "!    OUTPUT_DIR"
    WRITE(IPT,*) "!    GRID_FILE_IN"
    WRITE(IPT,*) "!    DEPTH_FILE_IN (OPTIONAL)"
    WRITE(IPT,*) "!    DEPTH_INPUT_UNITS (OPTIONAL*)"
    WRITE(IPT,*) "!    DEPTH_MINIUMUM (OPTIONAL)"
    WRITE(IPT,*) "!    GRID_FILE_OUT"
    WRITE(IPT,*) "!    DEPTH_FILE_OUT (OPTIONAL*)"
    WRITE(IPT,*) "!    CORIOLIS_FILE_OUT (OPTIONAL)"
    WRITE(IPT,*) "!    GRID_FILE_UNITS (OPTIONAL*)"
    WRITE(IPT,*) "!    PROJECTION_REFERENCE (OPTIONAL*)"
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    EXAMPLE NAMELIST:"
    write(UNIT=IPT,NML=NML_SMS2GRID)


    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT

  SUBROUTINE INITIALIZE_NML
    IMPLICIT NONE

    INPUT_DIR            = "/my/input/directory"
    OUTPUT_DIR           = "/my/input/directory"
    GRID_FILE_IN         = "default"
    DEPTH_FILE_IN        = "default"
    DEPTH_INPUT_UNITS     = "default"
    DEPTH_MINIMUM         = -999.9_SP
    GRID_FILE_OUT        = "default"
    DEPTH_FILE_OUT       = "default"
    CORIOLIS_FILE_OUT    = "default"
    GRID_INPUT_UNITS     = "default"
    GRID_OUTPUT_UNITS    = "default"
    PROJECTION_REFERENCE = "none"
    

  END SUBROUTINE INITIALIZE_NML

  SUBROUTINE READ_NAMELIST
    IMPLICIT NONE
    integer :: ios, i
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: Read_Name_List;"


    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)


    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read IO Information
    READ(UNIT=NMLUNIT, NML=NML_SMS2GRID,IOSTAT=ios)
    if(ios .NE. 0 )THEN
       write(UNIT=IPT,NML=NML_SMS2GRID)
       
       CALL FATAL_ERROR &
         &("Can Not Read NameList NML_SMS2GRID from file: "//trim(FNAME))
    END if
    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_SMS2GRID)

    CLOSE(NMLUNIT)


  END SUBROUTINE READ_NAMELIST


  SUBROUTINE OPEN_FILES
    IMPLICIT NONE

    integer :: ios
    logical :: fexist
    character(len=160) :: pathnfile
    
    !Open Grid File:
    IF(GRID_FILE_IN == 'default' .or.  GRID_FILE_OUT=='default') CALL FATAL_ERROR&
            &('GRID FILE NAME HAS BEEN SET FOR ONLY ONE OF INPUT AND OUTPUT?',&
            & 'GRID_FILE_IN:'//TRIM(GRID_FILE_IN),&
            & 'GRID_FILE_OUT:'//TRIM(GRID_FILE_OUT))

    pathnfile = trim(INPUT_DIR)//trim(GRID_FILE_IN)
    Call FOPEN(GRDINUNIT,trim(pathnfile),'cfr')

    pathnfile = trim(OUTPUT_DIR)//trim(GRID_FILE_OUT)
    Call FOPEN(GRDOUTUNIT,trim(pathnfile),'ofr')

    

    PROJ_ON = (GRID_INPUT_UNITS == 'degrees' .or. GRID_INPUT_UNITS == 'meters')

    IF(PROJ_ON) PROJ_ON = (GRID_OUTPUT_UNITS == 'degrees' .or. GRID_OUTPUT_UNITS == 'meters')

    IF(PROJ_ON) THEN
# if !defined (PROJ)
        Call Fatal_error&
            & ("Proj is not compiled - please recompile with PROJ (see make.inc)")
# endif
       IF(.not. HAVE_PROJ(PROJECTION_REFERENCE) ) Call Fatal_error&
            & ("Proj is not working properly, Please fix your PROJECTION_REFERENCE string")
    END IF
    
    ! DEPTH FILE
    DEPTH_ON = DEPTH_FILE_IN /= 'default' .and. DEPTH_FILE_OUT /= 'default'
    IF(DEPTH_ON) THEN

       IF(PROJ_ON) THEN
          
          IF(DEPTH_INPUT_UNITS == 'default') CALL FATAL_ERROR&
               &("Please specify the coordinates units in the Depth file!")
       END IF

       IF(DEPTH_FILE_IN == GRID_FILE_IN) THEN 
          DEPTH_IN_GRID = .true.

       else
          pathnfile = trim(INPUT_DIR)//trim(DEPTH_FILE_IN)
          Call FOPEN(DPTHINUNIT,trim(pathnfile),'cfr')
       end IF

       pathnfile = trim(OUTPUT_DIR)//trim(DEPTH_FILE_OUT)
       Call FOPEN(DPTHOUTUNIT,trim(pathnfile),'ofr')
       
    ELSE
       ! CHECK FOR ERROR
       IF(DEPTH_FILE_IN /=  DEPTH_FILE_OUT) CALL FATAL_ERROR&
            &('DEPTH FILE NAME HAS BEEN SET FOR ONLY ONE OF INPUT AND OUTPUT?',&
            & 'DEPTH_FILE_IN:'//TRIM(DEPTH_FILE_IN),&
            & 'DEPTH_FILE_OUT:'//TRIM(DEPTH_FILE_OUT))
       
    END IF



    COR_ON = CORIOLIS_FILE_OUT /= 'default'
    IF(COR_ON .and. .not. PROJ_ON) CALL FATAL_ERROR&
         &("IF YOU WANT TO MAKE YOUR CORRIOLIS",&
         & "FILE PLEASE COMPILE AND RUN WITH PROJ",&
         & "You must specify, the units of in put and output...")

    pathnfile = trim(OUTPUT_DIR)//trim(CORIOLIS_FILE_OUT)
    Call FOPEN(COROUTUNIT,trim(pathnfile),'ofr')



  END SUBROUTINE OPEN_FILES


  SUBROUTINE CONVERT_SMS2DAT
    IMPLICIT NONE
    
    real(Dp) :: ri1, ri2, ri3, ri4
    real(Dp) :: ro1, ro2, ro3, ro4, ro5, cc1, cc2

    integer :: ii1, ii2, ii3, ii4
    integer :: io1, io2, io3, io4
    
    CHARACTER(len=20), ALLOCATABLE :: lsplit(:)
    CHARACTER(LEN=120) :: line,pathnfile
    CHARACTER(LEN=20) :: CTEMP
    
    integer :: ios, Nodes, Cells, I,sizel,MYOS

    ! READ AND CONVERT THE GRID FILE

    IOS = 0
    nodes = 0 
    cells = 0
    NSCNT = 0
    NNS = 0

    ! COUNT THE NUMBER OF NODES AND ELEMENTS


    DO WHILE(IOS == 0)
       READ(GRDINUNIT,'(a)',IOSTAT=IOS) line
       
       IF(line(1:3) == 'E3T') THEN
          
          CELLS = CELLS +1
          
       ELSEIF(line(1:2) == 'ND') THEN
          
          NODES = NODES + 1
          
       END IF
    END DO
    REWIND(GRDINUNIT)


    IF(NODES < 3) CALL FATAL_ERROR("Invalid sms grid file?",&
         & "Less than three nodes in the file!")

    IF(CELLS < 1) CALL FATAL_ERROR("Invalid sms grid file?",&
         & "Less than one cell in the file!")


    WRITE(GRDOUTUNIT,*) "Node Number = ",NODES
    WRITE(GRDOUTUNIT,*) "Cell Number = ",CELLS

    IF(COR_ON) WRITE(COROUTUNIT,*) "Node Number = ",NODES

    IF(DEPTH_IN_GRID) WRITE(DPTHOUTUNIT,*) "Node Number = ",NODES

    IOS = 0
    DO WHILE(IOS == 0)
       line = ''
       READ(GRDINUNIT,'(a)',IOSTAT=IOS) line

       CALL SPLIT_STRING(LINE," ",lsplit)

       sizel = size(lsplit)

       IF(sizel == 0) THEN
          deallocate(lsplit)
          cycle
       END IF
          
       IF(lsplit(1) == 'E3T') THEN

          IF(sizel /= 6) Call fatal_error("Invalid line for CELL TYPE:",&
               & TRIM(line))

          ! READ A TRIANGLE
          !READ(lsplit(2),'(I)') ii1
          !READ(lsplit(3),'(I)') ii2
          !READ(lsplit(4),'(I)') ii3
          !READ(lsplit(5),'(I)') ii4
          READ(lsplit(2),*) ii1
          READ(lsplit(3),*) ii2
          READ(lsplit(4),*) ii3
          READ(lsplit(5),*) ii4


          WRITE(GRDOUTUNIT,'(I7,1X,I7,1X,I7,1X,I7)') ii1, ii2, ii3, ii4

       ELSEIF(lsplit(1) == 'ND') THEN
          ! READ NODE LOCATIONS

          IF(sizel /= 5) Call fatal_error("Invalid line for NODE TYPE:",&
               & TRIM(line))
          
          !READ(lsplit(2),'(I)') ii1
          !READ(lsplit(3),'(F)') ri2
          !READ(lsplit(4),'(F)') ri3
          !READ(lsplit(5),'(F)') ri4
          READ(lsplit(2),*) ii1
          READ(lsplit(3),*) ri2
          READ(lsplit(4),*) ri3
          READ(lsplit(5),*) ri4

          
          ro4 = ri4
          
          IF(GRID_INPUT_UNITS == GRID_OUTPUT_UNITS) THEN
             ro2 = ri2 
             ro3 = ri3
          ELSEIF(GRID_INPUT_UNITS == 'meters') THEN

             CALL METERS2DEGREES(ri2,ri3,PROJECTION_REFERENCE,ro2,ro3)

          ELSE
             
             CALL DEGREES2METERS(ri2,ri3,PROJECTION_REFERENCE,ro2,ro3)
          END IF

          WRITE(GRDOUTUNIT,'(I7,1X,ES19.12,1X,ES19.12,1x,ES19.12)') ii1, ro2, ro3, ro4


          ! DUMP THE CORIOLIS FILE
          IF(COR_ON) THEN
             
             IF(GRID_INPUT_UNITS == 'meters') THEN

                CALL METERS2DEGREES(ri2,ri3,PROJECTION_REFERENCE,cc1,cc2)
                WRITE(COROUTUNIT,'(ES19.12,1X,ES19.12,1x,ES19.12)') ro2, ro3, cc2
             ELSE

                WRITE(COROUTUNIT,'(ES19.12,1X,ES19.12,1x,ES19.12)') ro2, ro3, ri3

             END IF
                

          END IF
          
          IF (DEPTH_IN_GRID) THEN
             ! SET A MINIMUM DEPTH HERE!
             ro4 = max(ro4,DEPTH_MINIMUM)


             WRITE(DPTHOUTUNIT,'(ES19.12,1X,ES19.12,1x,ES19.12)') ro2, ro3, ro4
          END IF

   
       ELSEIF(lsplit(1) == 'NS') THEN

          IF(.not. NSOPEN) THEN

             nsopen = .true.
             NSCNT = NSCNT +1
             WRITE(CTEMP,'(I3.3)') nscnt
             NODESTRINGNAME= "SMSNODESTRING_"//trim(ctemp)
             
             pathnfile = trim(OUTPUT_DIR)//trim(NODESTRINGNAME)
             Call FOPEN(NODESTRINGUNIT,trim(pathnfile),'ofr')
             
             ! BURN ONE LINE
             write(NODESTRINGUNIT,*)
          END IF
          
          DO I = 2,sizel

             nns = nns +1
             !READ(lsplit(I),'(I)') ii1
             READ(lsplit(I),*) ii1
             
             IF (ii1 > 0) THEN
           
                WRITE(NODESTRINGUNIT,'(I7,1X,I7,1X,I7)') nns, ii1, 1
                
             ELSE

                ii1 = abs(ii1)
                WRITE(NODESTRINGUNIT,'(I7,1X,I7,1X,I7)') nns, ii1, 1
                nns = 0
                close(NODESTRINGUNIT)
                NSOPEN = .false.
             END IF
          END DO

       ELSEIF(lsplit(1) == 'BEGPARAMDEF') THEN
          ! END OF FILE
          exit

       END IF

       ! DEALLOCATE THE SPLIT STRING
       DEALLOCATE(LSPLIT)


    END DO
    
    CLOSE(GRDOUTUNIT)
    CLOSE(GRDINUNIT)
    IF(COR_ON) CLOSE(COROUTUNIT)


    ! CONVERT THE DEPTH FILE
    IF(DEPTH_ON .and. .not. DEPTH_IN_GRID) THEN
       
       WRITE(DPTHOUTUNIT,*) "Node Number = ",NODES
       
       I=0
       IOS=0
       DO WHILE(IOS == 0)
          
!          READ(DPTHINUNIT,*,IOSTAT=IOS,END=99) ri1, ri2, ri3

          line = ''
          READ(DPTHINUNIT,'(a)',IOSTAT=IOS) line 
          CALL SPLIT_STRING(LINE," ",lsplit)

          IF(size(lsplit)/=3) THEN
             deallocate(lsplit)
             CYCLE
          END IF

          !READ(lsplit(1),'(F)',IOSTAT=MYOS) ri1
          READ(lsplit(1),*,IOSTAT=MYOS) ri1
          IF(MYOS/=0) THEN
             deallocate(lsplit)
             CYCLE
          END IF

          !READ(lsplit(2),'(F)',IOSTAT=MYOS) ri2
          READ(lsplit(2),*,IOSTAT=MYOS) ri2
          IF(MYOS/=0) THEN
             deallocate(lsplit)
             CYCLE
          END IF

          !READ(lsplit(3),'(F)',IOSTAT=MYOS) ri3
          READ(lsplit(3),*,IOSTAT=MYOS) ri3
          IF(MYOS/=0) THEN
             deallocate(lsplit)
             CYCLE
          END IF
          
          I = I +1

          ! SET A MINIMUM DEPTH HERE!
          ro3 = max(ri3,DEPTH_MINIMUM)


          IF(PROJ_ON) THEN
             IF(DEPTH_INPUT_UNITS == GRID_OUTPUT_UNITS) THEN
                ro1 = ri1 
                ro2 = ri2
                
             ELSEIF(DEPTH_INPUT_UNITS == 'meters') THEN
                
                CALL METERS2DEGREES(ri1,ri2,PROJECTION_REFERENCE,ro1,ro2)
                
             ELSE
                
                CALL DEGREES2METERS(ri1,ri2,PROJECTION_REFERENCE,ro1,ro2)
             END IF
          ELSE
             
             ro1 = ri1
             ro2 = ri2
             
          END IF
          
          
          WRITE(DPTHOUTUNIT,'(ES19.12,1X,ES19.12,1x,ES19.12)') ro1, ro2, ro3
          
          
          
       END DO
          
       IF(NODES > I) THEN
          WRITE(IPT,*) "FOUND ",I,"; valid lines in depth file..."
          Call FATAL_ERROR('Unexpected end of Depth File: too few lines?',&
               & TRIM(DEPTH_FILE_IN))

       ELSEIF(NODES < I) THEN
          WRITE(IPT,*) "FOUND ",I,"; valid lines in depth file..."
          Call FATAL_ERROR('Unexpected end of Depth File: too many lines?',&
               & TRIM(DEPTH_FILE_IN))
       END IF


    END IF


  END SUBROUTINE CONVERT_SMS2DAT


end module mod_sms2grid
