MODULE MOD_nctest
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_INPUT
  USE MOD_TIME
  USE CONTROL
  USE LIMS
  USE ALL_VARS
  IMPLICIT NONE


CONTAINS


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
       if(MSR) WRITE(IPT,*) "You must specify: '--filename=<namelist>' "
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

          else if (opt_sng == "dbg_par" .or.opt_sng == "Dbg_Par"&
               & .or.opt_sng == "DBG_PAR") then

             dbg_par = .true.

          else if (opt_sng == "GRIDFILE" .or.opt_sng == "GridFile"&
               & .or.opt_sng == "gridfile") then

             call ftn_arg_get(arg_idx,arg_val,GRID_FILE) ! [sng] Input file
             GRID_FILE=GRID_FILE(1:ftn_strlen(GRID_FILE))
             ! Convert back to a fortran string!

          else if (opt_sng == "USE_MPI_IO" .or.opt_sng == "Use_Mpi_Io"&
               & .or.opt_sng == "use_mpi_io") then

             call ftn_arg_get(arg_idx,arg_val,USE_MPI_IO_MODE) ! [sng] Input file
!             USE_MPI_IO_MODE=USE_MPI_IO_MODE(1:ftn_strlen(USE_MPI_IO_MODE))
             ! Convert back to a fortran string!
             
          else if (opt_sng == "help" .or.opt_sng == "HELP" .or. opt_sng&
               & == "Help") then

             if(MSR) call MYHelpTxt


             call PSHUTDOWN
!!$   THIS DOES NOT SEEM PRACTICAL - MODIFY THE RUN FILE INSTEAD
!!$          else if (opt_sng == "CrashRestart") then
!!$             call ftn_arg_get(arg_idx,arg_val,CrashRestart) ! [lgc] Logical

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

          if(MSR) call MYHelpTxt
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

    WRITE(IPT,*) "! WELCOME TO THE NCTOOLS TEST SUIT"
    WRITE(IPT,*) "! OPTIONS:"
    WRITE(IPT,*) "! --GRIDFILE=<a valid FVCOM GRID FILE NAME>"
    WRITE(IPT,*) "! NOTES: This program is parallel!"
  END SUBROUTINE MYHELPTXT


  SUBROUTINE GET_FVCOM_GRID
    USE MOD_SETUP
    IMPLICIT NONE
    CHARACTER(LEN=80) FNAME
    INTEGER STATUS

    ! OPEN AND READ THE FVCOM GRID FILE
    IF (MSR) THEN
       FNAME = GRID_FILE
       WRITE(IPT,*) "OPENING GRIDFILE: "//TRIM(FNAME)
       Call FOPEN(GRIDUNIT,TRIM(FNAME),'cfr')
    END IF
    
    CALL LOAD_COLDSTART_GRID(NVG)
    KB = 5

    CALL SETUP_DOMAIN


  END SUBROUTINE GET_FVCOM_GRID





  SUBROUTINE TEST_NCLL
    IMPLICIT NONE

    CHARACTER(LEN=80) :: FNAME
    CHARACTER(LEN=6) :: attnum, dimnum,varnum,filenum

    integer:: i,j,k
    integer, Parameter:: test_size=987
    logical :: FOUND

    ! NC POINTERS
    TYPE(NCDIM), POINTER :: DIM

    TYPE(NCATT), POINTER :: ATT
    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCFILE), POINTER :: NCF

    nullify(DIM,ATT,VAR,NCF)


 !================================ ATTRIBUTES ===============================================

    IF (DBG_SET(DBG_LOG))    write(ipt,*) "Testing Att allocate and delete for Files"
     NCF => NEW_FILE()
     DO i = 1,test_size

        write(attnum,'(I6.6)') I

        ATT => NC_make_att(attnum,I)

        NCF => ADD(NCF,ATT)

     END DO
     IF (DBG_SET(DBG_LOG))   write(ipt,*) "Created N Atts: ",COUNT_ATT_LIST(NCF)

     k = test_size/2
     write(attnum,'(I6.6)') K
     ATT => FIND_ATT(NCF,attnum,FOUND)
     if (.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND ATT NAMED:"&
          &//trim(attnum))
     call print_att(att)


     IF (DBG_SET(DBG_LOG))  write(ipt,*) "deleting every other att"     
     DO i = 1,test_size,2

        IF ( I < test_size/2) THEN
           write(attnum,'(I6.6)') I
           CALL delete_att_link(NCF,attnum,FOUND)
           IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE ATT TO DE&
                &LETE?","Att name: "//trim(attnum))
        ELSE
           write(attnum,'(I6.6)') I
           CALL delete_att_link(NCF,I,FOUND)
           IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE ATT TO DE&
                &LETE?","AttID: "//trim(attnum))
        END IF

     END DO
     IF (DBG_SET(DBG_LOG))   write(ipt,*) "N Atts left: ",COUNT_ATT_LIST(NCF),"; Delete the list"
     call delete_att_list(NCF)
     IF (DBG_SET(DBG_LOG))   write(ipt,*) "N Atts left: ",COUNT_ATT_LIST(NCF)
     call KILL_FILE(NCF)


     IF (DBG_SET(DBG_LOG))   write(ipt,*) "Testing Att allocate and delete for Variables" 
    VAR => NEW_VAR()
    DO i = 1,test_size
       write(attnum,'(I6.6)') I
       ATT => NC_make_att(attnum,I)
       VAR => ADD(VAR,ATT)
    END DO
    IF (DBG_SET(DBG_LOG))   write(ipt,*) "Created N Atts: ",COUNT_ATT_LIST(VAR)

    k = test_size/2
    write(attnum,'(I6.6)') K
    ATT => FIND_ATT(VAR,attnum,FOUND)
    if (.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND ATT NAMED:"&
         &//trim(attnum))
    call print_att(att)

    IF (DBG_SET(DBG_LOG)) write(ipt,*) "deleting every other att" 
    DO i = 1,test_size,2
       
       IF ( I < test_size/2) THEN
          write(attnum,'(I6.6)') I
          CALL delete_att_link(VAR,attnum,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE ATT TO DE&
               &LETE?","Att name: "//trim(attnum))
       ELSE
          write(attnum,'(I6.6)') I
          CALL delete_att_link(VAR,I,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE ATT TO DE&
               &LETE?","AttID: "//trim(attnum))
       END IF

    END DO
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "N Atts left: ",COUNT_ATT_LIST(VAR),"; Delete the list"
    call delete_att_list(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "N Atts left: ",COUNT_ATT_LIST(VAR)
    CALL PRINT_VAR(VAR)

    call KILL_VAR(VAR)

    IF (DBG_SET(DBG_LOG)) write(ipt,*) "FINISHED Testing Att allocate and delete:" 

!================================ DIMENSIONS ===============================================

    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Testing DIM allocate and delete in FILES:" 
    NCF => NEW_FILE()
    DO i = 1,test_size
       write(dimnum,'(I6.6)') I
       DIM => NC_make_dim(dimnum,I)
       NCF => ADD(NCF,DIM)
    END DO
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Created N Dims: ",COUNT_DIM_LIST(NCF)

    k = test_size/2
    write(dimnum,'(I6.6)') K
    DIM => FIND_DIM(NCF,dimnum,FOUND)
    if (.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIM NAMED:"&
         &//trim(DIMnum))
    call print_DIM(DIM)

    IF (DBG_SET(DBG_LOG)) write(ipt,*) "deleting every other DIM"     
    DO i = 1,test_size,2
   
       IF ( I < test_size/2) THEN
          write(dimnum,'(I6.6)') I
          CALL delete_dim_link(NCF,dimnum,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE DIM TO DE&
               &LETE?","DIM name: "//trim(dimnum))
       ELSE
          write(dimnum,'(I6.6)') I
          CALL delete_DIM_link(NCF,I,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE DIM TO DE&
               &LETE?","DIMID: "//trim(dimnum))
       END IF

    END DO
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "N DIM left: ",COUNT_DIM_LIST(NCF),"; Delete the list"
    call delete_DIM_list(NCF)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "N DIM left: ",COUNT_DIM_LIST(NCF)
    call KILL_FILE(NCF)
    

    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Testing DIM allocate and delete in VARIABLES:" 
    NCF => NEW_FILE()
    VAR => NEW_VAR()
    VAR%VARNAME="test"
    DO i = 1,test_size
       write(dimnum,'(I6.6)') I
       DIM => NC_make_dim(dimnum,I)
       VAR => ADD(VAR,DIM)
    END DO
    
    ! TO SET THE DIMID OF EACH DIMENSION WE MUST ADD THE VAR TO A FILE
    NCF => ADD(NCF,VAR)

    VAR => FIND_VAR(NCF,"test",FOUND)
    if (.not. FOUND) CALL FATAL_ERROR("could not find the test variable?")

    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Dims IN VAR: ",COUNT_DIM_LIST(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Dims IN FILE: ",count_dim_list(ncf)

    k = test_size/2
    write(dimnum,'(I6.6)') K
    DIM => FIND_DIM(VAR,dimnum,FOUND)
    if (.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIM NAMED:"&
         &//trim(DIMnum))
    call print_DIM(DIM)

    IF (DBG_SET(DBG_LOG)) write(ipt,*) "deleting every other DIM"     
    DO i = 1,test_size,2
   
       IF ( I < test_size/4) THEN
          write(dimnum,'(I6.6)') I
          CALL delete_dim_link(VAR,dimnum,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE DIM TO DE&
               &LETE?","DIM name: "//trim(dimnum))
       ELSE IF ( I < test_size/2) THEN
          write(dimnum,'(I6.6)') I
          CALL delete_DIM_link(VAR,I,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE DIM TO DE&
               &LETE?","DIMID: "//trim(dimnum))

       ELSE 
          write(dimnum,'(I6.6)') I
          CALL delete_DIM_link(NCF,I,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE DIM TO DE&
               &LETE?","DIMID: "//trim(dimnum))
       END IF

    END DO
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Dims IN VAR: ",COUNT_DIM_LIST(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Dims IN FILE: ",count_dim_list(ncf)

    IF (DBG_SET(DBG_LOG))  write(ipt,*) "N DIM left: ",COUNT_DIM_LIST(VAR),"; Delete the list"
    call delete_DIM_list(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "N DIM left: ",COUNT_DIM_LIST(VAR)

    ! DO not delete the var, it belongs to the file!
    call print_file(NCF)
    call kill_FILE(NCF)
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "FINISHED Testing DIM allocate and delete in VARIABLES:" 

!================================ VARIABLES ===============================================
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Testing VAR allocate and delete in FILES:" 
    NCF => NEW_FILE()
    DO i = 1,test_size
       write(varnum,'(I6.6)') I
       VAR => NC_make_avar(varnum,I) ! a scalar integer variable
       NCF => ADD(NCF,VAR)
    END DO
    
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "VARS IN FILE: ",count_var_list(ncf)

    k = test_size/2
    write(varnum,'(I6.6)') K
    VAR => FIND_VAR(NCF,varnum,FOUND)
    if (.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VAR NAMED:"&
         &//trim(VARnum))
    call print_VAR(VAR)

    IF (DBG_SET(DBG_LOG))  write(ipt,*) "deleting every other VAR"     
    DO i = 1,test_size,2
   
       IF ( I < test_size/2) THEN
          write(varnum,'(I6.6)') I
          CALL delete_var_link(NCF,varnum,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE VAR TO DE&
               &LETE?","VAR name: "//trim(dimnum))
       ELSE
          write(varnum,'(I6.6)') I
          CALL delete_var_link(NCF,I,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE VAR TO DE&
               &LETE?","VARID: "//trim(varnum))
       END IF

    END DO
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "VARS IN FILE: ",count_var_list(ncf),"; Delete the list"

    call delete_VAR_list(NCF)
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "N VAR left: ",COUNT_VAR_LIST(NCF)

    ! DO not delete the var, it belongs to the file!
    call print_file(NCF)
    call kill_FILE(NCF)

!================================ FILES ===============================================
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Testing FILE allocate and delete in FILEHEAD:" 
    
    DO i = 1,test_size
       write(FILEnum,'(I6.6)') I
       CALL NC_INIT(NCF,FILENUM)
       FILEHEAD => ADD(FILEHEAD,NCF)
    END DO
    
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "VARS IN FILE: ",count_FILE_list(FILEHEAD)

    k = test_size/2
    write(FILEnum,'(I6.6)') K
    NCF => FIND_FILE(FILEHEAD,FILEnum,FOUND)
    if (.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND FILE NAMED:"&
         &//trim(FILEnum))
    call print_FILE(NCF)

    IF (DBG_SET(DBG_LOG))  write(ipt,*) "deleting every other FILE"     
    DO i = 1,test_size,2
   
       IF ( I < test_size/2) THEN
          write(FILEnum,'(I6.6)') I
          CALL delete_FILE_link(FILEHEAD,FILEnum,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE VAR TO DE&
               &LETE?","VAR name: "//trim(filenum))
       ELSE
          ! NONE OF THE FILES ARE OPEN, JUST DELETE THEM AS -1
          CALL delete_FILE_link(FILEHEAD,-1,FOUND)
          IF (.not. FOUND) Call fatal_error("COULD NOT FIND THE VAR TO DE&
               &LETE?","VARID: -1 (FIRST NOT OPEN FILE!)")
       END IF

    END DO
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "VARS IN FILE: ",count_file_list(filehead),"; Delete the list"

    call delete_FILE_list(FILEHEAD)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "N VAR left: ",COUNT_FILE_LIST(FILEHEAD)

    ! DO not delete the var, it belongs to the file!
    call print_file_LIST(FILEHEAD)

    CALL KILL_FILEHEAD(FILEHEAD)

!================================ FILES ===============================================
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "Testing COPY for DIMS, ATTS, VARS and FILES" 
  
    NCF => NEW_FILE()
    VAR => NEW_VAR()

    VAR%VARNAME="test"
    DO i = 1,3
       write(dimnum,'(I6.6)') I
       DIM => NC_make_dim(dimnum,I)
       VAR => ADD(VAR,DIM)

       ATT => NC_MAKE_ATT(dimnum,2*I)
       VAR => ADD(VAR,ATT)

    END DO
    ! DIM LIST
    call print_dim_list(VAR)
!    call print_dim_list(COPY_VAR(VAR))

    

     
    NCF => ADD(NCF,COPY_VAR(VAR))
    VAR%VARNAME="test2"

    NCF => ADD(NCF,COPY_VAR(VAR))
    VAR%VARNAME="test3"

    ! THIS TIME ADD THE MEMORY AND CLEAR
    NCF => ADD(NCF,VAR)

    call print_var_list(NCF)
    call print_dim_list(NCF)

    VAR => find_var(NCF,"test2",FOUND)
    call print_att_list(VAR)
    call print_dim_list(VAR)

    
    FILEHEAD => NEW_FILEHEAD()

    DO i = 1,3
       write(filenum,'(I6.6)') I
       NCF%FNAME=filenum
       FILEHEAD => ADD(FILEHEAD,COPY_FILE(NCF))
    END DO

    CALL PRINT_FILE_LIST(FILEHEAD)

    call kill_file(NCF)

    CALL KILL_FILEHEAD(FILEHEAD)


  END SUBROUTINE TEST_NCLL


  SUBROUTINE TEST_BASIC_IO
    IMPLICIT NONE

    CHARACTER(LEN=80) :: FNAME
    CHARACTER(LEN=4) :: procnum

    integer:: i,j,k
    logical :: FOUND

    ! NC POINTERS
    TYPE(NCDIM), POINTER :: DIM_1,DIM_1r
    TYPE(NCDIM), POINTER :: DIM_2,DIM_2r
    TYPE(NCDIM), POINTER :: DIM_3,DIM_3r
    TYPE(NCDIM), POINTER :: DIM_4,DIM_4r
    TYPE(NCDIM), POINTER :: DIM_5,DIM_5r

    TYPE(NCATT), POINTER :: ATT
    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCFILE), POINTER :: NCFR

    INTEGER, allocatable :: strt(:), cnt(:), strd(:)

    ! SOME DUMMY VARIABLES
    INTEGER :: int, intr
    REAL(SP) :: flt, fltr
    REAL(DP) :: dbl, dblr

    CHARACTER(LEN=80) :: str20,str10, str20r, str10r
    CHARACTER(LEN=80), ALLOCATABLE :: str_vec10(:),str_vec20(:), str_vec10r(:),str_vec20r(:)
    
    INTEGER, ALLOCATABLE :: int1(:),int1r(:)
    REAL(SP), ALLOCATABLE :: flt1(:),flt1r(:)
    REAL(DP), ALLOCATABLE :: dbl1(:),dbl1r(:)

    INTEGER, ALLOCATABLE :: int2(:,:),int2r(:,:)
    REAL(SP), ALLOCATABLE :: flt2(:,:),flt2r(:,:)
    REAL(DP), ALLOCATABLE :: dbl2(:,:),dbl2r(:,:)

    INTEGER, ALLOCATABLE :: int3(:,:,:),int3r(:,:,:)
    REAL(SP), ALLOCATABLE :: flt3(:,:,:),flt3r(:,:,:)
    REAL(DP), ALLOCATABLE :: dbl3(:,:,:),dbl3r(:,:,:)
    
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "TEST writing a file from the local processor"

    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    write(procnum,'(I4.4)')myid
    FNAME= "./results/FTEST1_"//trim(procnum)//".nc"
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "WRITING FILE: "//TRIM(FNAME)

    NCF%FNAME=FNAME

    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='title',values="A TEST FILE") 
    NCF => ADD(NCF,ATT)

    ! CAN SET A VALUE NOT KNOW BY IOPROC
    ATT => NC_MAKE_runtime_ATT_CHR(name='run_time_att',values="A TEST FILE") 
    NCF => ADD(NCF,ATT)

    ! NOTE SHOULD ADD INTERFACE TO MAKE RUNTIME ATTS OF ALL TYPES...
    ! ON THE TODO LIST

    ATT => NC_MAKE_ATT(name='Float',values=5.0) 
    NCF => ADD(NCF,ATT)

    allocate(flt1(2)); flt1=(/5.0_SP, 3.14159265358979323846_SP/)
    ATT => NC_MAKE_ATT(name='Floats',values=flt1) 
    NCF => ADD(NCF,ATT)
    deallocate(flt1)

    ATT => NC_MAKE_ATT(name='int',values=5) 
    NCF => ADD(NCF,ATT)
    
    allocate(int1(2)); int1=(/5, 100/)
    ATT => NC_MAKE_ATT(name='ints',values=int1) 
    NCF => ADD(NCF,ATT)
    deallocate(int1)

    ATT => NC_MAKE_ATT(name='Double',values=5.0_DP) 
    NCF => ADD(NCF,ATT)

    allocate(dbl1(2)); dbl1=(/5.0_DP, 3.14159265358979323846_DP/)
    ATT => NC_MAKE_ATT(name='Doubles',values=dbl1)
    NCF => ADD(NCF,ATT)
    deallocate(dbl1)


    DIM_1       => NC_MAKE_DIM(name='dim1',len=2)

    DIM_2       => NC_MAKE_DIM(name='dim2',len=3)

    DIM_3       => NC_MAKE_DIM(name='dim3',len=2)

    DIM_4       => NC_MAKE_DIM(name='dim4',len=NF90_UNLIMITED)

    ! CAN SET A VALUE NOT KNOW BY IOPROC
    DIM_5       => NC_MAKE_RUNTIME_DIM(name='dim5',len=20)

    ! SCALAR VARIABLES
    VAR  => NC_MAKE_AVAR(name='flt',values=3.14159265358979323846_SP)
    ATT  => NC_MAKE_ATT(name='variable_att',values='something')
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='dbl',values=3.14159265358979323846_DP)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='int',values=123456)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED SCALARS
    VAR  => NC_MAKE_AVAR(name='u_flt',values&
         &=3.14159265358979323846_SP,DIM1=DIM_4)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_dbl',values&
         &=3.14159265358979323846_DP,DIM1=DIM_4)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_int',values=int,DIM1=DIM_4)
    NCF  => ADD(NCF,VAR)

    ! VECTORS
    allocate(flt1(2)); flt1=2
    VAR  => NC_MAKE_AVAR(name='flt1',values=flt1,DIM1=DIM_1)
    NCF  => ADD(NCF,VAR)

    allocate(dbl1(3)); dbl1=(/0.0_DP,5.0_DP,3.14159265358979323846_DP/)
    VAR  => NC_MAKE_AVAR(name='dbl1',values=dbl1,DIM1=DIM_2)
    NCF  => ADD(NCF,VAR)
    
    allocate(int1(2)); int1=-100
    VAR  => NC_MAKE_AVAR(name='int1',values=int1,DIM1=DIM_1)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED VECTORS
    VAR  => NC_MAKE_AVAR(name='u_flt1',values=flt1,DIM1=DIM_1,DIM2=DIM_4)
    NCF  => ADD(NCF,VAR)
    
    VAR  => NC_MAKE_AVAR(name='u_dbl1',values=dbl1,DIM1=DIM_2,DIM2=DIM_4)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_int1',values=int1,DIM1=DIM_1,DIM2=DIM_4)
    NCF  => ADD(NCF,VAR)

    ! ARRAYS
    allocate(flt2(2,3)); flt2=1.5
    VAR  => NC_MAKE_AVAR(name='flt2',values=flt2,DIM1=DIM_1,DIM2=DIM_2)
    NCF  => ADD(NCF,VAR)

    ! TEST USING DIMENSION TWICE IN THE SAME VARIABLE!
    allocate(dbl2(3,3)); dbl2=3.2_DP
    VAR  => NC_MAKE_AVAR(name='dbl2',values=dbl2,DIM1=DIM_2,DIM2=DIM_2)
    NCF  => ADD(NCF,VAR)
    
    allocate(int2(2,2)); int2=-100
    VAR  => NC_MAKE_AVAR(name='int2',values=int2,DIM1=DIM_1,DIM2=DIM_3)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED ARRAYS
    VAR  => NC_MAKE_AVAR(name='u_flt2',values=flt2,DIM1=DIM_1,DIM2=DIM_2,DIM3=DIM_4)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_dbl2',values=dbl2,DIM1=DIM_2,DIM2=DIM_2,DIM3=DIM_4)
    NCF  => ADD(NCF,VAR)
    
    VAR  => NC_MAKE_AVAR(name='u_int2',values=int2,DIM1=DIM_1,DIM2=DIM_3,DIM3=DIM_4)
    NCF  => ADD(NCF,VAR)


    ! CUBES
    allocate(flt3(2,3,2)); flt3=-1.5
    VAR  => NC_MAKE_AVAR(name='flt3',values=flt3,DIM1=DIM_1,DIM2=DIM_2,DIM3=DIM_3)
    NCF  => ADD(NCF,VAR)

    allocate(dbl3(3,3,2)); dbl3=0.0_DP
    VAR  => NC_MAKE_AVAR(name='dbl3',values=dbl3,DIM1=DIM_2,DIM2=DIM_2,DIM3=DIM_3)
    NCF  => ADD(NCF,VAR)
    
    allocate(int3(3,2,2))
    DO i =1,3
       DO j = 1,2
          DO k = 1,2
             int3(i,j,k) = i*j*k
          END DO
       END DO
    END DO
    VAR  => NC_MAKE_AVAR(name='int3',values=int3,DIM1=DIM_2,DIM2=DIM_3,DIM3=DIM_3)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED CUBES
    VAR  => NC_MAKE_AVAR(name='u_flt3',values=flt3,DIM1=DIM_1,DIM2=DIM_2,DIM3=DIM_3,DIM4=DIM_4)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_dbl3',values=dbl3,DIM1=DIM_2,DIM2=DIM_2,DIM3=DIM_3,DIM4=DIM_4)
    NCF  => ADD(NCF,VAR)
    
    VAR  => NC_MAKE_AVAR(name='u_int3',values=int3,DIM1=DIM_2,DIM2=DIM_3,DIM3=DIM_3,DIM4=DIM_4)
    NCF  => ADD(NCF,VAR)

    ! CHARCTER STRINGS
    str10="String ten!"
    VAR  => NC_MAKE_AVAR(name='short_string',values=str10,DIM1=DIM_5)
    NCF  => ADD(NCF,VAR)

    str20="This one is twenty!!"
    VAR  => NC_MAKE_AVAR(name='string',values=str20,DIM1=DIM_5)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED CHARACTER STRINGS
    VAR  => NC_MAKE_AVAR(name='u_short_string',values=str10,DIM1=DIM_5,DIM2=DIM_4)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_string',values=str20,DIM1=DIM_5,DIM2=DIM_4)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(str_vec10(3),str_vec20(3))
    DO i=1,3
       str_vec10(i)=TRIM(str10)
       str_vec20(i)=TRIM(str20)
    END DO
    ! VECTOR CHARACTER STRINGS
    VAR  => NC_MAKE_AVAR(name='vec_short_string',values=str_vec10,DIM1=DIM_5,DIM2=DIM_2)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='vec_string',values=str_vec20,DIM1=DIM_5,DIM2=DIM_2)
    NCF  => ADD(NCF,VAR)

    ! WRITE THIS ONLY FROM THE LOCAL PROC...
    call NC_WRITE_FILE(NCF,LOCAL_ONLY=.TRUE.)

    ! ADD FTIME
    NCF%FTIME => NEW_FTIME()

    ! TEST WRITING THE UNLIMITED DATA
    int=1
    NCF%FTIME%NEXT_STKCNT = 1
    call NC_WRITE_FILE(NCF,LOCAL_ONLY=.TRUE.)

    int=2
    int2=-200
    NCF%FTIME%NEXT_STKCNT = 2
    call NC_WRITE_FILE(NCF,LOCAL_ONLY=.TRUE.)


    ! TEST WRITING THE UNLIMITED DATA IN A WIERD ORDER
    int=5
    int2=-500
    NCF%FTIME%NEXT_STKCNT = 5
    call NC_WRITE_FILE(NCF,LOCAL_ONLY=.TRUE.)

    int=4
    int2=-400
    NCF%FTIME%NEXT_STKCNT = 4
    call NC_WRITE_FILE(NCF,LOCAL_ONLY=.TRUE.)

    int=3
    int2=-300
    NCF%FTIME%NEXT_STKCNT = 3
    call NC_WRITE_FILE(NCF,LOCAL_ONLY=.TRUE.)


    ! NOW TEST WRITING USING START,COUNT,STRIDE
!!$    call print_file(NCF)
!!$
!!$    VAR => find_var(NCF,"u_int3",FOUND)
!!$    IF(.not.FOUND) call fatal_error("Can't find u_int3?")
!!$    int3=5
!!$
!!$    CALL NC_OPEN(NCF)
!!$
!!$    call print_var(var)
!!$
!!$    allocate(strt(4),cnt(4),strd(4))
!!$    strt=1
!!$    cnt=(/3,1,2,2/)
!!$    strd=1
!!$    write(ipt,*)cnt
!!$
!!$    call NC_WRITE_VAR(VAR,.TRUE.,.FALSE.,MYID,IOSTART=strt,IOCOUNT=cnt)
!!$
!!$    write(ipt,*) Var%cub_int
!!$    
    

    ! NOW TEST READING THE DATA BACK
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "TEST READING DATA FROM A local processor"
    NCFR=>NEW_FILE()
    NCFR%FNAME=NCF%FNAME
    CALL NC_OPEN(NCFR)

    CALL NC_LOAD(NCFR)

    CALL PRINT_FILE(NCFR)
    
    call print_att_list(NCFR)

    CALL print_dim_list(NCFR)

    DIM_1r=>FIND_DIM(NCFR,"dim1",FOUND)
    CALL PRINT_DIM(DIM_1R)

    DIM_2r=>FIND_DIM(NCFR,"dim2",FOUND)
    CALL PRINT_DIM(DIM_2R)

    DIM_3r=>FIND_DIM(NCFR,"dim3",FOUND)
    CALL PRINT_DIM(DIM_3R)

    DIM_4r=>FIND_DIM(NCFR,"dim4",FOUND)
    CALL PRINT_DIM(DIM_4R)

    DIM_5r=>FIND_DIM(NCFR,"dim5",FOUND)
    CALL PRINT_DIM(DIM_5R)


    ! SCALARS
    VAR => find_var(NCFR,"flt",found)
    CALL NC_CONNECT_AVAR(VAR,fltr)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "flt=",fltr

    call print_Att_list(VAR)
    call Print_dim_list(VAR)

    VAR => find_var(NCFR,"int",found)
    CALL NC_CONNECT_AVAR(VAR,intr)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "int=",intr

    VAR => find_var(NCFR,"dbl",found)
    CALL NC_CONNECT_AVAR(VAR,dblr)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "dbl=",dblr

    ! UNLIMITED SCLS
    VAR => find_var(NCFR,"u_flt",found)
    CALL NC_CONNECT_AVAR(VAR,fltr)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_flt=",fltr

    VAR => find_var(NCFR,"u_int",found)
    CALL NC_CONNECT_AVAR(VAR,intr)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_int=",intr

    VAR => find_var(NCFR,"u_dbl",found)
    CALL NC_CONNECT_AVAR(VAR,dblr)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_dbl=",dblr

    ! VECTORS
    allocate(flt1r(2))
    allocate(dbl1r(3))
    allocate(int1r(2))
    VAR => find_var(NCFR,"flt1",found)
    CALL NC_CONNECT_AVAR(VAR,flt1r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "flt1=",flt1r

    call print_Att_list(VAR)
    call Print_dim_list(VAR)

    VAR => find_var(NCFR,"int1",found)
    CALL NC_CONNECT_AVAR(VAR,int1r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "int1=",int1r

    VAR => find_var(NCFR,"dbl1",found)
    CALL NC_CONNECT_AVAR(VAR,dbl1r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "dbl1=",dbl1r

    ! UNLIMITED VECTORS
    VAR => find_var(NCFR,"u_flt1",found)
    CALL NC_CONNECT_AVAR(VAR,flt1r)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_flt1=",flt1r

    VAR => find_var(NCFR,"u_int1",found)
    CALL NC_CONNECT_AVAR(VAR,int1r)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_int1=",int1r

    VAR => find_var(NCFR,"u_dbl1",found)
    CALL NC_CONNECT_AVAR(VAR,dbl1r)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_dbl1=",dbl1r

    ! ARRAYS
    allocate(flt2r(2,3))
    allocate(dbl2r(3,3))
    allocate(int2r(2,2))
    VAR => find_var(NCFR,"flt2",found)
    CALL NC_CONNECT_AVAR(VAR,flt2r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "flt2=",flt2r

    VAR => find_var(NCFR,"int2",found)
    CALL NC_CONNECT_AVAR(VAR,int2r)
    call print_var(VAR)
    CALL PRINT_DIM_LIST(VAR)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "int2=",int2r

    VAR => find_var(NCFR,"dbl2",found)
    CALL NC_CONNECT_AVAR(VAR,dbl2r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "dbl2=",dbl2r

    ! UNLIMITED ARRAYS
    VAR => find_var(NCFR,"u_flt2",found)
    CALL NC_CONNECT_AVAR(VAR,flt2r)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_flt2=",flt2r

    VAR => find_var(NCFR,"u_int2",found)
    CALL NC_CONNECT_AVAR(VAR,int2r)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_int2=",int2r

    VAR => find_var(NCFR,"u_dbl2",found)
    CALL NC_CONNECT_AVAR(VAR,dbl2r)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_dbl2=",dbl2r


    ! CUBES
    allocate(flt3r(2,3,2))
    allocate(dbl3r(3,3,2))
    allocate(int3r(3,2,2))
    VAR => find_var(NCFR,"flt3",found)
    CALL NC_CONNECT_AVAR(VAR,flt3r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "flt3=",flt3r

    VAR => find_var(NCFR,"int3",found)
    CALL NC_CONNECT_AVAR(VAR,int3r)
    call print_var(VAR)
    CALL PRINT_DIM_LIST(VAR)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "int3=",int3r

    VAR => find_var(NCFR,"dbl3",found)
    CALL NC_CONNECT_AVAR(VAR,dbl3r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "dbl3=",dbl3r

    ! UNLIMITED CUBES
    VAR => find_var(NCFR,"u_flt3",found)
    CALL NC_CONNECT_AVAR(VAR,flt3r)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_flt3=",flt3r

    VAR => find_var(NCFR,"u_int3",found)
    CALL NC_CONNECT_AVAR(VAR,int3r)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_int3=",int3r

    VAR => find_var(NCFR,"u_dbl3",found)
    CALL NC_CONNECT_AVAR(VAR,dbl3r)
    call nc_read_var(VAR,STKCNT=2,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_dbl3=",dbl3r


    ! CHARCTER STRINGS
    VAR => find_var(NCFR,"short_string",found)
    CALL NC_CONNECT_AVAR(VAR,str10r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "str10=",trim(str10r)

    VAR => find_var(NCFR,"string",found)
    CALL NC_CONNECT_AVAR(VAR,str20r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "str20=",trim(str20r)

    ALLOCATE(str_vec10r(3),str_vec20r(3))
    VAR => find_var(NCFR,"vec_short_string",found)
    CALL NC_CONNECT_AVAR(VAR,str_vec10r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    DO I = 1,3
       IF (DBG_SET(DBG_LOG)) write(ipt,*) "str_vec10="&
            &,trim(str_vec10r(I))
    END DO

    VAR => find_var(NCFR,"vec_string",found)
    CALL NC_CONNECT_AVAR(VAR,str_vec20r)
    call nc_read_var(VAR,PARALLEL=.FALSE.)
    DO I = 1,3
       IF (DBG_SET(DBG_LOG)) write(ipt,*) "str_vec20=",trim(str_vec20r(I))
    END DO


  END SUBROUTINE TEST_BASIC_IO
    
  SUBROUTINE TEST_UNSTRUCTURED_IO
    IMPLICIT NONE
    CHARACTER(LEN=80) :: FNAME
    CHARACTER(LEN=4) :: procnum

    integer:: i,j,k,ierr
    logical :: FOUND

    ! NC POINTERS
    TYPE(NCDIM), POINTER :: DIM_1,DIM_1r
    TYPE(NCDIM), POINTER :: DIM_2,DIM_2r
    TYPE(NCDIM), POINTER :: DIM_3,DIM_3r
    TYPE(NCDIM), POINTER :: DIM_4,DIM_4r
    TYPE(NCDIM), POINTER :: DIM_5,DIM_5r
    TYPE(NCDIM), POINTER :: DIM_6,DIM_6r

    TYPE(NCATT), POINTER :: ATT
    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCFILE), POINTER :: NCFR

    INTEGER, allocatable :: strt(:), cnt(:), strd(:)

    ! SOME DUMMY VARIABLES
    INTEGER :: int, intr
    REAL(SP) :: flt, fltr
    REAL(DP) :: dbl, dblr

    CHARACTER(LEN=80) :: str20,str10, str20r, str10r
    CHARACTER(LEN=80), ALLOCATABLE :: str_vec10(:),str_vec20(:), str_vec10r(:),str_vec20r(:)
    
    INTEGER, ALLOCATABLE :: int1(:),int1r(:)
    REAL(SP), ALLOCATABLE :: flt1(:),flt1r(:)
    REAL(DP), ALLOCATABLE :: dbl1(:),dbl1r(:)

    INTEGER, ALLOCATABLE :: int2(:,:),int2r(:,:)
    REAL(SP), ALLOCATABLE :: flt2(:,:),flt2r(:,:)
    REAL(DP), ALLOCATABLE :: dbl2(:,:),dbl2r(:,:)

    INTEGER, ALLOCATABLE :: int3(:,:,:),int3r(:,:,:)
    REAL(SP), ALLOCATABLE :: flt3(:,:,:),flt3r(:,:,:)
    REAL(DP), ALLOCATABLE :: dbl3(:,:,:),dbl3r(:,:,:)
    
    IF(PAR .and. DBG_SET(DBG_LOG))  write(ipt,*) "TEST PARALLEL READ/WRITE"
    
    ! SEPERATE THE IO_PROC AND ALLOCATE MEMORY DIFFERENTLY!
    IF(IOPROCID == MYID) THEN

       allocate(flt1(mgl));
       allocate(dbl1(ngl));
       allocate(int1(ngl));
       
       allocate(flt2(mgl,3));
       allocate(dbl2(ngl,3));
       allocate(int2(ngl,4)); 
       
       allocate(flt3(ngl,3,3));
       allocate(dbl3(ngl,3,3));
       allocate(int3(mgl,3,4))

    ELSE
       
       flt = 3.14159265358979323846_SP
       dbl = 3.14159265358979323846_DP
       int = nprocs
       
       allocate(flt1(mt)); flt1=2.0_SP
       allocate(dbl1(n)); dbl1=1.0_DP
       allocate(int1(0:N)); int1=MYID
       
       allocate(flt2(0:MT,3)); flt2=1.5
       allocate(dbl2(0:N,3)); dbl2=3.2_DP
       allocate(int2(1:NT,4)); int2=-100
       
       allocate(flt3(0:NT,3,3)); flt3=-1.5
       allocate(dbl3(1:N,3,3)); dbl3=3.0_DP
       allocate(int3(0:MT,3,4))
       DO i =1,MT
          DO j = 1,3
             DO k = 1,4
                int3(i,j,k) = NGID_X(i)*j*k
             END DO
          END DO
       END DO
       
       
    END IF
    



    NCF => NEW_FILE()

    FNAME= "./results/FTEST2.nc"
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "WRITING FILE: "//TRIM(FNAME)
    NCF%FNAME=FNAME

    ATT => NC_MAKE_ATT(name='title',values="A TEST FILE") 
    NCF => ADD(NCF,ATT)


    DIM_1       => NC_MAKE_DIM(name='dim1',len=MGL)

    DIM_2       => NC_MAKE_DIM(name='dim2',len=NGL)

    DIM_3       => NC_MAKE_DIM(name='dim3',len=3)

    DIM_4       => NC_MAKE_DIM(name='dim4',len=4)

    DIM_5       => NC_MAKE_DIM(name='dim5',len=NF90_UNLIMITED)

    DIM_6       => NC_MAKE_DIM(name='dim6',len=20)

    ! SCALAR VARIABLES
    VAR  => NC_MAKE_AVAR(name='flt',values=flt)
    ATT  => NC_MAKE_ATT(name='variable_att',values='something')
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='dbl',values=dbl)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='int',values=int)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED SCALARS
    VAR  => NC_MAKE_AVAR(name='u_flt',values=flt,DIM1=DIM_5)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_dbl',values=dbl,DIM1=DIM_5)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_int',values=int,DIM1=DIM_5)
    NCF  => ADD(NCF,VAR)


    ! VECTORS
    VAR  => NC_MAKE_AVAR(name='flt1',values=flt1,DIM1=DIM_1)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='dbl1',values=dbl1,DIM1=DIM_2)
    NCF  => ADD(NCF,VAR)
    
    VAR  => NC_MAKE_AVAR(name='int1',values=int1,DIM1=DIM_2)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED VECTORS
    VAR  => NC_MAKE_AVAR(name='u_flt1',values=flt1,DIM1=DIM_1,DIM2=DIM_5)
    NCF  => ADD(NCF,VAR)
    
    VAR  => NC_MAKE_AVAR(name='u_dbl1',values=dbl1,DIM1=DIM_2,DIM2=DIM_5)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_int1',values=int1,DIM1=DIM_2,DIM2=DIM_5)
    NCF  => ADD(NCF,VAR)

    ! ARRAYS
    VAR  => NC_MAKE_AVAR(name='flt2',values=flt2,DIM1=DIM_1,DIM2=DIM_3)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='dbl2',values=dbl2,DIM1=DIM_2,DIM2=DIM_3)
    NCF  => ADD(NCF,VAR)
    
    VAR  => NC_MAKE_AVAR(name='int2',values=int2,DIM1=DIM_2,DIM2=DIM_4)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED ARRAYS
    VAR  => NC_MAKE_AVAR(name='u_flt2',values=flt2,DIM1=DIM_1,DIM2=DIM_3,DIM3=DIM_5)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_dbl2',values=dbl2,DIM1=DIM_2,DIM2=DIM_3,DIM3=DIM_5)
    NCF  => ADD(NCF,VAR)
    
    VAR  => NC_MAKE_AVAR(name='u_int2',values=int2,DIM1=DIM_2,DIM2=DIM_4,DIM3=DIM_5)
    NCF  => ADD(NCF,VAR)

    ! CUBES
    VAR  => NC_MAKE_AVAR(name='flt3',values=flt3,DIM1=DIM_2,DIM2=DIM_3,DIM3=DIM_3)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='dbl3',values=dbl3,DIM1=DIM_2,DIM2=DIM_3,DIM3=DIM_3)
    NCF  => ADD(NCF,VAR)
    
    VAR  => NC_MAKE_AVAR(name='int3',values=int3,DIM1=DIM_1,DIM2=DIM_3,DIM3=DIM_4)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED CUBES
    VAR  => NC_MAKE_AVAR(name='u_flt3',values=flt3,DIM1=DIM_2,DIM2=DIM_3,DIM3=DIM_3,DIM4=DIM_5)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_dbl3',values=dbl3,DIM1=DIM_2,DIM2=DIM_3,DIM3=DIM_3,DIM4=DIM_5)
    NCF  => ADD(NCF,VAR)
    
    VAR  => NC_MAKE_AVAR(name='u_int3',values=int3,DIM1=DIM_1,DIM2=DIM_3,DIM3=DIM_4,DIM4=DIM_5)
    NCF  => ADD(NCF,VAR)

    ! CHARCTER STRINGS
    str10="String ten!"
    VAR  => NC_MAKE_AVAR(name='short_string',values=str10,DIM1=DIM_6)
    NCF  => ADD(NCF,VAR)

    str20="This one is twenty!!"
    VAR  => NC_MAKE_AVAR(name='string',values=str20,DIM1=DIM_6)
    NCF  => ADD(NCF,VAR)

    ! UNLIMITED CHARACTER STRINGS
    VAR  => NC_MAKE_AVAR(name='u_short_string',values=str10,DIM1=DIM_6,DIM2=DIM_5)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='u_string',values=str20,DIM1=DIM_6,DIM2=DIM_5)
    NCF  => ADD(NCF,VAR)

    ALLOCATE(str_vec10(3),str_vec20(3))
    DO i=1,3
       str_vec10(i)=TRIM(str10)
       str_vec20(i)=TRIM(str20)
    END DO
    ! VECTOR CHARACTER STRINGS
    VAR  => NC_MAKE_AVAR(name='vec_short_string',values=str_vec10,DIM1=DIM_6,DIM2=DIM_3)
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='vec_string',values=str_vec20,DIM1=DIM_6,DIM2=DIM_3)
    NCF  => ADD(NCF,VAR)

    ! WRITE THIS ONLY FROM THE LOCAL PROC...
    call NC_WRITE_FILE(NCF)
    ! ADD FTIME
    NCF%FTIME => NEW_FTIME()

    ! TEST WRITING THE UNLIMITED DATA
    int=1
    NCF%FTIME%NEXT_STKCNT = 1
    call NC_WRITE_FILE(NCF)

    int=2
    NCF%FTIME%NEXT_STKCNT = 2
    call NC_WRITE_FILE(NCF)

# if defined(MULTIPROCESSOR)
    call mpi_barrier(mpi_comm_world,ierr)
# endif
    ! THE IOPROC SHOULD NOT TRY TO READ DATA!
    IF( MYID == IOPROCID) RETURN

    ! NOW TEST READING THE DATA BACK
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "TEST READING DATA BACK IN PARALLEL"
    NCFR=>NEW_FILE()
    NCFR%FNAME=NCF%FNAME
    CALL NC_OPEN(NCFR)

    CALL NC_LOAD(NCFR)

    CALL PRINT_FILE(NCFR)

    DIM_1r=>FIND_DIM(NCFR,"dim1",FOUND)
    CALL PRINT_DIM(DIM_1R)

    DIM_2r=>FIND_DIM(NCFR,"dim2",FOUND)
    CALL PRINT_DIM(DIM_2R)

    DIM_3r=>FIND_DIM(NCFR,"dim3",FOUND)
    CALL PRINT_DIM(DIM_3R)

    DIM_4r=>FIND_DIM(NCFR,"dim4",FOUND)
    CALL PRINT_DIM(DIM_4R)

    DIM_5r=>FIND_DIM(NCFR,"dim5",FOUND)
    CALL PRINT_DIM(DIM_5R)

    DIM_6r=>FIND_DIM(NCFR,"dim6",FOUND)
    CALL PRINT_DIM(DIM_6R)

    ! SCALARS
    VAR => find_var(NCFR,"flt",found)
    CALL NC_CONNECT_AVAR(VAR,fltr)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "flt=",fltr

    VAR => find_var(NCFR,"dbl",found)
    CALL NC_CONNECT_AVAR(VAR,dblr)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "dbl=",dblr

    VAR => find_var(NCFR,"int",found)
    CALL NC_CONNECT_AVAR(VAR,intr)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "int=",intr

    ! UNLIMITED SCLS
    VAR => find_var(NCFR,"u_flt",found)
    CALL NC_CONNECT_AVAR(VAR,fltr)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_flt=",fltr

    VAR => find_var(NCFR,"u_dbl",found)
    CALL NC_CONNECT_AVAR(VAR,dblr)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_dbl=",dblr

    VAR => find_var(NCFR,"u_int",found)
    CALL NC_CONNECT_AVAR(VAR,intr)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_int=",intr

    ! VECTORS
    allocate(flt1r(MT)); fltr=0.0_SP
    allocate(dbl1r(n)); dblr=0.0_DP
    allocate(int1r(0:N)); intr=0

    VAR => find_var(NCFR,"flt1",found)
    CALL NC_CONNECT_AVAR(VAR,flt1r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "flt1(MT)=",flt1r(MT)

    VAR => find_var(NCFR,"dbl1",found)
    CALL NC_CONNECT_AVAR(VAR,dbl1r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "dbl1(1:5)=",dbl1r(1:5)

    VAR => find_var(NCFR,"int1",found)
    CALL NC_CONNECT_AVAR(VAR,int1r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "int1(0:5)=",int1r(0:5)

    ! UNLIMITED VECTORS
    VAR => find_var(NCFR,"u_flt1",found)
    CALL NC_CONNECT_AVAR(VAR,flt1r)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_flt1(1)=",flt1r(1)

    VAR => find_var(NCFR,"u_dbl1",found)
    CALL NC_CONNECT_AVAR(VAR,dbl1r)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_dbl1(1:5)=",dbl1r(1:5)

    VAR => find_var(NCFR,"u_int1",found)
    CALL NC_CONNECT_AVAR(VAR,int1r)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_int1(0:5)=",int1r(0:5)

    ! ARRAYS
    allocate(flt2r(0:MT,3)); flt2r=0.0_sp
    allocate(dbl2r(0:N,3)); dbl2r=0.0_dp
    allocate(int2r(1:NT,4)); int2r=0

    VAR => find_var(NCFR,"flt2",found)
    CALL NC_CONNECT_AVAR(VAR,flt2r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "flt2(0:5,1)=",flt2r(0:5,1)

    VAR => find_var(NCFR,"int2",found)
    CALL NC_CONNECT_AVAR(VAR,int2r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "int2(0:5,1)=",int2r(0:5,1)

    VAR => find_var(NCFR,"dbl2",found)
    CALL NC_CONNECT_AVAR(VAR,dbl2r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "dbl2(NT,1)=",dbl2r(N,1)

    ! UNLIMITED ARRAYS
    VAR => find_var(NCFR,"u_flt2",found)
    CALL NC_CONNECT_AVAR(VAR,flt2r)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_flt2(0,:)=",flt2r(0,:)

    VAR => find_var(NCFR,"u_int2",found)
    CALL NC_CONNECT_AVAR(VAR,int2r)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_int2(1,:)=",int2r(1,:)

    VAR => find_var(NCFR,"u_dbl2",found)
    CALL NC_CONNECT_AVAR(VAR,dbl2r)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_dbl2(N,:)=",dbl2r(N,:)

    ! CUBES
    allocate(flt3r(0:NT,3,3))
    allocate(dbl3r(1:N,3,3))
    allocate(int3r(0:MT,3,4))
    VAR => find_var(NCFR,"flt3",found)
    CALL NC_CONNECT_AVAR(VAR,flt3r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "flt3(1,:,:)=",flt3r(1,:,:)

    VAR => find_var(NCFR,"dbl3",found)
    CALL NC_CONNECT_AVAR(VAR,dbl3r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "dbl3(N,:,:)=",dbl3r(N,:,:)

    VAR => find_var(NCFR,"int3",found)
    CALL NC_CONNECT_AVAR(VAR,int3r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "int3(MT,:,:)=",int3r(MT,:,:)

    ! UNLIMITED CUBES
    VAR => find_var(NCFR,"u_flt3",found)
    CALL NC_CONNECT_AVAR(VAR,flt3r)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_flt3(1,:,:)=",flt3r(1,:,:)

    VAR => find_var(NCFR,"u_dbl3",found)
    CALL NC_CONNECT_AVAR(VAR,dbl3r)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_dbl3(N,:,:)=",dbl3r(N,:,:)

    VAR => find_var(NCFR,"u_int3",found)
    CALL NC_CONNECT_AVAR(VAR,int3r)
    call nc_read_var(VAR,STKCNT=2)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "u_int3(MT,:,:)=",int3r(MT,:,:)

    ! CHARCTER STRINGS
    VAR => find_var(NCFR,"short_string",found)
    CALL NC_CONNECT_AVAR(VAR,str10r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "str10=",trim(str10r)

    VAR => find_var(NCFR,"string",found)
    CALL NC_CONNECT_AVAR(VAR,str20r)
    call nc_read_var(VAR)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "str20=",trim(str20r)

    ALLOCATE(str_vec10r(3),str_vec20r(3))
    VAR => find_var(NCFR,"vec_short_string",found)
    CALL NC_CONNECT_AVAR(VAR,str_vec10r)
    call nc_read_var(VAR)
    DO I = 1,3
       IF (DBG_SET(DBG_LOG)) write(ipt,*) "str_vec10="&
            &,trim(str_vec10r(I))
    END DO

    VAR => find_var(NCFR,"vec_string",found)
    CALL NC_CONNECT_AVAR(VAR,str_vec20r)
    call nc_read_var(VAR)
    DO I = 1,3
       IF (DBG_SET(DBG_LOG)) write(ipt,*) "str_vec20=",trim(str_vec20r(I))
    END DO



  END SUBROUTINE TEST_UNSTRUCTURED_IO


  SUBROUTINE TEST_START_COUNT_STRIDE
    USE NETCDF
    IMPLICIT NONE

       CHARACTER(LEN=80) :: FNAME
    CHARACTER(LEN=4) :: procnum

    integer:: i,j,k,ierr,status
    logical :: FOUND

    ! NC POINTERS
    TYPE(NCDIM), POINTER :: DIM_1,DIM_1r
    TYPE(NCDIM), POINTER :: DIM_2,DIM_2r
    TYPE(NCDIM), POINTER :: DIM_3,DIM_3r
    TYPE(NCDIM), POINTER :: DIM_4,DIM_4r
    TYPE(NCDIM), POINTER :: DIM_5,DIM_5r
    TYPE(NCDIM), POINTER :: DIM_6,DIM_6r
    TYPE(NCDIM), POINTER :: DIM_7,DIM_7r

    TYPE(NCATT), POINTER :: ATT
    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCFILE), POINTER :: NCFR

    INTEGER, allocatable :: strt(:), cnt(:), strd(:)

    ! SOME DUMMY VARIABLES
    INTEGER :: int, intr
    REAL(SP) :: flt, fltr
    REAL(DP) :: dbl, dblr

    CHARACTER(LEN=80) :: str20,str10, str20r, str10r
    CHARACTER(LEN=80), ALLOCATABLE :: str_vec10(:),str_vec20(:), str_vec10r(:),str_vec20r(:)
    
    INTEGER, ALLOCATABLE :: int1(:),int1r(:)
    REAL(SP), ALLOCATABLE :: flt1(:),flt1r(:)
    REAL(DP), ALLOCATABLE :: dbl1(:),dbl1r(:)

    INTEGER, ALLOCATABLE :: int2(:,:),int2r(:,:)
    REAL(SP), ALLOCATABLE :: flt2(:,:),flt2r(:,:)
    REAL(DP), ALLOCATABLE :: dbl2(:,:),dbl2r(:,:)

    INTEGER, ALLOCATABLE :: int3(:,:,:),int3r(:,:,:)
    REAL(SP), ALLOCATABLE :: flt3(:,:,:),flt3r(:,:,:)
    REAL(DP), ALLOCATABLE :: dbl3(:,:,:),dbl3r(:,:,:)
 


     IF(IOPROCID == MYID) THEN

       allocate(flt1(mgl));
       allocate(dbl1(ngl));
       allocate(int1(ngl));
       
       allocate(flt2(mgl,3));
       allocate(dbl2(ngl,3));
       allocate(int2(ngl,4)); 
       
       allocate(flt3(ngl,3,3));
       allocate(dbl3(ngl,3,3));
       allocate(int3(mgl,3,4))

    ELSE
       
       flt = 3.14159265358979323846_SP
       dbl = 3.14159265358979323846_DP
       int = nprocs
       
       allocate(flt1(mt)); flt1=2.0_SP
       allocate(dbl1(n)); dbl1=1.0_DP
       allocate(int1(0:N)); int1=MYID
       
       allocate(flt2(0:MT,2)); flt2=1.5
       allocate(dbl2(0:N,2)); dbl2=3.2_DP
       allocate(int2(1:NT,2)); int2=-100
       
       allocate(flt3(0:NT,2,2)); flt3=-1.5
       allocate(dbl3(1:N,2,2)); dbl3=3.0_DP
       allocate(int3(0:MT,2,2))
       DO i =1,MT
          DO j = 1,2
             DO k = 1,2
                int3(i,j,k) = NGID_X(i)*j*k
             END DO
          END DO
       END DO
       

    END IF


    NCF => NEW_FILE()

    FNAME= "./results/FTEST3.nc"
    IF (DBG_SET(DBG_LOG))  write(ipt,*) "WRITING FILE: "//TRIM(FNAME)
    NCF%FNAME=FNAME

    ATT => NC_MAKE_ATT(name='title',values="A TEST FILE for start coun&
         &t stride data") 
    NCF => ADD(NCF,ATT)


    DIM_1       => NC_MAKE_DIM(name='dim1',len=MGL)

    DIM_2       => NC_MAKE_DIM(name='dim2',len=NGL)

    DIM_3       => NC_MAKE_DIM(name='dim3',len=3)

    DIM_4       => NC_MAKE_DIM(name='dim4',len=4)

    DIM_5       => NC_MAKE_DIM(name='dim5',len=NF90_UNLIMITED)

    DIM_6       => NC_MAKE_DIM(name='dim6',len=20)

    DIM_7       => NC_MAKE_DIM(name='dim7',len=2)


    ! SCALAR VARIABLES
    VAR  => NC_MAKE_AVAR(name='flt',values=flt)
    VAR => add(VAR,COPY_DIM(DIM_3))
    VAR => add(VAR,COPY_DIM(DIM_4))
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='dbl',values=dbl)
    VAR => add(VAR,COPY_DIM(DIM_3))
    VAR => add(VAR,COPY_DIM(DIM_4))
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='int',values=int)
    VAR => add(VAR,COPY_DIM(DIM_3))
    VAR => add(VAR,COPY_DIM(DIM_7))
    NCF  => ADD(NCF,VAR)

    ! VECTOR VARIABLES
    VAR  => NC_MAKE_AVAR(name='flt1',values=flt1,DIM1=DIM_1)
    VAR => add(VAR,COPY_DIM(DIM_3))
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='dbl1',values=dbl1,DIM1=DIM_2)
    VAR => add(VAR,COPY_DIM(DIM_3))
    VAR => add(VAR,COPY_DIM(DIM_4))
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='int1',values=int1,DIM1=DIM_2)
    VAR => add(VAR,COPY_DIM(DIM_3))
    VAR => add(VAR,COPY_DIM(DIM_7))
    NCF  => ADD(NCF,VAR)

    ! ARRAY VARIABLES
    VAR  => NC_MAKE_AVAR(name='flt2',values=flt2,DIM1=DIM_1,DIM2=DIM_7)
    VAR => add(VAR,COPY_DIM(DIM_3))
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='dbl2',values=dbl2,DIM1=DIM_2,DIM2=DIM_7)
    VAR => add(VAR,COPY_DIM(DIM_7))
    VAR => add(VAR,COPY_DIM(DIM_3))
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='int2',values=int2,DIM1=DIM_2,DIM2=DIM_7)
    VAR => add(VAR,COPY_DIM(DIM_7))
    VAR => add(VAR,COPY_DIM(DIM_7))
    NCF  => ADD(NCF,VAR)

    ! CUBE VARIABLES
    VAR  => NC_MAKE_AVAR(name='flt3',values=flt3,DIM1=DIM_2,DIM2=DIM_7,DIM3=DIM_7)
    VAR => add(VAR,COPY_DIM(DIM_3))
    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='dbl3',values=dbl3,DIM1=DIM_2,DIM2=DIM_7,DIM3=DIM_7)
    VAR => add(VAR,COPY_DIM(DIM_7))

    NCF  => ADD(NCF,VAR)

    VAR  => NC_MAKE_AVAR(name='int3',values=int3,DIM1=DIM_1,DIM2=DIM_7,DIM3=DIM_7)
    VAR => add(VAR,COPY_DIM(DIM_7))
    NCF  => ADD(NCF,VAR)
    
    IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! CREATING FILE AND WRITING DATA:"

    IF(MSR)THEN
       CALL NC_CREATE(NCF)
       CALL NC_SAVE(NCF)
    END IF
    

    ALLOCATE(STRT(2),CNT(2),STRD(2))

    strt=1
    cnt=1
    strd=1
    VAR => FIND_VAR(NCF,"flt",FOUND)
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

    strt(1)=3
    strt(2)=3
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)
    
    VAR => FIND_VAR(NCF,"dbl",FOUND)
    strt(1)=3
    strt(2)=4
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

    VAR => FIND_VAR(NCF,"int",FOUND)
    strt(1)=1
    strt(2)=1
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

    VAR => FIND_VAR(NCF,"flt1",FOUND)
    strt(1)=1
    strt(2)=2
    cnt(1)=MGL
    cnt(2)=1
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)
   

    DEALLOCATE(STRT,CNT,STRD)
    ALLOCATE(STRT(3),CNT(3),STRD(3))

    strd=1

    VAR => FIND_VAR(NCF,"dbl1",FOUND)
    strt(1)=1
    strt(2)=1
    strt(3)=1
    cnt(1)=NGL
    cnt(2)=1
    cnt(3)=1
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

    VAR => FIND_VAR(NCF,"int1",FOUND)
    strt(1)=1
    strt(2)=1
    strt(3)=2
    cnt(1)=NGL
    cnt(2)=1
    cnt(3)=1
    ! TESTING NO STRIDE ARGUMENT
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT)
    strt(2)=3
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT)

    VAR => FIND_VAR(NCF,"flt2",FOUND)
    STRD(1)=1
    STRD(2)=2
    STRD(3)=2
    strt(1)=1
    strt(2)=2
    strt(3)=1
    cnt(1)=MGL
    cnt(2)=1
    cnt(3)=2

    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

    DEALLOCATE(STRT,CNT,STRD)
    ALLOCATE(STRT(4),CNT(4),STRD(4))

    VAR => FIND_VAR(NCF,"dbl2",FOUND)
    STRD(1)=1
    STRD(2)=1
    STRD(3)=1
    STRD(4)=2
   
    strt(1)=1
    strt(2)=1
    strt(3)=1
    strt(4)=1
   
    cnt(1)=NGL
    cnt(2)=1
    cnt(3)=1
    cnt(4)=2

    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

    VAR => FIND_VAR(NCF,"int2",FOUND)
    STRD(1)=1
    STRD(2)=1
    STRD(3)=1
    STRD(4)=1
   
    strt(1)=1
    strt(2)=1
    strt(3)=1
    strt(4)=1
   
    cnt(1)=NGL
    cnt(2)=2
    cnt(3)=1
    cnt(4)=1

    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

    VAR => FIND_VAR(NCF,"flt3",FOUND)
    STRD(1)=1
    STRD(2)=1
    STRD(3)=1
    STRD(4)=1
   
    strt(1)=1
    strt(2)=2
    strt(3)=1
    strt(4)=1
   
    cnt(1)=NGL
    cnt(2)=1
    cnt(3)=2
    cnt(4)=2

    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

    VAR => FIND_VAR(NCF,"dbl3",FOUND)
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

    VAR => FIND_VAR(NCF,"int3",FOUND)
    cnt(1)=MGL
    CALL NC_WRITE_VAR(VAR,MSR,PAR,MSRID,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD)

# if defined(MULTIPROCESSOR)
    call mpi_barrier(mpi_FVCOM_group,status)
# endif
    CALL NC_OPEN(NCF)


    deallocate(int2)
    Allocate(int2(2,2)); int2=0

    STRD(1)=1
    STRD(2)=1
    STRD(3)=1
    STRD(4)=1
   
    strt(1)=5
    strt(2)=2
    strt(3)=1
    strt(4)=1
   
    cnt(1)=1
    cnt(2)=1
    cnt(3)=2
    cnt(4)=2

    CALL NC_CONNECT_AVAR(VAR,int2)
    call nc_read_var(VAR,IOSTART=STRT,IOCOUNT=CNT,IOSTRIDE=STRD,DEALERID=MYID,PARALLEL=.FALSE.)
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "int2=",int2


    ! SHUOLD ADD MORE TESTING HERE, BUT THIS IS A PAIN IN THE NECK!
    IF (DBG_SET(DBG_LOG)) write(ipt,*) "! SHUOLD ADD MORE TESTING HERE, BUT THIS IS A PAIN IN THE NECK!" 

  END SUBROUTINE TEST_START_COUNT_STRIDE

END MODULE MOD_NCTEST
