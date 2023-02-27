MODULE MOD_PARALLEL
  USE MOD_UTILS
  USE MOD_INPUT
  USE MOD_PAR
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

          if (opt_sng == "GRIDFILE" .or.opt_sng == "GridFile"&
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


    ! Special Setting:
    dbg_lvl=5
    dbg_par=.true.

    CALL dbg_init(IPT_BASE,.false.)

  END SUBROUTINE GET_COMMANDLINE

  SUBROUTINE MYHELPTXT
    IMPLICIT NONE

    WRITE(IPT,*) "! WELCOME TO THE NCTOOLS TEST SUIT"
    WRITE(IPT,*) "! OPTIONS:"
    WRITE(IPT,*) "! --GRIDFILE=<a valid FVCOM GRID FILE NAME>"
    WRITE(IPT,*) "! --USE_MPI_IO"
    WRITE(IPT,*) "!   This option test paralle communication used in MPI_IO"
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

    close(GRIDUNIT)

    CALL SETUP_DOMAIN

  END SUBROUTINE GET_FVCOM_GRID


  SUBROUTINE PAR_TEST
    IMPLICIT NONE

#if defined(MULTIPROCESSOR)
    INTEGER ::  SENDID,I,RECVID,J

    REAL(SPA), POINTER, DIMENSION(:) :: vec_flt_GBL
    REAL(SPA), POINTER, DIMENSION(:) :: vec_flt_LCL

    REAL(SPA), POINTER, DIMENSION(:,:) :: arr_flt_GBL
    REAL(SPA), POINTER, DIMENSION(:,:) :: arr_flt_LCL

    REAL(DP), POINTER, DIMENSION(:) :: vec_dbl_GBL
    REAL(DP), POINTER, DIMENSION(:) :: vec_dbl_LCL

    REAL(DP), POINTER, DIMENSION(:,:) :: arr_dbl_GBL
    REAL(DP), POINTER, DIMENSION(:,:) :: arr_dbl_LCL

    INTEGER, POINTER, DIMENSION(:) :: vec_int_GBL
    INTEGER, POINTER, DIMENSION(:) :: vec_int_LCL

    INTEGER, POINTER, DIMENSION(:,:) :: arr_int_GBL
    INTEGER, POINTER, DIMENSION(:,:) :: arr_int_LCL


    !==============================================================================|
    !   TEST DEAL AND COLLECT  - debug mode only!                                  | 
    !==============================================================================|

    WRITE(IPT,*)  '!  NODE INT_VEC_DEAL TEST        :     '

    allocate(vec_int_GBL(0:MGL))
    DO I=1,MGL
       vec_int_GBL(I)=I
    END DO
    allocate(vec_int_LCL(0:M+NHN)); vec_int_LCL=0

    SENDID =1
    if (USE_MPI_IO_MODE) SENDID = IOPROCID ! TEST DEAL FROM IOPROC
    CALL PDEAL(MYID,SENDID,NPROCS,NXMAP,VEC_INT_GBL,VEC_INT_LCL)


    if(.NOT. IOPROC) then ! TEST
       DO I = 1,M
          IF(VEC_INT_GBL(NGID(I)) .NE. VEC_INT_LCL(I))&
               & CALL FATAL_ERROR("VEC_INT_DEAL TEST: INTERNAL NODES D&
               &O NOT MATCH")
       END DO
       DO I = 1,NHN
          IF(VEC_INT_GBL(HN_LST(I)) .NE. VEC_INT_LCL(I+M))&
               &CALL FATAL_ERROR("VEC_INT_DEAL HALO NODES IS BROKEN")
       END DO
    end if


    WRITE(IPT,*)  '!  NODE INT_VEC_DEAL TEST        :    PASSED    '
    WRITE(IPT,*)  '!  NODE INT_VEC_COLLECT TEST     :     '

    ! NOW PASS IT BACK TO THE GLOBAL ON EACH PROCESSOR

    VEC_INT_GBL = 0
    DO I = 1, NPROCS_total
       RECVID= I
       if( .not. IOPROC .or. (RECVID .EQ. MYID)) &
            CALL PCOLLECT(MYID,RECVID,NPROCS,NXMAP,VEC_INT_LCL,VEC_INT_GBL)
    END DO

    write(ipt,*)"====================================="
    DO I = 1,MGL

       IF( I .NE. VEC_INT_GBL(I)) &
            & CALL FATAL_ERROR("VEC_INT_COLLECT NODE TEST: FAILED")
    END DO

    DEALLOCATE(VEC_INT_GBL)
    DEALLOCATE(VEC_INT_LCL)

    WRITE(IPT,*)  '!  NODE INT_VEC_COLLECT TEST     :    PASSED  '
    WRITE(IPT,*)  '!  ELEMENT INT_VEC_DEL TEST      :     '


    allocate(vec_int_GBL(NGL))
    vec_int_GBL=EL_PID
    allocate(vec_int_LCL(N+NHE))

    SENDID =1
    if (USE_MPI_IO_MODE) SENDID = IOPROCID ! TEST DEAL FROM IOPROC
    CALL PDEAL(MYID,SENDID,NPROCS,EXMAP,VEC_INT_GBL,VEC_INT_LCL)

    if(.NOT. IOPROC) then ! TEST
       DO I = 1,N
          IF(EL_PID(EGID(I)) .NE. VEC_INT_LCL(I))&
               & CALL FATAL_ERROR("VEC_INT_DEAL TEST: INTERNAL ELEMENTS D&
               &O NOT MATCH")
       END DO
       DO I = 1,NHE
          IF(EL_PID(HE_LST(I)) .NE. VEC_INT_LCL(I+N))&
               & CALL FATAL_ERROR("VEC_INT_DEAL HALO ELEMENTS IS BROKEN")
       END DO
    end if


    WRITE(IPT,*)  '!  ELEMENT INT_VEC_DEAL TEST     :    PASSED    '
    WRITE(IPT,*)  '!  ELEMENT INT_VEC_COLLECT TEST  :     '

    ! NOW PASS IT BACK TO THE GLOBAL ON EACH PROCESSOR
    VEC_INT_GBL = 0
    DO I = 1, NPROCS_total
       RECVID= I
       if( .not. IOPROC .or. (RECVID .EQ. MYID)) &
            & CALL PCOLLECT(MYID,RECVID,NPROCS,EXMAP,VEC_INT_LCL,VEC_INT_GBL)
    END DO


    DO I = 1,NGL
       IF(EL_PID(I) .NE. VEC_INT_GBL(I)) &
            & CALL FATAL_ERROR("VEC_INT_COLLECT TEST: FAILED")
    END DO

    DEALLOCATE(VEC_INT_GBL)
    DEALLOCATE(VEC_INT_LCL)
    WRITE(IPT,*)  '!  ELEMENT INT_VEC_COLLECT TEST  :    PASSED '
    WRITE(IPT,*)  '!  NODE INT_ARR_DEAL TEST        :     '


    allocate(ARR_int_GBL(MGL,KB))
    DO I = 1,MGL
       DO J=1,KB
          ARR_int_GBL(I,J)=J*I
       END DO
    END DO

    allocate(ARR_int_LCL(M+NHN,KB))

    SENDID =1
    if (USE_MPI_IO_MODE) SENDID = IOPROCID ! TEST DEAL FROM IOPROC
    CALL PDEAL(MYID,SENDID,NPROCS,NXMAP,ARR_INT_GBL,ARR_INT_LCL)

    IF (.NOT. IOPROC) THEN ! TEST
       DO I = 1,M
          DO J =1,KB
             IF(ARR_INT_GBL(NGID(I),J) .NE. ARR_INT_LCL(I,J))&
                  & CALL FATAL_ERROR("ARR_INT_DEAL TEST: INTERNAL ELEMENTS D&
                  &O NOT MATCH")
          END DO
       END DO
       DO I = 1,NHN
          DO J = 1,KB
             IF(ARR_INT_GBL(HN_LST(I),J) .NE. ARR_INT_LCL(I+M,J))&
                  & CALL FATAL_ERROR("ARR_INT_DEAL HALO ELEMENTS IS BROKEN")
          END DO
       END DO
    END IF


    WRITE(IPT,*)  '!  NODE INT_ARR_DEAL TEST        :    PASSED    '
    WRITE(IPT,*)  '!  NODE INT_ARR_COLELCT TEST     :    '

    ! NOW PASS IT BACK TO THE GLOBAL ON EACH PROCESSOR
    ARR_INT_GBL = 0
    DO I = 1, NPROCS_total
       RECVID = I
       if( .not. IOPROC .or. (RECVID .EQ. MYID)) &
            &CALL PCOLLECT(MYID,RECVID,NPROCS,NXMAP,ARR_INT_LCL,ARR_INT_GBL)
    END DO

    DO I = 1,MGL
       DO J = 1,KB
          IF(I*J .NE. ARR_INT_GBL(I,j)) &
               & CALL FATAL_ERROR("ARR_INT_COLLECT TEST: FAILED")
       END DO
    END DO


    DEALLOCATE(ARR_INT_GBL)
    DEALLOCATE(ARR_INT_LCL)

    WRITE(IPT,*)  '!  NODE INT_ARR_COLLECT TEST     :    PASSED '
    WRITE(IPT,*)  '!  ELEMENT INT_ARR_DEAL TEST     :     '

    allocate(ARR_int_GBL(NGL,KB))
    DO I = 1,KB
       ARR_int_GBL(:,I)=EL_PID*I
    END DO

    allocate(ARR_int_LCL(N+NHE,KB))

    SENDID =1
    if (USE_MPI_IO_MODE) SENDID = IOPROCID ! TEST DEAL FROM IOPROC
    CALL PDEAL(MYID,SENDID,NPROCS,EXMAP,ARR_INT_GBL,ARR_INT_LCL)

    IF (.NOT. IOPROC) THEN ! TEST
       DO I = 1,N
          DO J =1,KB
             IF(EL_PID(EGID(I))*J .NE. ARR_INT_LCL(I,J))&
                  & CALL FATAL_ERROR("ARR_INT_DEAL TEST: INTERNAL ELEMENTS D&
                  &O NOT MATCH")
          END DO
       END DO
       DO I = 1,NHE
          DO J = 1,KB
             IF(EL_PID(HE_LST(I))*J .NE. ARR_INT_LCL(I+N,J))&
                  & CALL FATAL_ERROR("ARR_INT_DEAL HALO ELEMENTS IS BROKEN")
          END DO
       END DO
    END IF


    WRITE(IPT,*)  '!  ELEMENT INT_ARR_DEAL TEST     :    PASSED    '
    WRITE(IPT,*)  '!  ELEMENT INT_ARR_COLLECT TEST  :   '

    ! NOW PASS IT BACK TO THE GLOBAL ON EACH PROCESSOR
    ARR_INT_GBL = 0
    DO I = 1, NPROCS_total
       RECVID = I
       if( .not. IOPROC .or. (RECVID .EQ. MYID)) &
            &CALL PCOLLECT(MYID,RECVID,NPROCS,EXMAP,ARR_INT_LCL,ARR_INT_GBL)
    END DO

    DO I = 1,NGL
       DO J = 1,KB
          IF(EL_PID(I)*J .NE. ARR_INT_GBL(I,j)) &
               & CALL FATAL_ERROR("ARR_INT_COLLECT TEST: FAILED")
       END DO
    END DO


    DEALLOCATE(ARR_INT_GBL)
    DEALLOCATE(ARR_INT_LCL)

    WRITE(IPT,*)  '!  ELEMENT INT_ARR_COLLECT TEST  :    PASSED '


    ! TEST ARRAY COLLECT AND DEAL for FLOATS and DOUBLES!

    WRITE(IPT,*)  '!  ELEMENT FLT_VEC_DEAL TEST     :  '

    allocate(vec_int_GBL(NGL))
    vec_int_GBL=EL_PID
    allocate(vec_int_LCL(N+NHE))

    allocate(vec_FLT_GBL(NGL))
    vec_FLT_GBL=REAL(VEC_INT_GBL,SPA)
    allocate(vec_FLT_LCL(N+NHE))

    SENDID =1
    if (USE_MPI_IO_MODE) SENDID = IOPROCID ! TEST DEAL FROM IOPROC
    CALL PDEAL(MYID,SENDID,NPROCS,EXMAP,VEC_FLT_GBL,VEC_FLT_LCL)

    if(.NOT. IOPROC) then ! TEST
       VEC_INT_LCL=ANINT(VEC_FLT_LCL)
       DO I = 1,N
          IF(EL_PID(EGID(I)) .NE. VEC_INT_LCL(I))&
               & CALL FATAL_ERROR("VEC_FLOAT_DEAL TEST: INTERNAL ELEMENTS D&
               &O NOT MATCH")
       END DO
       DO I = 1,NHE
          IF(EL_PID(HE_LST(I)) .NE. VEC_INT_LCL(I+N))&
               & CALL FATAL_ERROR("VEC_FLOAT_DEAL HALO ELEMENTS IS BROKEN")
       END DO
    end if

    WRITE(IPT,*)  '!  ELEMENT FLT_VEC_DEAL          :    PASSED    '
    WRITE(IPT,*)  '!  ELEMENT FLT_VEC_COLLECT       :    '
    !   IF(DBG_SET(DBG_LOG))WRITE(IPT,*)  '!  ELEMENT COLLECT TEST  :     '

    ! NOW PASS IT BACK TO THE GLOBAL ON EACH PROCESSOR
    VEC_FLT_GBL = 0.0
    VEC_INT_GBL = 0
    DO I = 1, NPROCS_total
       RECVID= I
       if( .not. IOPROC .or. (RECVID .EQ. MYID)) &
            & CALL PCOLLECT(MYID,RECVID,NPROCS,EXMAP,VEC_FLT_LCL,VEC_FLT_GBL)
    END DO

    VEC_INT_GBL = ANINT(VEC_FLT_GBL)

    DO I = 1,NGL
       IF(EL_PID(I) .NE. VEC_INT_GBL(I)) &
            & CALL FATAL_ERROR("VEC_FLOAT_COLLECT TEST: FAILED")
    END DO


    DEALLOCATE(VEC_INT_GBL)
    DEALLOCATE(VEC_INT_LCL)
    DEALLOCATE(VEC_FLT_GBL)
    DEALLOCATE(VEC_FLT_LCL)


    WRITE(IPT,*)  '!  ELEMENT FLT_VEC_COLLECT TEST  :    PASSED '
    WRITE(IPT,*)  '!  ELEMENT FLT_ARR_DEAL TEST     :     '


    allocate(ARR_int_GBL(NGL,KB))
    DO I = 1,KB
       ARR_int_GBL(:,I)=EL_PID*I
    END DO
    allocate(ARR_int_LCL(N+NHE,KB))

    allocate(ARR_FLT_GBL(NGL,KB))
    DO I = 1,KB
       ARR_FLT_GBL(:,I)=REAL(ARR_int_GBL(:,I),SPA)
    END DO
    allocate(ARR_FLT_LCL(N+NHE,KB))


    SENDID =1
    if (USE_MPI_IO_MODE) SENDID = IOPROCID ! TEST DEAL FROM IOPROC
    CALL PDEAL(MYID,SENDID,NPROCS,EXMAP,ARR_FLT_GBL,ARR_FLT_LCL)

    IF (.NOT. IOPROC) THEN ! TEST
       ARR_INT_LCL=ANINT(ARR_FLT_LCL)
       DO I = 1,N
          DO J =1,KB
             IF(EL_PID(EGID(I))*J .NE. ARR_INT_LCL(I,J))&
                  & CALL FATAL_ERROR("ARR_FLOAT_DEAL TEST: INTERNAL ELEMENTS D&
                  &O NOT MATCH")
          END DO
       END DO
       DO I = 1,NHE
          DO J = 1,KB
             IF(EL_PID(HE_LST(I))*J .NE. ARR_INT_LCL(I+N,J))&
                  & CALL FATAL_ERROR("ARR_FLOAT_DEAL HALO ELEMENTS IS BROKEN")
          END DO
       END DO
    END IF


    WRITE(IPT,*)  '!  ELEMENT FLT_ARR_DEAL TEST     :    PASSED    '
    WRITE(IPT,*)  '!  ELEMENT FLT_ARR_COLLECT TEST  :  '

    ! NOW PASS IT BACK TO THE GLOBAL ON EACH PROCESSOR
    ARR_INT_GBL = 0
    ARR_FLT_GBL = 0.0
    DO I = 1, NPROCS_total
       RECVID = I
       if( .not. IOPROC .or. (RECVID .EQ. MYID)) &
            &CALL PCOLLECT(MYID,RECVID,NPROCS,EXMAP,ARR_FLT_LCL,ARR_FLT_GBL)
    END DO
    ARR_INT_GBL =ANINT(ARR_FLT_GBL)

    DO I = 1,NGL
       DO J = 1,KB
          IF(EL_PID(I)*J .NE. ARR_INT_GBL(I,j)) &
               & CALL FATAL_ERROR("ARR_FLOAT_COLLECT TEST: FAILED")
       END DO
    END DO


    DEALLOCATE(ARR_INT_GBL)
    DEALLOCATE(ARR_INT_LCL)
    DEALLOCATE(ARR_FLT_GBL)
    DEALLOCATE(ARR_FLT_LCL)

    WRITE(IPT,*)  '!  ELEMENT FLT_ARR_COLLECT TEST  :    PASSED '
    WRITE(IPT,*)  '!  ELEMENT FLT_DBL_DEAL TEST     :  '

    allocate(vec_int_GBL(NGL))
    vec_int_GBL=EL_PID
    allocate(vec_int_LCL(N+NHE))

    allocate(vec_DBL_GBL(NGL))
    vec_DBL_GBL=REAL(VEC_INT_GBL,DP)
    allocate(vec_DBL_LCL(N+NHE))

    SENDID =1
    if (USE_MPI_IO_MODE) SENDID = IOPROCID ! TEST DEAL FROM IOPROC
    CALL PDEAL(MYID,SENDID,NPROCS,EXMAP,VEC_DBL_GBL,VEC_DBL_LCL)

    if(.NOT. IOPROC) then ! TEST
       VEC_INT_LCL=ANINT(VEC_DBL_LCL)
       DO I = 1,N
          IF(EL_PID(EGID(I)) .NE. VEC_INT_LCL(I))&
               & CALL FATAL_ERROR("VEC_DBL_DEAL TEST: INTERNAL ELEMENTS D&
               &O NOT MATCH")
       END DO
       DO I = 1,NHE
          IF(EL_PID(HE_LST(I)) .NE. VEC_INT_LCL(I+N))&
               & CALL FATAL_ERROR("VEC_DBL_DEAL HALO ELEMENTS IS BROKEN")
       END DO
    end if


    WRITE(IPT,*)  '!  ELEMENT DBL_VEC_DEAL          :    PASSED    '
    WRITE(IPT,*)  '!  ELEMENT DBL_VEC_COLLECT       :    '

    ! NOW PASS IT BACK TO THE GLOBAL ON EACH PROCESSOR
    VEC_DBL_GBL = 0.0
    VEC_INT_GBL = 0
    DO I = 1, NPROCS_total
       RECVID= I
       if( .not. IOPROC .or. (RECVID .EQ. MYID)) &
            & CALL PCOLLECT(MYID,RECVID,NPROCS,EXMAP,VEC_DBL_LCL,VEC_DBL_GBL)
    END DO

    VEC_INT_GBL = ANINT(VEC_DBL_GBL)

    DO I = 1,NGL
       IF(EL_PID(I) .NE. VEC_INT_GBL(I)) &
            & CALL FATAL_ERROR("VEC_DBL_COLLECT TEST: FAILED")
    END DO


    DEALLOCATE(VEC_INT_GBL)
    DEALLOCATE(VEC_INT_LCL)
    DEALLOCATE(VEC_DBL_GBL)
    DEALLOCATE(VEC_DBL_LCL)

    WRITE(IPT,*)  '!  ELEMENT DBL_VEC_COLLECT TEST  :    PASSED '
    WRITE(IPT,*)  '!  ELEMENT DBL_ARR_DEAL TEST     :     '


    allocate(ARR_int_GBL(NGL,KB))
    DO I = 1,KB
       ARR_int_GBL(:,I)=EL_PID*I
    END DO
    allocate(ARR_int_LCL(N+NHE,KB))

    allocate(ARR_DBL_GBL(NGL,KB))
    DO I = 1,KB
       ARR_DBL_GBL(:,I)=REAL(ARR_int_GBL(:,I),DP)
    END DO
    allocate(ARR_DBL_LCL(N+NHE,KB))


    SENDID =1
    if (USE_MPI_IO_MODE) SENDID = IOPROCID ! TEST DEAL FROM IOPROC
    CALL PDEAL(MYID,SENDID,NPROCS,EXMAP,ARR_DBL_GBL,ARR_DBL_LCL)

    IF (.NOT. IOPROC) THEN ! TEST
       ARR_INT_LCL=ANINT(ARR_DBL_LCL)
       DO I = 1,N
          DO J =1,KB
             IF(EL_PID(EGID(I))*J .NE. ARR_INT_LCL(I,J))&
                  & CALL FATAL_ERROR("ARR_DBL_DEAL TEST: INTERNAL ELEMENTS D&
                  &O NOT MATCH")
          END DO
       END DO
       DO I = 1,NHE
          DO J = 1,KB
             IF(EL_PID(HE_LST(I))*J .NE. ARR_INT_LCL(I+N,J))&
                  & CALL FATAL_ERROR("ARR_DBL_DEAL HALO ELEMENTS IS BROKEN")
          END DO
       END DO
    END IF


    WRITE(IPT,*)  '!  ELEMENT DBL_ARR_DEAL TEST     :    PASSED    '
    WRITE(IPT,*)  '!  ELEMENT DBL_ARR_COLLECT TEST  :  '

    ! NOW PASS IT BACK TO THE GLOBAL ON EACH PROCESSOR
    ARR_INT_GBL = 0
    ARR_DBL_GBL = 0.0
    DO I = 1, NPROCS_total
       RECVID = I
       if( .not. IOPROC .or. (RECVID .EQ. MYID)) &
            &CALL PCOLLECT(MYID,RECVID,NPROCS,EXMAP,ARR_DBL_LCL,ARR_DBL_GBL)
    END DO
    ARR_INT_GBL =ANINT(ARR_DBL_GBL)

    DO I = 1,NGL
       DO J = 1,KB
          IF(EL_PID(I)*J .NE. ARR_INT_GBL(I,j)) &
               & CALL FATAL_ERROR("ARR_DBL_COLLECT TEST: FAILED")
       END DO
    END DO


    DEALLOCATE(ARR_INT_GBL)
    DEALLOCATE(ARR_INT_LCL)
    DEALLOCATE(ARR_DBL_GBL)
    DEALLOCATE(ARR_DBL_LCL)


    WRITE(IPT,*)  '!  ELEMENT DBL_ARR_COLLECT TEST  :    PASSED '

#else
    write(IPT,*) "!==============================================="
    write(IPT,*) "! MULTIPROCESSOR IS TURNED OFF IN THE MAKE FILE"
    write(IPT,*) "!==============================================="
# endif    
    
  END SUBROUTINE PAR_TEST
  


END MODULE MOD_PARALLEL
