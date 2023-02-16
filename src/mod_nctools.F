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

MODULE MOD_NCTOOLS
  USE MOD_NCLL
  USE LIMS
  USE MOD_PAR
  USE CONTROL, only: IPT, USE_REAL_WORLD_TIME,timeprec,datestrlen !, MSR, USE_MPI_IO_MODE, PAR, SERIAL, IOPROC
  implicit none
  save


  ! AT PRESENT THIS LIBRARY DOES NOT PROPERLY SUPPORT UDUNITS TIME
  ! CONVENTIONS. THE FOLLOWING ARE THE UNDERSTOOD UNITS ATTRIBUTES
  character(len=80), parameter :: seconds_units  ='seconds'
  character(len=80), parameter :: days_units  ='days since 0.0'
  character(len=80), parameter :: mjd_units  ='days since 1858-11-17 00:00:00'
  character(len=80), parameter :: msec_units ='msec since 00:00:00'
  character(len=80), parameter :: fmat ='modified julian day (MJD)'
  character(len=80), parameter :: rfmat ='defined reference date'


  INTERFACE NC_GET_ATT
     MODULE PROCEDURE NC_GET_VATT
     MODULE PROCEDURE NC_GET_GATT
  END INTERFACE

  INTERFACE NC_MAKE_RUNTIME_ATT_CHR
     MODULE PROCEDURE NC_MAKE_RUNTIME_ATT_CHR_SCL
     MODULE PROCEDURE NC_MAKE_RUNTIME_ATT_CHR_VEC
  END INTERFACE

  INTERFACE NC_MAKE_ATT
     MODULE PROCEDURE NC_MAKE_ATT_CHR_SCL
     MODULE PROCEDURE NC_MAKE_ATT_CHR_VEC
     MODULE PROCEDURE NC_MAKE_ATT_INT
     MODULE PROCEDURE NC_MAKE_ATT_INT_VEC
     MODULE PROCEDURE NC_MAKE_ATT_FLT
     MODULE PROCEDURE NC_MAKE_ATT_FLT_VEC
     MODULE PROCEDURE NC_MAKE_ATT_DBL
     MODULE PROCEDURE NC_MAKE_ATT_DBL_VEC
  END INTERFACE

  INTERFACE NC_MAKE_AVAR
     MODULE PROCEDURE NC_MAKE_AVAR_SCL_CHR
     MODULE PROCEDURE NC_MAKE_AVAR_VEC_CHR

     MODULE PROCEDURE NC_MAKE_AVAR_SCL_DBL
     MODULE PROCEDURE NC_MAKE_AVAR_VEC_DBL
     MODULE PROCEDURE NC_MAKE_AVAR_ARR_DBL
     MODULE PROCEDURE NC_MAKE_AVAR_CUB_DBL
     MODULE PROCEDURE NC_MAKE_AVAR_FDA_DBL

     MODULE PROCEDURE NC_MAKE_AVAR_SCL_FLT
     MODULE PROCEDURE NC_MAKE_AVAR_VEC_FLT
     MODULE PROCEDURE NC_MAKE_AVAR_ARR_FLT
     MODULE PROCEDURE NC_MAKE_AVAR_CUB_FLT
     MODULE PROCEDURE NC_MAKE_AVAR_FDA_FLT

     MODULE PROCEDURE NC_MAKE_AVAR_SCL_INT
     MODULE PROCEDURE NC_MAKE_AVAR_VEC_INT
     MODULE PROCEDURE NC_MAKE_AVAR_ARR_INT
     MODULE PROCEDURE NC_MAKE_AVAR_CUB_INT
     MODULE PROCEDURE NC_MAKE_AVAR_FDA_INT
  END INTERFACE

  INTERFACE NC_MAKE_PVAR
     MODULE PROCEDURE NC_MAKE_PVAR_SCL_CHR
     MODULE PROCEDURE NC_MAKE_PVAR_VEC_CHR

     MODULE PROCEDURE NC_MAKE_PVAR_SCL_DBL
     MODULE PROCEDURE NC_MAKE_PVAR_VEC_DBL
     MODULE PROCEDURE NC_MAKE_PVAR_ARR_DBL
     MODULE PROCEDURE NC_MAKE_PVAR_CUB_DBL
     MODULE PROCEDURE NC_MAKE_PVAR_FDA_DBL

     MODULE PROCEDURE NC_MAKE_PVAR_SCL_FLT
     MODULE PROCEDURE NC_MAKE_PVAR_VEC_FLT
     MODULE PROCEDURE NC_MAKE_PVAR_ARR_FLT
     MODULE PROCEDURE NC_MAKE_PVAR_CUB_FLT
     MODULE PROCEDURE NC_MAKE_PVAR_FDA_FLT

     MODULE PROCEDURE NC_MAKE_PVAR_SCL_INT
     MODULE PROCEDURE NC_MAKE_PVAR_VEC_INT
     MODULE PROCEDURE NC_MAKE_PVAR_ARR_INT
     MODULE PROCEDURE NC_MAKE_PVAR_CUB_INT
     MODULE PROCEDURE NC_MAKE_PVAR_FDA_INT
  END INTERFACE

  ! NC_DISCONNECT
  INTERFACE NC_CONNECT_AVAR
     MODULE PROCEDURE NC_CONNECT_VAR_SCL_CHR
     MODULE PROCEDURE NC_CONNECT_AVAR_VEC_CHR

     MODULE PROCEDURE NC_CONNECT_AVAR_SCL_DBL
     MODULE PROCEDURE NC_CONNECT_AVAR_VEC_DBL
     MODULE PROCEDURE NC_CONNECT_AVAR_ARR_DBL
     MODULE PROCEDURE NC_CONNECT_AVAR_CUB_DBL
     MODULE PROCEDURE NC_CONNECT_AVAR_FDA_DBL

     MODULE PROCEDURE NC_CONNECT_AVAR_SCL_FLT
     MODULE PROCEDURE NC_CONNECT_AVAR_VEC_FLT
     MODULE PROCEDURE NC_CONNECT_AVAR_ARR_FLT
     MODULE PROCEDURE NC_CONNECT_AVAR_CUB_FLT
     MODULE PROCEDURE NC_CONNECT_AVAR_FDA_FLT

     MODULE PROCEDURE NC_CONNECT_AVAR_SCL_INT
     MODULE PROCEDURE NC_CONNECT_AVAR_VEC_INT
     MODULE PROCEDURE NC_CONNECT_AVAR_ARR_INT
     MODULE PROCEDURE NC_CONNECT_AVAR_CUB_INT
     MODULE PROCEDURE NC_CONNECT_AVAR_FDA_INT
  END INTERFACE

  INTERFACE NC_CONNECT_PVAR
     MODULE PROCEDURE NC_CONNECT_PVAR_SCL_CHR
     MODULE PROCEDURE NC_CONNECT_PVAR_VEC_CHR

     MODULE PROCEDURE NC_CONNECT_PVAR_SCL_DBL
     MODULE PROCEDURE NC_CONNECT_PVAR_VEC_DBL
     MODULE PROCEDURE NC_CONNECT_PVAR_ARR_DBL
     MODULE PROCEDURE NC_CONNECT_PVAR_CUB_DBL
     MODULE PROCEDURE NC_CONNECT_PVAR_FDA_DBL

     MODULE PROCEDURE NC_CONNECT_PVAR_SCL_FLT
     MODULE PROCEDURE NC_CONNECT_PVAR_VEC_FLT
     MODULE PROCEDURE NC_CONNECT_PVAR_ARR_FLT
     MODULE PROCEDURE NC_CONNECT_PVAR_CUB_FLT
     MODULE PROCEDURE NC_CONNECT_PVAR_FDA_FLT

     MODULE PROCEDURE NC_CONNECT_PVAR_SCL_INT
     MODULE PROCEDURE NC_CONNECT_PVAR_VEC_INT
     MODULE PROCEDURE NC_CONNECT_PVAR_ARR_INT
     MODULE PROCEDURE NC_CONNECT_PVAR_CUB_INT
     MODULE PROCEDURE NC_CONNECT_PVAR_FDA_INT
  END INTERFACE

  INTERFACE NC_POINT_VAR
     MODULE PROCEDURE NC_POINT_VAR_SCL_CHR
     MODULE PROCEDURE NC_POINT_VAR_VEC_CHR

     MODULE PROCEDURE NC_POINT_VAR_SCL_DBL
     MODULE PROCEDURE NC_POINT_VAR_VEC_DBL
     MODULE PROCEDURE NC_POINT_VAR_ARR_DBL
     MODULE PROCEDURE NC_POINT_VAR_CUB_DBL
     MODULE PROCEDURE NC_POINT_VAR_FDA_DBL

     MODULE PROCEDURE NC_POINT_VAR_SCL_FLT
     MODULE PROCEDURE NC_POINT_VAR_VEC_FLT
     MODULE PROCEDURE NC_POINT_VAR_ARR_FLT
     MODULE PROCEDURE NC_POINT_VAR_CUB_FLT
     MODULE PROCEDURE NC_POINT_VAR_FDA_FLT

     MODULE PROCEDURE NC_POINT_VAR_SCL_INT
     MODULE PROCEDURE NC_POINT_VAR_VEC_INT
     MODULE PROCEDURE NC_POINT_VAR_ARR_INT
     MODULE PROCEDURE NC_POINT_VAR_CUB_INT
     MODULE PROCEDURE NC_POINT_VAR_FDA_INT
  END INTERFACE


  INTERFACE GET_FILE_TIME
     MODULE PROCEDURE GET_FILE_TIME_NCF
     MODULE PROCEDURE GET_FILE_TIME_NCFTIME
  END INTERFACE
  
  ! DECLARE A HEAD FILE FOR INIT_NCF SO WE CAN KEEP TRACK OF OPEN FILES
  TYPE(NCFILELIST), POINTER, SAVE :: FILEHEAD

  
CONTAINS
!====================================================================
!====================================================================
SUBROUTINE handle_ncerr(status,programer_msg)
  IMPLICIT NONE
  INTEGER,            intent(in) :: status
  CHARACTER(len=*),   intent(in) :: programer_msg
  CHARACTER(len=150)               :: msg
  LOGICAL :: ERROR=.false.
  if(status /=nf90_noerr)&
       & CALL FATAL_ERROR(trim(programer_msg),"NF90ERROR::"//trim(nf90_strerror(status)))
END SUBROUTINE handle_ncerr
!====================================================================
!====================================================================
SUBROUTINE NC_INIT(NCF,NAME)
  IMPLICIT NONE
  TYPE(NCFILE),       POINTER   :: NCF
  CHARACTER(LEN=*), INTENT(IN):: NAME
  LOGICAL FOUND

  IF (.NOT. ASSOCIATED(FILEHEAD))THEN
     FILEHEAD => NEW_FILEHEAD()
  END IF

  NCF => FIND_FILE(FILEHEAD,TRIM(NAME),FOUND)
  IF (FOUND) RETURN

  NCF => NEW_FILE()
  NCF%FNAME=trim(NAME)

END SUBROUTINE NC_INIT
!====================================================================
!====================================================================
SUBROUTINE NC_OPEN(NCF)
# if defined (PIO)
  USE CONTROL, ONLY : USE_MPI_IO_MODE
# endif  
  IMPLICIT NONE
  TYPE(NCFILE), intent(inout) :: NCF
  CHARACTER(LEN=120)          :: errmsg
  integer :: status
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "STARTING NC_OPEN"

  if (NCF%WRITABLE) then
     errmsg="File named: "//NCF%fname//"; Can not be opend by nf90_open"
#if defined (PIO)
     IF(NPROCS_FVCOM > 1 .and. USE_MPI_IO_MODE) THEN
        status = nf90_open(trim(NCF%fname),IOR(NF90_WRITE,NF90_MPIIO),NCF%ncid, comm = MPI_FVCOM_GROUP,info = MPI_INFO_NULL)
     ELSE
        status = nf90_open(trim(NCF%fname),NF90_WRITE,NCF%ncid)
     ENDIF
#else
     status = nf90_open(trim(NCF%fname), NF90_WRITE, NCF%ncid)
#endif
     CALL HANDLE_NCERR(status,trim(errmsg))
     
  else  ! default is open with nf90_nowrite
     errmsg="File: "//TRIM(NCF%fname)//"; Can not be opend by nf90_open"
     status = nf90_open(TRIM(NCF%fname), nf90_noWrite, NCF%ncid)
     CALL HANDLE_NCERR(status,TRIM(errmsg))
  end if

  NCF%OPEN = .TRUE.

  if(DBG_SET(dbg_io))  &
         & write(IPT,*) "Opened File: ",trim(NCF%FNAME)

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_OPEN"
END SUBROUTINE NC_OPEN
!====================================================================
!====================================================================
SUBROUTINE NC_CREATE(NCF)
# if defined (PIO)
  USE CONTROL, ONLY : USE_MPI_IO_MODE
# endif  
  IMPLICIT NONE
  TYPE(NCFILE), intent(inout) :: NCF
  CHARACTER(LEN=120)          :: errmsg,libnetcdf
  integer :: status
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "STARTING NC_CREATE"

  libnetcdf = trim(nf90_inq_libvers())
  errmsg="File named: "//NCF%fname//"; Can not be created by nf90_create"
# if defined (PIO)
# if defined (USE_NETCDF4)
  if(libnetcdf(1:1) /= '4') CALL FATAL_ERROR("The version of NetCDF is old.",&
               & "Please use NetCDF4.")
  IF(NPROCS_FVCOM > 1 .and. USE_MPI_IO_MODE) THEN
    status = nf90_create(trim(NCF%fname),IOR(NF90_NETCDF4,NF90_MPIIO), &
                         NCF%ncid,comm=MPI_FVCOM_GROUP,info=MPI_INFO_NULL)
  ELSE
    status = nf90_create(trim(NCF%fname),NF90_NETCDF4,NCF%ncid)
  ENDIF
# else
!  if(libnetcdf(1:1) == '4') CALL FATAL_ERROR("The version of NetCDF is 4.",&
!               & "Please use NetCDF3.")
  IF(NPROCS_FVCOM > 1 .and. USE_MPI_IO_MODE) THEN
    status = nf90_create(trim(NCF%fname),IOR(NF90_64BIT_OFFSET,NF90_MPIIO), &
                         NCF%ncid,comm=MPI_FVCOM_GROUP,info=MPI_INFO_NULL)
  ELSE
    status = nf90_create(trim(NCF%fname),NF90_64BIT_OFFSET,NCF%ncid)
  ENDIF
# endif
# else
# if defined (USE_NETCDF4)
  if(libnetcdf(1:1) /= '4') CALL FATAL_ERROR("The version of NetCDF is old.",&
               & "Please use NetCDF4.")
  status = nf90_create(trim(NCF%fname),NF90_NETCDF4,NCF%ncid)
# else
!  if(libnetcdf(1:1) == '4') CALL FATAL_ERROR("The version of NetCDF is 4.",&
!               & "Please use NetCDF3.")
  status = nf90_create(trim(NCF%fname),NF90_64BIT_OFFSET,NCF%ncid)
# endif
# endif
  if(status/= nf90_eexist)then
     CALL HANDLE_NCERR(status,trim(errmsg))
  else
     CALL FATAL_ERROR("The file: "//trim(NCF%fname)//"; already exists",&
          & "FVCOM will not overwrite old output files. You must move&
          & or delete them first")
  end if
  NCF%writable = .true.
  NCF%OPEN = .TRUE.
  NCF%INDEFMODE = .TRUE.
  NCF%CONNECTED=.FALSE.

  IF (ASSOCIATED(NCF%FTIME) ) THEN
     NCF%FTIME%NEXT_STKCNT = 0
     NCF%FTIME%STK_LEN=0
  END IF


  if(DBG_SET(dbg_io))  &
       & write(IPT,*) "Created File: ",trim(NCF%FNAME)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_CREATE"
END SUBROUTINE NC_CREATE
!====================================================================
!====================================================================
SUBROUTINE NC_CLOSE(NCF)
  IMPLICIT NONE
  TYPE(NCFILE), intent(INOUT)    :: NCF
  CHARACTER(LEN=120)             :: errmsg
  integer :: status
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "STARTING NC_CLOSE"
  errmsg="File:"//TRIM(NCF%fname)//"; Not open or Already closed file"
  status = nf90_close(NCF%ncid)
  CALL handle_ncerr(status,trim(errmsg))
  NCF%ncid=-1
  NCF%OPEN=.FALSE.
  NCF%INDEFMODE = .FALSE.
  if(DBG_SET(dbg_io))  &
         & write(IPT,*) "Closed File: ",trim(NCF%FNAME)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_CLOSE"
END SUBROUTINE NC_CLOSE
!====================================================================
!====================================================================
SUBROUTINE NC_REFRESH(NCF)
  use control
  IMPLICIT NONE
  TYPE(NCFILE), POINTER :: NCF
  TYPE(NCDIM), pointer  :: DIM
  integer :: status
  LOGICAL :: FOUND
  CHARACTER(LEN=120)     :: errmsg

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "STARTING NC_REFRESH"

  IF(.NOT.ASSOCIATED(NCF)) CALL FATAL_ERROR&
       ("NC_REFRESH: NCF NOT ASSICATED!")

  IF(dbg_set(dbg_io)) write(ipt,*) "====== REFRESHING FILE NAME: "//TRIM(NCF%FNAME)
  

  IF(NCF%OPEN) CALL NC_CLOSE(NCF)
  CALL NC_OPEN(NCF)
  DIM => FIND_UNLIMITED(NCF,FOUND)
  IF(.not. FOUND) RETURN

  status = nf90_inquire_dimension(NCF%NCID,DIM%DIMID,DIM%DIMNAME, DIM%DIM)
  errmsg="Can not get dimensions: "//trim(NCF%FNAME)
  call handle_ncerr(status,errmsg)

  IF(ASSOCIATED(NCF%FTIME)) THEN
     NCF%FTIME%STK_LEN = DIM%DIM
  ELSE
     CALL FATAL_ERROR("NC_REFRESH: FTIME NOT ASSOCIATED FOR FILE:"//&
          & TRIM(NCF%FNAME))
  END IF


  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_REFRESH"

END SUBROUTINE NC_REFRESH
!====================================================================
!====================================================================
SUBROUTINE NC_LOAD(NCF)
  use control
  IMPLICIT NONE
  TYPE(NCFILE), POINTER :: NCF
  TYPE(NCVAR), pointer  :: VAR
  TYPE(NCATT), pointer  :: ATT
  TYPE(NCDIM), pointer  :: DIM
  CHARACTER(LEN=120)             :: errmsg
  integer :: status,i,j, len, nvars, ndims, natts, unlimDimid
  integer, dimension(NF90_MAX_VAR_DIMS):: dimids
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "STARTING NC_LOAD"

  IF(dbg_set(dbg_io)) write(ipt,*) "====== LOADING FILE NAME: "//TRIM(NCF%FNAME)

  IF(NCF%CONNECTED) CALL FATAL_ERROR("CAN NOT LOAD A FILE WHEN ALREADY&
       & CONNECTED!", "FILE NAME: "//TRIM(NCF%FNAME))

! BASIC FILE INFO
  status = nf90_inquire(NCF%NCID, ndims, nVars, nAtts, UnlimDimid)
  errmsg="Can not get file contents: "//trim(NCF%FNAME)
  call handle_ncerr(status,errmsg)
  

! DIMENSIONS
  if(dbg_set(dbg_io)) write(ipt,*) "====== READING FILE DIMENSIONS:"
  do i=1,nDims
     DIM => NC_GET_DIM(NCF,i)
     IF(DIM%DIMID .EQ. UNLIMDIMID) DIM%UNLIMITED = .TRUE.
     if(dbg_set(dbg_io)) write(ipt,*) "      "//trim(DIM%DIMNAME)
     NCF => ADD(NCF,DIM)
  end do

  if(UNLIMDIMID .NE. NCF%UNLIMDIMID) then
     CALL PRINT_FILE(NCF)
     CALL PRINT_DIM_LIST(NCF)
     CALL FATAL_ERROR &
       &("NC_LOAD: UNLIMITED DIMENSION ID FROM nf90_inquire does not m&
       &atch the file objects UNLIMDIMID?")
  end if

  if(nDims /= count_dim_list(NCF) ) then
     if(dbg_set(dbg_log)) call print_dim_list(NCF)
     call fatal_error("The number of dimensions in the file does not m&
          &atch the number loaded in the file object.")
  end if

! call print_dim_list(NCF)

! ATTRIBUTES
  if(dbg_set(dbg_io)) write(ipt,*) "====== READING FILE ATTRIBUTES:"
  do i=1,nAtts
     ATT => NC_GET_ATT(NCF,i)
     if(dbg_set(dbg_io)) write(ipt,*) "      "//trim(ATT%ATTNAME)
     NCF => ADD(NCF,ATT)
  end do

  if(nAtts /= count_att_list(NCF)) then
     if(dbg_set(dbg_log)) call print_att_list(NCF)
     call fatal_error("The number of attributes in the file does not m&
          &atch the number loaded in the file object.")
  end if

! call print_att_list(NCF)


! VARIABLES
  if(dbg_set(dbg_io)) write(ipt,*) "====== READING FILE VARIABLES:"
  do i=1,nVars
     VAR => NC_GET_VAR(NCF,i)
!      call print_dim_list(VAR)
!      call print_att_list(VAR)
     NCF => ADD(NCF,VAR)
  end do

  if(nVars /= count_var_list(NCF)) then
     if(dbg_set(dbg_log)) call print_var_list(NCF)
     call fatal_error("The number of variables in the file does not m&
          &atch the number loaded in the file object.")
  end if

! call print_var_list(NCF)


  STATUS = SET_FILE_TIME_TYPE(NCF)
  ! SEE 'set_file_time_type' for possible results
  !  Since this is a generic routine, take no action based on status...

  NCF%CONNECTED = .TRUE.

 if(DBG_SET(dbg_io)) write(ipt,*) "====== FINISHED LOADING FILE NAME: "//TRIM(NCF%FNAME)

 if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_LOAD"
END SUBROUTINE NC_LOAD
!====================================================================
!====================================================================
SUBROUTINE NC_SAVE(NCF)
  IMPLICIT NONE
  TYPE(NCFILE)         :: NCF
  TYPE(NCVAR), pointer :: VAR
  TYPE(NCATT), pointer :: ATT
  TYPE(NCDIM), pointer :: DIM
  integer              :: attid, dimid, varid
  CHARACTER(LEN=120)             :: errmsg
  integer :: status,i,j, len
  LOGICAL :: FOUND

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "STARTING NC_SAVE"

  IF(NCF%CONNECTED) CALL FATAL_ERROR("CAN NOT SAVE A FILE WHEN ALREADY&
       & CONNECTED!", "FILE NAME: "//TRIM(NCF%FNAME))

  IF(.not. NCF%WRITABLE) CALL FATAL_ERROR("CAN NOT SAVE A FILE WHEN ALREADY&
       & CONNECTED!", "FILE NAME: "//TRIM(NCF%FNAME))

  IF (.NOT. NCF%INDEFMODE) THEN ! FILE MAY COME FROM CREATE OR OPEN...
     status = nf90_redef(NCF%ncid)
     if(status /= NF90_NOERR) &
          & CALL FATAL_ERROR("EXPECTED FILE: "//trim(NCF%FNAME)//"; to be available &
          &for REDEF from an open statement.")
     NCF%INDEFMODE = .TRUE.
  END IF


! DIMENSIONS
  if(dbg_set(dbg_io)) write(ipt,*) "====== WRITING FILE DIMENSIONS:"
  do i=1,count_dim_list(NCF)
     CALL NC_DEF_DIM(NCF,i)
  end do


  if(dbg_set(dbg_io)) write(ipt,*) "====== WRITING FILE GLOBAL ATTRIBUTES:"
  Do i = 1,COUNT_ATT_LIST(NCF)
     ATT => FIND_ATT(NCF,i,FOUND)
     IF (.NOT. FOUND) THEN
        if (DBG_SET(dbg_log)) call print_att_list(NCF)
        CALL FATAL_ERROR&
             &("NC_SAVE: COULD NOT FIND THE GLOBAL ATTRIBUTE WITH CORRECT ATTID W&
             &HILE PUTTING THE ATTRIBUTE IN THE FILE")
     END IF
     
     CALL WRITE_ATT_TYPE(NCF%NCID,NF90_GLOBAL,ATT)
  End Do

  
  if(dbg_set(dbg_io)) write(ipt,*) "====== DEFINE FILE VARIABLES:"
  do i=1,COUNT_VAR_LIST(NCF)
      CALL NC_DEF_VAR(NCF,i)
  end do

  status = NF90_ENDDEF(NCF%NCID)
  errmsg="Can not ENDDEF MODE for file: "//trim(NCF%FNAME)
  call handle_ncerr(status,errmsg)

  NCF%INDEFMODE = .FALSE.
  NCF%CONNECTED = .TRUE.

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_SAVE"
END SUBROUTINE NC_SAVE
!====================================================================
!====================================================================
FUNCTION NC_GET_DIM(NCF,dimid) RESULT(DIM)
  implicit none
  TYPE(NCFILE), INTENT(IN) :: NCF
  TYPE(NCDIM),  POINTER    :: DIM
  INTEGER,      intent(in) :: dimid
  integer :: status
  CHARACTER(LEN=120)     :: errmsg

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_GET_DIM"
  nullify(dim)

  DIM => NEW_DIM()
  DIM%DIMID=DIMID
  status = nf90_inquire_dimension(NCF%NCID,DIMID,DIM%DIMNAME, DIM%DIM)
  errmsg="Can not get dimensions: "//trim(NCF%FNAME)
  call handle_ncerr(status,errmsg)
  
!  if(dimid == NCF%UNLIMDIMID) DIM%UNLIMITED = .TRUE.

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_GET_DIM"

END FUNCTION NC_GET_DIM
!====================================================================
!====================================================================
FUNCTION NC_MAKE_DIM(NAME,LEN) RESULT(DIM)
  IMPLICIT NONE
  TYPE(NCDIM), POINTER :: DIM
  INTEGER, INTENT(IN)  :: LEN
!  logical, intent(in) :: UNLIMITED
  character(len=*), intent(in) :: name
  DIM => NEW_DIM()
  DIM%dimid=-1
  DIM%dimname=NAME
  DIM%dim=LEN
  IF(LEN == NF90_UNLIMITED) DIM%UNLIMITED=.true.
END FUNCTION NC_MAKE_DIM
!====================================================================
!====================================================================
FUNCTION NC_MAKE_RUNTIME_DIM(NAME,LEN) RESULT(DIM)
  USE CONTROL
  IMPLICIT NONE
  TYPE(NCDIM), POINTER :: DIM
  INTEGER, INTENT(IN)  :: LEN
!  logical, intent(in) :: UNLIMITED
  INTEGER, PARAMETER :: TAG = 40003
  INTEGER :: SOURCE, DEST, IERR
# if defined(MULTIPROCESSOR)
  INTEGER :: STAT(MPI_STATUS_SIZE)
# endif

  character(len=*), intent(in) :: name
  DIM => NEW_DIM()
  DIM%dimid=-1
  DIM%dimname=NAME

  IF (.NOT. IOPROC) THEN
     DIM%dim=LEN
     IF(LEN == NF90_UNLIMITED) DIM%UNLIMITED=.true.
  END IF

  IF (USE_MPI_IO_MODE) THEN
# if defined(MULTIPROCESSOR)
!!$     IF(MSR) THEN
     IF(MYID_iogroup == 1) THEN

        DEST = IOPROCID - 1
        CALL MPI_SEND(len,1,MPI_INTEGER,DEST,TAG,MPI_IO_GROUP,IERR)

     ELSE IF (IOPROC) THEN
        
!!$        SOURCE = MSRID -1
        SOURCE = 0
        CALL MPI_RECV(DIM%DIM,1,MPI_INTEGER,SOURCE,TAG,MPI_IO_GROUP,STAT,IERR)

     END IF
# endif
  END IF

END FUNCTION NC_MAKE_RUNTIME_DIM
!====================================================================
!====================================================================
SUBROUTINE NC_DEF_DIM(NCF,DIMID)
  implicit none
  INTEGER,      INTENT(IN)   :: DIMID
  TYPE(NCFILE), INTENT(INOUT):: NCF ! MUST BE ALLOWED TO SET UNLIMDIMID
  TYPE(NCDIM),  POINTER      :: DIM
  INTEGER                    :: status, tmp
  LOGICAL FOUND
  
  DIM => FIND_DIM(NCF,DIMID,FOUND)

  IF (.NOT. FOUND) THEN
     if (DBG_SET(dbg_log)) call print_dim_list(NCF)
     CALL FATAL_ERROR&
          &("NC_DEF_DIM: COULD NOT FIND THE FILE DIMENSION WITH CORRECT DIMID W&
          &HILE DEFINING THE DIMENSION IN THE FILE")
  END IF

  IF(DIM%UNLIMITED) THEN
     status =  nf90_def_dim(NCF%ncid,DIM%dimname, NF90_UNLIMITED,tmp)
  ELSE
     status =  nf90_def_dim(NCF%ncid,DIM%dimname, DIM%dim, tmp)
  END IF

  CALL HANDLE_NCERR(status,"ERROR DURING DEF_DIM, DIMNAME:"//TRIM(DIM%DIMNAME))

  if (tmp .NE. Dim%dimid) CALL FATAL_ERROR &
       &("NC_DEF_DIM: NF90_DEF_DIM returned a dimension id which",&
       & "is different from that set in the dimension object!",&
          & trim(NCF%FNAME)//" : "//TRIM(DIM%DIMNAME))

END SUBROUTINE NC_DEF_DIM
!====================================================================
!====================================================================
FUNCTION NC_MAKE_ATT_CHR_SCL(NAME,VALUES) RESULT(ATT)
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: VALUES

  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_CHAR
  ATT%LEN     = len_trim(VALUES)
  ALLOCATE(ATT%CHR(1))
  ATT%CHR(1)     = VALUES

END FUNCTION NC_MAKE_ATT_CHR_SCL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_ATT_CHR_VEC(NAME,VALUES) RESULT(ATT)
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  character(len=*),ALLOCATABLE, intent(in) :: VALUES(:)

  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_CHAR
  ATT%LEN     = -1
  ALLOCATE(ATT%CHR(size(values)))
  ATT%CHR     = VALUES
END FUNCTION NC_MAKE_ATT_CHR_VEC
!====================================================================
!====================================================================
FUNCTION NC_MAKE_RUNTIME_ATT_CHR_SCL(NAME,VALUES) RESULT(ATT)
  ! NAME IS FIXED, BUT VALUES ARE SENT BY MPI IF NEEDED
  USE CONTROL
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: VALUES
  integer :: ierr, len, DEST,SOURCE
  integer, parameter :: tag = 40004
# if defined(MULTIPROCESSOR)
  INTEGER :: STAT(MPI_STATUS_SIZE)
# endif

  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_CHAR

  IF (.NOT. IOPROC) THEN
     ALLOCATE(ATT%CHR(1))
     ATT%CHR = TRIM(VALUES)
     ATT%LEN     = LEN_TRIM(VALUES)
  END IF

  IF (USE_MPI_IO_MODE) THEN
# if defined(MULTIPROCESSOR)
!!$     IF(MSR) THEN
     IF(MYID_iogroup == 1) THEN

        len = len_trim(VALUES)
        DEST = IOPROCID - 1
        CALL MPI_SEND(len,1,MPI_INTEGER,DEST,TAG,MPI_IO_GROUP,IERR)
        
        CALL MPI_SEND(values(1:len),len,MPI_CHARACTER,DEST,TAG,MPI_IO_GROUP,IERR)

     ELSE IF (IOPROC) THEN
        
        ALLOCATE(ATT%CHR(1))
!!$        SOURCE = MSRID -1
        SOURCE = 0
        CALL MPI_RECV(len,1,MPI_INTEGER,SOURCE,TAG,MPI_IO_GROUP,STAT,IERR)
        ATT%LEN     = len

        CALL MPI_RECV(ATT%CHR(1),len,MPI_CHARACTER,SOURCE,TAG,MPI_IO_GROUP,STAT,IERR)
        
        ATT%CHR(1) = ATT%CHR(1)(1:len)

     END IF
# endif
  END IF

  IF(ATT%len==0) CALL KILL_ATT(ATT)

  
END FUNCTION NC_MAKE_RUNTIME_ATT_CHR_SCL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_RUNTIME_ATT_CHR_VEC(NAME,VALUES) RESULT(ATT)
  ! NAME IS FIXED, BUT VALUES ARE SENT BY MPI IF NEEDED
  USE CONTROL
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  character(len=*),ALLOCATABLE, intent(in) :: VALUES(:)
  integer :: ierr, len, DEST,SOURCE,I,csize
  integer, parameter :: tag = 40004
# if defined(MULTIPROCESSOR)
  INTEGER :: STAT(MPI_STATUS_SIZE)
# endif

!  if(len_trim(VALUES) .GT. Char_max_attlen) & 
!       & Call fatal_error("Can not make attribute: "//trim(NAME),&
!       &"attribute string is too long. Increase 'char_max_attlen' in mod_ncll.F")
  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_CHAR

  IF (.NOT. IOPROC .and. allocated(values)) THEN
     ALLOCATE(ATT%CHR(SIZE(VALUES)))
     ATT%CHR = VALUES
     ATT%LEN     = -1
  END IF

  IF (USE_MPI_IO_MODE) THEN
# if defined(MULTIPROCESSOR)
!!$     IF(MSR) THEN
     IF(MYID_iogroup == 1) THEN
        
        DEST = IOPROCID - 1
        csize=0
        if (allocated(Values)) csize = SIZE(VALUES)
        CALL MPI_SEND(csize,1,MPI_INTEGER,DEST,TAG,MPI_IO_GROUP,IERR)

        DO I = 1,csize
           len = len_trim(VALUES(I))
           CALL MPI_SEND(len,1,MPI_INTEGER,DEST,TAG,MPI_IO_GROUP,IERR)
           
           CALL MPI_SEND(values(I)(1:len),len,MPI_CHARACTER,DEST,TAG,MPI_IO_GROUP,IERR)
        END DO

     ELSE IF (IOPROC) THEN
        
        
!!$        SOURCE = MSRID -1
        SOURCE = 0
        CALL MPI_RECV(csize,1,MPI_INTEGER,SOURCE,TAG,MPI_IO_GROUP,STAT,IERR)
        ATT%LEN     = -1

        if(csize >0) THEN
           
           ALLOCATE(ATT%CHR(csize))
           
           DO I = 1,csize
              CALL MPI_RECV(len,1,MPI_INTEGER,SOURCE,TAG,MPI_IO_GROUP,STAT,IERR)
              
              CALL MPI_RECV(ATT%CHR(I),len,MPI_CHARACTER,SOURCE,TAG,MPI_IO_GROUP,STAT,IERR)
              
              ATT%CHR(I)  = ATT%CHR(I)(1:len)
              
           END DO
        END if

     END IF
# endif
  END IF

  IF(.not.allocated(att%chr)) call Kill_att(ATT)

  
END FUNCTION NC_MAKE_RUNTIME_ATT_CHR_VEC
!====================================================================
!====================================================================
FUNCTION NC_MAKE_ATT_INT_VEC(NAME,values) RESULT(ATT)
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  INTEGER, allocatable, intent(in) :: values(:)

  if(.not. allocated(values)) & 
       & Call fatal_error("Can not make attribute: "//trim(NAME),&
       &"argument 'values' passed must be allocated and contain data")

  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_INT
  ATT%LEN     = size(values)
  allocate(att%int(att%len))
  ATT%int     = values

END FUNCTION NC_MAKE_ATT_INT_VEC
!====================================================================
!====================================================================
FUNCTION NC_MAKE_ATT_INT(NAME,values) RESULT(ATT)
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  INTEGER, intent(in) :: values

  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_INT
  ATT%LEN     = 1
  allocate(att%int(att%len))
  ATT%int     = values

END FUNCTION NC_MAKE_ATT_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_ATT_FLT(NAME,values) RESULT(ATT)
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  REAL(SPA), intent(in) :: values

  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_FLOAT
  ATT%LEN     = 1
  allocate(att%flt(att%len))
  ATT%flt     = values
  
END FUNCTION NC_MAKE_ATT_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_ATT_FLT_VEC(NAME,values) RESULT(ATT)
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  REAL(SPA), allocatable, intent(in) :: values(:)

  if(.not. allocated(values)) & 
       & Call fatal_error("Can not make attribute: "//trim(NAME),&
       &"argument 'values' passed must be allocated and contain data")

  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_FLOAT
  ATT%LEN     = size(values)
  allocate(att%flt(att%len))
  ATT%flt     = values
  
END FUNCTION NC_MAKE_ATT_FLT_VEC
!====================================================================
!====================================================================
FUNCTION NC_MAKE_ATT_DBL(NAME,values) RESULT(ATT)
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  REAL(DP), intent(in) :: values

  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_DOUBLE
  ATT%LEN     = 1
  allocate(att%dbl(att%len))
  ATT%dbl     = values
  
END FUNCTION NC_MAKE_ATT_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_ATT_DBL_VEC(NAME,values) RESULT(ATT)
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  character(len=*), intent(in) :: name
  REAL(DP), allocatable, intent(in) :: values(:)

  if(.not. allocated(values)) & 
       & Call fatal_error("Can not make attribute: "//trim(NAME),&
       &"argument 'values' passed must be allocated and contain data")

  ATT => NEW_ATT()

  ATT%attid   = -1
  ATT%ATTname = TRIM(NAME)
  ATT%xtype   = NF90_DOUBLE
  ATT%LEN     = size(values)
  allocate(att%dbl(att%len))
  ATT%dbl     = values
  
END FUNCTION NC_MAKE_ATT_DBL_VEC
!====================================================================
!====================================================================
FUNCTION NC_GET_GATT(NCF,attid) RESULT(ATT)
  implicit none
  TYPE(NCFILE), INTENT(IN) :: NCF
  TYPE(NCATT),  pointer    :: ATT
  integer,      intent(in) :: attid
  integer                  :: status
  CHARACTER(LEN=120)       :: errmsg
  CHARACTER(LEN=NF90_MAX_NAME+1)   :: NAME

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_GET_GATT"

  status=nf90_inq_attname(NCF%NCID,NF90_GLOBAL,ATTID,NAME)
  errmsg="Can not get a file's global attribute name: "//trim(NCF%FNAME)
  call handle_ncerr(status,errmsg)

  ATT => NEW_ATT()
  ATT%attname=trim(name)
  ATT%attID=attid
  
  status = nf90_inquire_attribute &
       & (NCF%NCID,NF90_GLOBAL,trim(ATT%ATTNAME),ATT%XTYPE,ATT%LEN)

  call READ_ATT_TYPE(NCF%NCID,NF90_GLOBAL,ATT)

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_GET_GATT"
END FUNCTION NC_GET_GATT
!====================================================================
!====================================================================
FUNCTION NC_GET_VATT(VAR,attid) RESULT(ATT)
  implicit none
  TYPE(NCVAR),  INTENT(IN) :: VAR
  TYPE(NCATT),  pointer    :: ATT
  integer,      intent(in) :: attid
  integer                  :: status
  CHARACTER(LEN=120)       :: errmsg
  CHARACTER(LEN=NF90_MAX_NAME+1) :: NAME

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_GET_VATT"

  status=nf90_inq_attname(VAR%NCID,VAR%VARID,ATTID,NAME)
  errmsg="Can not get variable attribute name: "//trim(VAR%VARNAME)
  call handle_ncerr(status,errmsg)

  ATT => NEW_ATT()
  ATT%attname=trim(name)
  ATT%attID=attid

  status = nf90_inquire_attribute &
       & (VAR%NCID,VAR%VARID,trim(ATT%ATTNAME),ATT%XTYPE,ATT%LEN)

  call READ_ATT_TYPE(VAR%NCID,VAR%VARID,ATT)

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_GET_VATT"
END FUNCTION NC_GET_VATT
!====================================================================
!====================================================================
SUBROUTINE READ_ATT_TYPE(NCID,VARID,ATT)
  implicit none
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  type(ncatt),pointer :: ATT
  integer len, status
  CHARACTER(LEN=120)     :: errmsg
  CHARACTER(LEN=4)     :: clen

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START READ_ATT_TYPE"

  len=ATT%LEN
  write(clen,'(I4.4)') len
  status = 1
  select case(ATT%XTYPE)
  case(NF90_BYTE)
     allocate(ATT%int(len),stat=status)
     if(status/=0) CALL Fatal_error("READ_ATT_TYPE could not allocate integer("//clen//")")
     status = nf90_get_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%int )
     errmsg="Can not get variable attribute (byte):"&
          & //trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_SHORT)
     allocate(ATT%int(len),stat=status)
     if(status/=0) CALL Fatal_error("READ_ATT_TYPE could not allocate integer("//clen//")")
     status = nf90_get_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%int )
     errmsg="Can not get variable attribute&
          & (short):"//trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_INT)
     allocate(ATT%int(len),stat=status)
     if(status/=0) CALL Fatal_error("READ_ATT_TYPE could not allocate integer("//clen//")")
     status = nf90_get_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%int )
     errmsg="Can not get variable attribute (int):"&
          & //trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_FLOAT)
     allocate(ATT%flt(len),stat=status)
     if(status/=0) CALL Fatal_error("READ_ATT_TYPE could not allocate float("//clen//")")
     status = nf90_get_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%flt )
     errmsg="Can not get variable attribute&
          & (float):"//trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_DOUBLE)
     allocate(ATT%DBL(len),stat=status)
     if(status/=0) CALL Fatal_error("READ_ATT_TYPE could not allocate double("//clen//")")
     status = nf90_get_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%dbl )
     errmsg="Can not get variable attribute (double):"&
          & //trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_CHAR)

     CALL CHAR_ATT_READ_HELPER(NCID,VARID,ATT%ATTNAME,ATT%chr,len)
     
  case default
     if(status/=0) CALL Fatal_error("READ_ATT_TYPE hit default case: b&
          &ad att type")
  end select  

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END READ_ATT_TYPE"

END SUBROUTINE READ_ATT_TYPE
!====================================================================
!====================================================================
SUBROUTINE CHAR_ATT_READ_HELPER(NCID,VARID,ATTNAME,CHR,LEN)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: NCID, VARID, LEN
  CHARACTER(LEN=*), INTENT(IN) :: ATTNAME
!  CHARACTER(LEN=LEN+1) :: TEMP
  CHARACTER(LEN=LEN) :: TEMP
  CHARACTER(LEN=*), ALLOCATABLE :: CHR(:)
  CHARACTER(LEN=120)     :: errmsg
  INTEGER :: STATUS

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START CHAR_ATT_READ_HELPER"

  status = nf90_get_att(NCID,VARID,ATTNAME, TEMP )
  errmsg="Can not get variable attribute (char):" &
       & //trim(ATTNAME)
  call handle_ncerr(status,errmsg)

  ! ACHAR(10) is a carage return!
  CALL SPLIT_STRING(TEMP,ACHAR(10), CHR)

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END CHAR_ATT_READ_HELPER"

END SUBROUTINE CHAR_ATT_READ_HELPER
!====================================================================
!====================================================================
SUBROUTINE WRITE_ATT_TYPE(ncid,varid,ATT)
  implicit none
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  TYPE(NCATT), intent(in) :: ATT
  integer len, status, I,slen
  CHARACTER(LEN=120)     :: errmsg

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START WRITE_ATT_TYPE"

  len=ATT%LEN
  status=1
  select case(ATT%XTYPE)
  case(NF90_BYTE)
     status = nf90_put_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%int )
     errmsg="Can not set variable attribute (byte):"&
          & //trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_SHORT)
     status = nf90_put_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%int )
     errmsg="Can not set variable attribute&
          & (short):"//trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_INT)
     status = nf90_put_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%int )
     errmsg="Can not set variable attribute (int):"&
          & //trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_FLOAT)
     status = nf90_put_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%flt )
     errmsg="Can not set variable attribute&
          & (float):"//trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_DOUBLE)
     status = nf90_put_att(NCID,VARID,&
          & ATT%ATTNAME, &
          & ATT%dbl )
     errmsg="Can not set variable attribute (double):"&
          & //trim(ATT%ATTNAME)
     call handle_ncerr(status,errmsg)
     
  case(NF90_CHAR)

     slen = 0
     DO I = 1,size(ATT%chr)
        slen = slen + len_trim(adjustl(ATT%chr(i))) + 1
     END DO

     CALL CHAR_ATT_WRITE_HELPER(NCID,VARID, ATT%ATTNAME,ATT%chr,slen)
     
  case default
     if(status/=0) CALL Fatal_error("WRITE_ATT_TYPE hit default case: b&
          &ad att type")
  end select  

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END WRITE_ATT_TYPE"

END SUBROUTINE WRITE_ATT_TYPE
!====================================================================
!====================================================================
SUBROUTINE CHAR_ATT_WRITE_HELPER(NCID,VARID,ATTNAME,CHR,LEN)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: NCID, VARID, LEN
  CHARACTER(LEN=*), INTENT(IN) :: ATTNAME
  CHARACTER(LEN=LEN) :: TEMP
  CHARACTER(LEN=*), ALLOCATABLE :: CHR(:)
  CHARACTER(LEN=120)     :: errmsg
  INTEGER :: STATUS, I, csize

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START CHAR_ATT_WRITE_HELPER"
  TEMP=""

  csize =size(chr) 
  DO I = 1,size(chr)
     TEMP = TRIM(TEMP)//trim(adjustl(chr(i)))
     
     ! DO not att charage return to last string
     IF(I<Csize) TEMP = TRIM(TEMP)//ACHAR(10)
  END DO
  
  status = nf90_put_att(NCID,VARID,ATTNAME, TRIM(TEMP) )
  errmsg="Can not put variable attribute (char):" &
       & //trim(ATTNAME)
  call handle_ncerr(status,errmsg)

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END CHAR_ATT_WRITE_HELPER"

END SUBROUTINE CHAR_ATT_WRITE_HELPER
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_SCL_CHR(NAME,VALUES,DIM1,DIM2) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM),  POINTER :: DIM1 ! MUST BE STR LENGTH
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM2 ! IF PRESENT THIS MUST BE UNLIMITED
  character(len=*), intent(in) :: name
  character(len=80), target, intent(in) :: values
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_SCL_CHR"

  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
 

  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
    & ("NC_MAKE_AVAR: WHEN MAKING CHARACTER DATA, THE FIRST DIMENSION",&
       & "MUST BE THE STING LENGTH. THE LAST DIMENSION (IF PRESENT) MUST BE UNLIMITED",&
       &"VARIABLE NAME: "//TRIM(NAME))

  VAR => NEW_VAR()
  
  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_CHAR
  VAR%scl_chr => values
  VAR => add(VAR,COPY_DIM(DIM1))
  IF(present(DIM2)) THEN
     if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))

     IF(.NOT. DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: WHEN MAKING CHARACTER DATA, THE LAST DIMENSION",&
          & "MUST BE THE UNLIMITED DIMENSION",&
          & "VARIABLE NAME: "//TRIM(NAME))

      VAR => add(VAR,COPY_DIM(DIM2))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_SCL_CHR"
END FUNCTION NC_MAKE_AVAR_SCL_CHR
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_VEC_CHR(NAME,VALUES,DIM1,DIM2,DIM3) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM),  POINTER :: DIM1 ! MUST BE STR LENGTH
  TYPE(NCDIM),  POINTER :: DIM2 ! Number of strings
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM3 ! IF PRESENT THIS MUST BE UNLIMITED
  character(len=*), intent(in) :: name
  character(len=80), target,allocatable, intent(in) :: values(:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_VEC_CHR"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))

  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))

  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: WHEN MAKING CHARACTER DATA, THE FIRST DIMENSION",&
       & "MUST BE THE STRING LENGTH. THE LAST DIMENSION (IF PRESENT) MUST BE UNLIMITED",&
       &"VARIABLE NAME: "//TRIM(NAME))
  

  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))

!  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
!       & ("NC_MAKE_AVAR: WHEN MAKING CHARACTER DATA, THE FIRST DIMENSION",&
!       & "MUST BE THE STRING LENGTH. THE LAST DIMENSION (IF PRESENT) MUST BE UNLIMITED",&
!       &"VARIABLE NAME: "//TRIM(NAME))

  VAR => NEW_VAR()
  
  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_CHAR
  VAR%vec_chr => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  IF(present(DIM3)) THEN
     if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))

     IF(.NOT. DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: WHEN MAKING CHARACTER DATA, THE LAST DIMENSION",&
          & "MUST BE THE UNLIMITED DIMENSION",&
          & "VARIABLE NAME: "//TRIM(NAME))

      VAR => add(VAR,COPY_DIM(DIM3))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_VEC_CHR"
END FUNCTION NC_MAKE_AVAR_VEC_CHR
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_SCL_DBL(NAME,VALUES,DIM1) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM1
  character(len=*), intent(in) :: name
  REAL(DP), target, intent(in) :: values
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_SCL_DBL"

  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%scl_dbl => values
  IF(present(DIM1)) THEN
     if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM1))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_SCL_DBL"
END FUNCTION NC_MAKE_AVAR_SCL_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_VEC_DBL(NAME,VALUES,DIM1,DIM2) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM2
  character(len=*), intent(in) :: name
  REAL(DP), allocatable, target, intent(in) :: values(:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_VEC_DBL"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%vec_dbl => values
  VAR => add(VAR,COPY_DIM(DIM1))
  if(present(DIM2)) then
     if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM2))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_VEC_DBL"
END FUNCTION NC_MAKE_AVAR_VEC_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_ARR_DBL(NAME,VALUES,DIM1,DIM2,DIM3) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM3
  character(len=*), intent(in) :: name
  REAL(DP), allocatable, target, intent(in) :: values(:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_ARR_DBL"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%arr_dbl => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  if(present(DIM3)) then
     if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM3))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_ARR_DBL"
END FUNCTION NC_MAKE_AVAR_ARR_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_CUB_DBL(NAME,VALUES,DIM1,DIM2,DIM3,DIM4) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM4
  character(len=*), intent(in) :: name
  REAL(DP), allocatable, target, intent(in) :: values(:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_CUB_DBL"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%cub_dbl => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  if(present(DIM4)) then
     if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM4%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM4))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_CUB_DBL"
END FUNCTION NC_MAKE_AVAR_CUB_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_FDA_DBL(NAME,VALUES,DIM1,DIM2,DIM3,DIM4,DIM5) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), POINTER :: DIM4
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM5
  character(len=*), intent(in) :: name
  REAL(DP), allocatable, target, intent(in) :: values(:,:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_FDA_DBL"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM4%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%fda_dbl => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  VAR => add(VAR,COPY_DIM(DIM4))
  if(present(DIM5)) then
     if(.NOT. ASSOCIATED(DIM5)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM5%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM5))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_FDA_DBL"
END FUNCTION NC_MAKE_AVAR_FDA_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_SCL_FLT(NAME,VALUES,DIM1) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM1
  character(len=*), intent(in) :: name
  REAL(SPA), target, intent(in) :: values
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_SCL_FLT"

  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%scl_flt => values
  IF(present(DIM1)) THEN
     if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM1))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_SCL_FLT"
END FUNCTION NC_MAKE_AVAR_SCL_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_VEC_FLT(NAME,VALUES,DIM1,DIM2) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM2
  character(len=*), intent(in) :: name
  REAL(SPA), allocatable, target, intent(in) :: values(:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_VEC_FLT"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%vec_flt => values(1:)
  VAR => add(VAR,COPY_DIM(DIM1))
  if(present(DIM2)) then
     if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM2))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_VEC_FLT"
END FUNCTION NC_MAKE_AVAR_VEC_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_ARR_FLT(NAME,VALUES,DIM1,DIM2,DIM3) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM3
  character(len=*), intent(in) :: name
  REAL(SPA), allocatable, target, intent(in) :: values(:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_ARR_FLT"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%arr_flt => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  if(present(DIM3)) then
     if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM3))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_ARR_FLT"
END FUNCTION NC_MAKE_AVAR_ARR_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_CUB_FLT(NAME,VALUES,DIM1,DIM2,DIM3,DIM4) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM4
  character(len=*), intent(in) :: name
  REAL(SPA), allocatable, target, intent(in) :: values(:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_CUB_FLT"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%cub_flt => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  if(present(DIM4)) then
     if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM4%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM4))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_CUB_FLT"
END FUNCTION NC_MAKE_AVAR_CUB_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_FDA_FLT(NAME,VALUES,DIM1,DIM2,DIM3,DIM4,DIM5) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), POINTER :: DIM4
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM5
  character(len=*), intent(in) :: name
  REAL(SPA), allocatable, target, intent(in) :: values(:,:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_FDA_FLT"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM4%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%fda_flt => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  VAR => add(VAR,COPY_DIM(DIM4))
  if(present(DIM5)) then
     if(.NOT. ASSOCIATED(DIM5)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM5%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM5))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_FDA_FLT"
END FUNCTION NC_MAKE_AVAR_FDA_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_SCL_INT(NAME,VALUES,DIM1) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM1
  character(len=*), intent(in) :: name
  INTEGER, target, intent(in) :: values
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_SCL_INT"

  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%scl_int => values
  IF(present(DIM1)) THEN
     if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM1))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_SCL_INT"
END FUNCTION NC_MAKE_AVAR_SCL_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_VEC_INT(NAME,VALUES,DIM1,DIM2) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM2
  character(len=*), intent(in) :: name
  INTEGER, allocatable, target, intent(in) :: values(:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_VEC_INT"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%vec_int => values
  VAR => add(VAR,COPY_DIM(DIM1))
  if(present(DIM2)) then
     if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM2))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_VEC_INT"
END FUNCTION NC_MAKE_AVAR_VEC_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_ARR_INT(NAME,VALUES,DIM1,DIM2,DIM3) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM3
  character(len=*), intent(in) :: name
  INTEGER, allocatable, target, intent(in) :: values(:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_ARR_INT"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%arr_int => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  if(present(DIM3)) then
     if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM3))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_ARR_INT"
END FUNCTION NC_MAKE_AVAR_ARR_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_CUB_INT(NAME,VALUES,DIM1,DIM2,DIM3,DIM4) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM4
  character(len=*), intent(in) :: name
  INTEGER, allocatable, target, intent(in) :: values(:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_CUB_INT"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%cub_int => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  if(present(DIM4)) then
     if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM4%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM4))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_CUB_INT"
END FUNCTION NC_MAKE_AVAR_CUB_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_AVAR_FDA_INT(NAME,VALUES,DIM1,DIM2,DIM3,DIM4,DIM5) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), POINTER :: DIM4
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM5
  character(len=*), intent(in) :: name
  INTEGER, allocatable, target, intent(in) :: values(:,:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_AVAR_FDA_INT"

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
       & ("NC_MAKE_AVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM4%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_AVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%fda_int => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  VAR => add(VAR,COPY_DIM(DIM4))
  if(present(DIM5)) then
     if(.NOT. ASSOCIATED(DIM5)) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM5%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_AVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM5))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_AVAR_FDA_INT"
END FUNCTION NC_MAKE_AVAR_FDA_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_SCL_CHR(NAME,VALUES,DIM1,DIM2) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM),  POINTER :: DIM1 ! MUST BE STR LENGTH
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM2 ! IF PRESENT THIS MUST BE UNLIMITED
  character(len=*), intent(in) :: name
  character(len=80), POINTER, intent(in) :: values
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_SCL_CHR"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))

  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
 

  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
    & ("NC_MAKE_PVAR: WHEN MAKING CHARACTER DATA, THE FIRST DIMENSION",&
       & "MUST BE THE STING LENGTH. THE LAST DIMENSION (IF PRESENT) MUST BE UNLIMITED",&
       &"VARIABLE NAME: "//TRIM(NAME))

  VAR => NEW_VAR()
  
  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_CHAR
  VAR%scl_chr => values
  VAR => add(VAR,COPY_DIM(DIM1))
  IF(present(DIM2)) THEN
     if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))

     IF(.NOT. DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: WHEN MAKING CHARACTER DATA, THE LAST DIMENSION",&
          & "MUST BE THE UNLIMITED DIMENSION",&
          & "VARIABLE NAME: "//TRIM(NAME))

      VAR => add(VAR,COPY_DIM(DIM2))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_SCL_CHR"
END FUNCTION NC_MAKE_PVAR_SCL_CHR
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_VEC_CHR(NAME,VALUES,DIM1,DIM2,DIM3) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM),  POINTER :: DIM1 ! MUST BE STR LENGTH
  TYPE(NCDIM),  POINTER :: DIM2 ! Number of strings
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM3 ! IF PRESENT THIS MUST BE UNLIMITED
  character(len=*), intent(in) :: name
  character(len=80), POINTER, intent(in) :: values(:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_VEC_CHR"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))

  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))

  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: WHEN MAKING CHARACTER DATA, THE FIRST DIMENSION",&
       & "MUST BE THE STRING LENGTH. THE LAST DIMENSION (IF PRESENT) MUST BE UNLIMITED",&
       &"VARIABLE NAME: "//TRIM(NAME))
  
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))

!  if(DIM2%UNLIMITED) CALL WARNING &
!       & ("NC_MAKE_PVAR: WHEN MAKING CHARACTER DATA, THE FIRST DIMENSION",&
!       & "MUST BE THE STRING LENGTH. THE LAST DIMENSION (IF PRESENT) MUST BE UNLIMITED",&
!       &"VARIABLE NAME: "//TRIM(NAME))

  VAR => NEW_VAR()
  
  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_CHAR
  VAR%vec_chr => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  IF(present(DIM3)) THEN
     if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))

     IF(.NOT. DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: WHEN MAKING CHARACTER DATA, THE LAST DIMENSION",&
          & "MUST BE THE UNLIMITED DIMENSION",&
          & "VARIABLE NAME: "//TRIM(NAME))

      VAR => add(VAR,COPY_DIM(DIM3))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_VEC_CHR"
END FUNCTION NC_MAKE_PVAR_VEC_CHR
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_SCL_DBL(NAME,VALUES,DIM1) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM1
  character(len=*), intent(in) :: name
  REAL(DP), POINTER, intent(in) :: values
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_SCL_DBL"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%scl_dbl => values
  IF(present(DIM1)) THEN
     if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM1))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_SCL_DBL"
END FUNCTION NC_MAKE_PVAR_SCL_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_VEC_DBL(NAME,VALUES,DIM1,DIM2) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM2
  character(len=*), intent(in) :: name
  REAL(DP), POINTER, intent(in) :: values(:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_VEC_DBL"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
!       & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  

  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%vec_dbl => values
  VAR => add(VAR,COPY_DIM(DIM1))
  if(present(DIM2)) then
     if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM2))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_VEC_DBL"
END FUNCTION NC_MAKE_PVAR_VEC_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_ARR_DBL(NAME,VALUES,DIM1,DIM2,DIM3) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM3
  character(len=*), intent(in) :: name
  REAL(DP), POINTER, intent(in) :: values(:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_ARR_DBL"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%arr_dbl => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  if(present(DIM3)) then
     if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM3))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_ARR_DBL"
END FUNCTION NC_MAKE_PVAR_ARR_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_CUB_DBL(NAME,VALUES,DIM1,DIM2,DIM3,DIM4) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM4
  character(len=*), intent(in) :: name
  REAL(DP), POINTER, intent(in) :: values(:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_CUB_DBL"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%cub_dbl => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  if(present(DIM4)) then
     if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM4%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM4))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_CUB_DBL"
END FUNCTION NC_MAKE_PVAR_CUB_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_FDA_DBL(NAME,VALUES,DIM1,DIM2,DIM3,DIM4,DIM5) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), POINTER :: DIM4
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM5
  character(len=*), intent(in) :: name
  REAL(DP), POINTER, intent(in) :: values(:,:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_FDA_DBL"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM4%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_DOUBLE
  VAR%fda_dbl => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  VAR => add(VAR,COPY_DIM(DIM4))
  if(present(DIM5)) then
     if(.NOT. ASSOCIATED(DIM5)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM5%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM5))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_FDA_DBL"
END FUNCTION NC_MAKE_PVAR_FDA_DBL
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_SCL_FLT(NAME,VALUES,DIM1) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM1
  character(len=*), intent(in) :: name
  REAL(SPA), POINTER, intent(in) :: values
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_SCL_FLT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))

  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%scl_flt => values
  IF(present(DIM1)) THEN
     if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM1))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_SCL_FLT"
END FUNCTION NC_MAKE_PVAR_SCL_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_VEC_FLT(NAME,VALUES,DIM1,DIM2) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM2
  character(len=*), intent(in) :: name
  REAL(SPA), POINTER, intent(in) :: values(:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_VEC_FLT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM1%UNLIMITED) CALL WARNING &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%vec_flt => values(1:)
  VAR => add(VAR,COPY_DIM(DIM1))
  if(present(DIM2)) then
     if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM2))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_VEC_FLT"
END FUNCTION NC_MAKE_PVAR_VEC_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_ARR_FLT(NAME,VALUES,DIM1,DIM2,DIM3) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM3
  character(len=*), intent(in) :: name
  REAL(SPA), POINTER, intent(in) :: values(:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_ARR_FLT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%arr_flt => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  if(present(DIM3)) then
     if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM3))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_ARR_FLT"
END FUNCTION NC_MAKE_PVAR_ARR_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_CUB_FLT(NAME,VALUES,DIM1,DIM2,DIM3,DIM4) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM4
  character(len=*), intent(in) :: name
  REAL(SPA), POINTER, intent(in) :: values(:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_CUB_FLT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%cub_flt => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  if(present(DIM4)) then
     if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM4%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM4))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_CUB_FLT"
END FUNCTION NC_MAKE_PVAR_CUB_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_FDA_FLT(NAME,VALUES,DIM1,DIM2,DIM3,DIM4,DIM5) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), POINTER :: DIM4
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM5
  character(len=*), intent(in) :: name
  REAL(SPA), POINTER, intent(in) :: values(:,:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_FDA_FLT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM4%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_FLOAT
  VAR%fda_flt => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  VAR => add(VAR,COPY_DIM(DIM4))
  if(present(DIM5)) then
     if(.NOT. ASSOCIATED(DIM5)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM5%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM5))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_FDA_FLT"
END FUNCTION NC_MAKE_PVAR_FDA_FLT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_SCL_INT(NAME,VALUES,DIM1) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM1
  character(len=*), intent(in) :: name
  INTEGER, POINTER, intent(in) :: values
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_SCL_INT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))

  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%scl_int => values
  IF(present(DIM1)) THEN
     if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM1))
  END IF

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_SCL_INT"
END FUNCTION NC_MAKE_PVAR_SCL_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_VEC_INT(NAME,VALUES,DIM1,DIM2) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM2
  character(len=*), intent(in) :: name
  INTEGER, POINTER, intent(in) :: values(:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_VEC_INT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM1%UNLIMITED) CALL WARNING &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%vec_int => values
  VAR => add(VAR,COPY_DIM(DIM1))
  if(present(DIM2)) then
     if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM2))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_VEC_INT"
END FUNCTION NC_MAKE_PVAR_VEC_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_ARR_INT(NAME,VALUES,DIM1,DIM2,DIM3) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM3
  character(len=*), intent(in) :: name
  INTEGER, POINTER, intent(in) :: values(:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_ARR_INT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%arr_int => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  if(present(DIM3)) then
     if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM3))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_ARR_INT"
END FUNCTION NC_MAKE_PVAR_ARR_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_CUB_INT(NAME,VALUES,DIM1,DIM2,DIM3,DIM4) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM4
  character(len=*), intent(in) :: name
  INTEGER, POINTER, intent(in) :: values(:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_CUB_INT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%cub_int => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  if(present(DIM4)) then
     if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM4%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM4))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_CUB_INT"
END FUNCTION NC_MAKE_PVAR_CUB_INT
!====================================================================
!====================================================================
FUNCTION NC_MAKE_PVAR_FDA_INT(NAME,VALUES,DIM1,DIM2,DIM3,DIM4,DIM5) RESULT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM1
  TYPE(NCDIM), POINTER :: DIM2
  TYPE(NCDIM), POINTER :: DIM3
  TYPE(NCDIM), POINTER :: DIM4
  TYPE(NCDIM), OPTIONAL, POINTER :: DIM5
  character(len=*), intent(in) :: name
  INTEGER, POINTER, intent(in) :: values(:,:,:,:)
  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_MAKE_PVAR_FDA_INT"

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: VALUES ARGUMENT MUST ALREADY BE ASSOCIATED!",&
       & "ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM1)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM1%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM2)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM2%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM3)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(DIM3%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  if(.NOT. ASSOCIATED(DIM4)) CALL FATAL_ERROR &
       & ("NC_MAKE_PVAR: THE DIMENSION ARGUMENT MUST BE ASSOCIATED",&
       &"ERROR MAKING VARIABLE: "//TRIM(NAME))
!  if(DIM4%UNLIMITED) CALL FATAL_ERROR &
!          & ("NC_MAKE_PVAR: NON-OPTIONAL DIMENSION ARGUMENT MUST BE NOT BE UNLIMITED (TIME)",&
!          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
  VAR => NEW_VAR()

  VAR%VARID   = -1
  VAR%VARNAME = TRIM(NAME)
  VAR%xtype   = NF90_INT
  VAR%fda_int => values
  VAR => add(VAR,COPY_DIM(DIM1))
  VAR => add(VAR,COPY_DIM(DIM2))
  VAR => add(VAR,COPY_DIM(DIM3))
  VAR => add(VAR,COPY_DIM(DIM4))
  if(present(DIM5)) then
     if(.NOT. ASSOCIATED(DIM5)) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE ASSOCIATED",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     if(.NOT. DIM5%UNLIMITED) CALL FATAL_ERROR &
          & ("NC_MAKE_PVAR: OPTIONAL DIMENSION ARGUMENT MUST BE UNLIMITED (TIME)",&
          &"ERROR MAKING VARIABLE: "//TRIM(NAME))
     VAR => add(VAR,COPY_DIM(DIM5))
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_MAKE_PVAR_FDA_INT"
END FUNCTION NC_MAKE_PVAR_FDA_INT
!====================================================================
!====================================================================
SUBROUTINE NC_DISCONNECT(VAR)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_DISCONNECT: THE VARIABLE ARGUMENT MUST BE ASSOCIATED")

    nullify(VAR%scl_int)
    nullify(VAR%vec_int)
    nullify(VAR%arr_int)
    nullify(VAR%cub_int)
    nullify(VAR%fda_int)

    nullify(VAR%scl_flt)
    nullify(VAR%vec_flt)
    nullify(VAR%arr_flt)
    nullify(VAR%cub_flt)
    nullify(VAR%fda_flt)

    nullify(VAR%scl_dbl)
    nullify(VAR%vec_dbl)
    nullify(VAR%arr_dbl)
    nullify(VAR%cub_dbl)
    nullify(VAR%fda_dbl)

    nullify(var%scl_chr)
    nullify(var%vec_chr)

END SUBROUTINE NC_DISCONNECT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_SCL_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA), target, intent(in) :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_FLT)")
  
  VAR%SCL_FLT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_SCL_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_VEC_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA), allocatable, target, intent(in) :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_FLT)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  
  VAR%VEC_FLT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_VEC_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_ARR_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA), allocatable, target, intent(in) :: values(:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (ARR_FLT)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%ARR_FLT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_ARR_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_CUB_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA), allocatable, target, intent(in) :: values(:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (CUB_FLT)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%CUB_FLT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_CUB_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_FDA_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA), allocatable, target, intent(in) :: values(:,:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (FDA_FLT)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%FDA_FLT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_FDA_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_SCL_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER, target, intent(in) :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_INT)")
  
  VAR%SCL_INT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_SCL_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_VEC_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER, allocatable, target, intent(in) :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_INT)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%VEC_INT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_VEC_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_ARR_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER, allocatable, target, intent(in) :: values(:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (ARR_INT)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%ARR_INT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_ARR_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_CUB_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER, allocatable, target, intent(in) :: values(:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (CUB_INT)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%CUB_INT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_CUB_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_FDA_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER, allocatable, target, intent(in) :: values(:,:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (FDA_INT)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%FDA_INT => VALUES

END SUBROUTINE NC_CONNECT_AVAR_FDA_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_SCL_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(DP), target, intent(in) :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_DBL)")
  
  VAR%SCL_DBL => VALUES

END SUBROUTINE NC_CONNECT_AVAR_SCL_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_VEC_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  real(DP), allocatable, target, intent(in) :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_DBL)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%VEC_DBL => VALUES

END SUBROUTINE NC_CONNECT_AVAR_VEC_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_ARR_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  real(DP), allocatable, target, intent(in) :: values(:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (ARR_DBL)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%ARR_DBL => VALUES

END SUBROUTINE NC_CONNECT_AVAR_ARR_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_CUB_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  real(DP), allocatable, target, intent(in) :: values(:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (CUB_DBL)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%CUB_DBL => VALUES

END SUBROUTINE NC_CONNECT_AVAR_CUB_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_FDA_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  real(DP), allocatable, target, intent(in) :: values(:,:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (FDA_DBL)")

  IF(.NOT. ALLOCATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%FDA_DBL => VALUES

END SUBROUTINE NC_CONNECT_AVAR_FDA_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_SCL_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA), POINTER :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_FLT)")

  IF(.NOT. ASSOCIATED(values)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_FLT)")
  
  VAR%SCL_FLT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_SCL_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_VEC_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA), POINTER :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_FLT)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  
  VAR%VEC_FLT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_VEC_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_ARR_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA),POINTER  :: values(:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (ARR_FLT)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%ARR_FLT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_ARR_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_CUB_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA),POINTER  :: values(:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (CUB_FLT)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%CUB_FLT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_CUB_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_FDA_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(SPA),POINTER  :: values(:,:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (FDA_FLT)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%FDA_FLT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_FDA_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_SCL_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER, POINTER :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_INT)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))

  VAR%SCL_INT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_SCL_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_VEC_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER,POINTER :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_INT)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%VEC_INT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_VEC_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_ARR_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER, POINTER :: values(:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (ARR_INT)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%ARR_INT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_ARR_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_CUB_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER, POINTER :: values(:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (CUB_INT)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%CUB_INT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_CUB_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_FDA_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  INTEGER, POINTER :: values(:,:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (FDA_INT)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%FDA_INT => VALUES

END SUBROUTINE NC_CONNECT_PVAR_FDA_INT
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_SCL_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  REAL(DP), POINTER :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_DBL)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))

  VAR%SCL_DBL => VALUES

END SUBROUTINE NC_CONNECT_PVAR_SCL_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_VEC_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  real(DP), POINTER :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_DBL)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
    
  VAR%VEC_DBL => VALUES

END SUBROUTINE NC_CONNECT_PVAR_VEC_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_ARR_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  real(DP), POINTER :: values(:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (ARR_DBL)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%ARR_DBL => VALUES

END SUBROUTINE NC_CONNECT_PVAR_ARR_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_CUB_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  real(DP), POINTER :: values(:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (CUB_DBL)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%CUB_DBL => VALUES

END SUBROUTINE NC_CONNECT_PVAR_CUB_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_FDA_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  real(DP), POINTER :: values(:,:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (FDA_DBL)")

  IF(.NOT. ASSOCIATED(VALUES)) CALL FATAL_ERROR &
       & ("NC_CONNECT_PVAR: VALUES ARGUMENT MUST ALREADY BE ALLOCATED!",&
       & "ERROR CONNECTING VARIABLE: "//TRIM(VAR%VARNAME))
  
  VAR%FDA_DBL => VALUES

END SUBROUTINE NC_CONNECT_PVAR_FDA_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_VAR_SCL_CHR(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  CHARACTER(LEN=80), target, intent(in) :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_CHR)")

  ! DO NO TYPE OR DIMENSION MATCHING HERE... TO COMPLICATED
  
  VAR%SCL_CHR => VALUES


END SUBROUTINE NC_CONNECT_VAR_SCL_CHR
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_AVAR_VEC_CHR(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  CHARACTER(LEN=80), target,ALLOCATABLE, intent(in) :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_CHR)")

  IF(.NOT. ALLOCATED(values)) CALL FATAL_ERROR &
       & ("THE CHARACTER ALLOCATABLE MUST BE ALLOCATED BEFORE BEING PASSE&
       &D TO NC_CONNECT_VAR:", VAR%VARNAME)

  ! DO NO TYPE OR DIMENSION MATCHING HERE... TO COMPLICATED
  
  VAR%VEC_CHR => VALUES


END SUBROUTINE NC_CONNECT_AVAR_VEC_CHR
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_SCL_CHR(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  CHARACTER(LEN=80), pointer, intent(in) :: values

  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_CHR)")

  IF(.NOT. ASSOCIATED(values)) CALL FATAL_ERROR &
       & ("THE CHARACTER POINTER MUST BE ASSOCIATED BEFORE BEING PASSE&
       &D TO NC_CONNECT_VAR:", VAR%VARNAME)

  ! DO NO TYPE OR DIMENSION MATCHING HERE... TO COMPLICATED
  
  VAR%SCL_CHR => VALUES


END SUBROUTINE NC_CONNECT_PVAR_SCL_CHR
!====================================================================
!====================================================================
SUBROUTINE NC_CONNECT_PVAR_VEC_CHR(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  CHARACTER(LEN=80), pointer, intent(in) :: values(:)

  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_CONNECT_AVAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_CHR)")

  IF(.NOT. ASSOCIATED(values)) CALL FATAL_ERROR &
       & ("THE CHARACTER POINTER MUST BE ASSOCIATED BEFORE BEING PASSE&
       &D TO NC_CONNECT_VAR:", VAR%VARNAME)

  ! DO NO TYPE OR DIMENSION MATCHING HERE... TO COMPLICATED
  
  VAR%VEC_CHR => VALUES


END SUBROUTINE NC_CONNECT_PVAR_VEC_CHR
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_SCL_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(SPA), POINTER :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_FLT)")

  IF(.NOT. ASSOCIATED(VAR%SCL_FLT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (SCL_FLT)")
  END IF

  VALUES => VAR%SCL_FLT

END SUBROUTINE NC_POINT_VAR_SCL_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_VEC_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(SPA), POINTER :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_FLT)")
  
  IF(.NOT. ASSOCIATED(VAR%VEC_FLT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (VEC_FLT)")
  END IF

 VALUES =>  VAR%VEC_FLT

END SUBROUTINE NC_POINT_VAR_VEC_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_ARR_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(SPA), POINTER :: values(:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (ARR_FLT)")
  
  IF(.NOT. ASSOCIATED(VAR%ARR_FLT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (ARR_FLT)")
  END IF

  VALUES =>  VAR%ARR_FLT

END SUBROUTINE NC_POINT_VAR_ARR_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_CUB_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(SPA), POINTER :: values(:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (CUB_FLT)")
  
  IF(.NOT. ASSOCIATED(VAR%CUB_FLT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (CUB_FLT)")
  END IF

  VALUES =>  VAR%CUB_FLT

END SUBROUTINE NC_POINT_VAR_CUB_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_FDA_FLT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(SPA), POINTER :: values(:,:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (FDA_FLT)")
  
  IF(.NOT. ASSOCIATED(VAR%FDA_FLT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (FDA_FLT)")
  END IF

  VALUES =>  VAR%FDA_FLT

END SUBROUTINE NC_POINT_VAR_FDA_FLT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_SCL_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(DP), POINTER :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_DBL)")

  IF(.NOT. ASSOCIATED(VAR%SCL_DBL)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (SCL_DBL)")
  END IF

  VALUES => VAR%SCL_DBL

END SUBROUTINE NC_POINT_VAR_SCL_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_VEC_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(DP), POINTER :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_DBL)")

  IF(.NOT. ASSOCIATED(VAR%VEC_DBL)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (VEC_DBL)")
  END IF

  VALUES =>  VAR%VEC_DBL
  
END SUBROUTINE NC_POINT_VAR_VEC_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_ARR_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(DP), POINTER :: values(:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (ARR_DBL)")
  
  IF(.NOT. ASSOCIATED(VAR%ARR_DBL)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (ARR_DBL)")
  END IF

  VALUES =>  VAR%ARR_DBL
  
END SUBROUTINE NC_POINT_VAR_ARR_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_CUB_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(DP), POINTER :: values(:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (CUB_DBL)")
  
    IF(.NOT. ASSOCIATED(VAR%CUB_DBL)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (CUB_DBL)")
  END IF

  VALUES =>  VAR%CUB_DBL
  
END SUBROUTINE NC_POINT_VAR_CUB_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_FDA_DBL(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  REAL(DP), POINTER :: values(:,:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (FDA_DBL)")
  
    IF(.NOT. ASSOCIATED(VAR%FDA_DBL)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (FDA_DBL)")
  END IF

  VALUES =>  VAR%FDA_DBL
  
END SUBROUTINE NC_POINT_VAR_FDA_DBL
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_SCL_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  INTEGER, POINTER :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_INT)")
  
  IF(.NOT. ASSOCIATED(VAR%SCL_INT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (SCL_INT)")
  END IF

  VALUES => VAR%SCL_INT

END SUBROUTINE NC_POINT_VAR_SCL_INT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_VEC_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  INTEGER, POINTER :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_INT)")
 
  IF(.NOT. ASSOCIATED(VAR%VEC_INT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (VEC_INT)")
  END IF
  
  VALUES =>  VAR%VEC_INT
  
END SUBROUTINE NC_POINT_VAR_VEC_INT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_ARR_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  INTEGER, POINTER :: values(:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (ARR_INT)")
  
  IF(.NOT. ASSOCIATED(VAR%ARR_INT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (ARR_INT)")
  END IF
  
  VALUES =>  VAR%ARR_INT
  
END SUBROUTINE NC_POINT_VAR_ARR_INT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_CUB_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  INTEGER, POINTER :: values(:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (CUB_INT)")
  
  IF(.NOT. ASSOCIATED(VAR%CUB_INT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (CUB_INT)")
  END IF

  VALUES =>  VAR%CUB_INT
  
END SUBROUTINE NC_POINT_VAR_CUB_INT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_FDA_INT(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(IN) :: VAR
  INTEGER, POINTER :: values(:,:,:,:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (FDA_INT)")
  
  IF(.NOT. ASSOCIATED(VAR%FDA_INT)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (FDA_INT)")
  END IF

  VALUES =>  VAR%FDA_INT
  
END SUBROUTINE NC_POINT_VAR_FDA_INT
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_SCL_CHR(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  CHARACTER(LEN=80), POINTER :: values
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (SCL_CHR)")

  IF(.NOT. ASSOCIATED(VAR%SCL_CHR)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (SCL_CHR)")
  END IF

   VALUES => VAR%SCL_CHR

 END SUBROUTINE NC_POINT_VAR_SCL_CHR
!====================================================================
!====================================================================
SUBROUTINE NC_POINT_VAR_VEC_CHR(VAR,VALUES)
  IMPLICIT NONE
  TYPE(NCVAR), POINTER, INTENT(INOUT) :: VAR
  CHARACTER(LEN=80), POINTER :: values(:)
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_POINT_VAR: THE VARIABLE ARGUMENT MUST BE ASSOCIATED (VEC_CHR)")

  IF(.NOT. ASSOCIATED(VAR%VEC_CHR)) THEN
     call print_var(var)
     CALL FATAL_ERROR &
          & ("NC_POINT_VAR: THE VARIABLE DATA MUST BE ASSOCIATED (VEC_CHR)")
  END IF

   VALUES => VAR%VEC_CHR

 END SUBROUTINE NC_POINT_VAR_VEC_CHR
!====================================================================
!====================================================================
FUNCTION NC_GET_VAR(NCF,varid) RESULT(VAR)
  implicit none
  TYPE(NCFILE), intent(in) :: NCF
  integer, intent(in) :: varid
  TYPE(NCVAR),pointer :: VAR
  TYPE(NCATT),pointer :: ATT
  TYPE(NCDIM),pointer :: DIM
  integer :: nDims, nAtts, xtype
  integer, allocatable :: dimids(:)
  integer :: status,i
  CHARACTER(LEN=120)     :: errmsg
  LOGICAL FOUND
  CHARACTER(Len=NF90_MAX_NAME+1) :: NAME

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_GET_VAR"

  status=nf90_inquire_variable(NCF%NCID,VARID,& 
       & NAME = NAME,&
       & XTYPE = XTYPE, &
       & NDIMS= nDims, &
       & NATTS= nAtts )
  errmsg="Can not get variable info: "//trim(NCF%FNAME)
  call handle_ncerr(status,errmsg)


  VAR => NEW_VAR()
  var%varname=trim(name)
  VAR%varid=varid
  VAR%NCID => NCF%NCID
  VAR%XTYPE=xtype
  

  if(dbg_set(dbg_io)) write(ipt,*)"====== ADDING VARIABLE ATTRIBUTES: "//trim(VAR%VARNAME)


  allocate(DIMIDS(nDims),stat=status)
  if(status /= 0) CALL FATAL_ERROR("NC_GET_VAR: Can not allocate DIMIDS")

  status=nf90_inquire_variable(NCF%NCID,VAR%VARID, DIMIDS = DIMIDS)
  errmsg="Can not get variable dimids: "//trim(NCF%FNAME)//":"//TRIM(VAR%VARNAME)
  call handle_ncerr(status,errmsg)

  ! ADD THE VARIABLES DIMENSIONS TO ITS LINKED LIST FROM THE FILES DIMENSIONS
  if(dbg_set(dbg_io)) write(ipt,*) "====== ADDING VARIABLE DIMENSIONS:"
  do i=1,nDims
     DIM => FIND_DIM(NCF,DIMIDS(i),FOUND)
     IF(.not. FOUND) THEN
        CALL PRINT_DIM_LIST(NCF)
        CALL FATAL_ERROR("NC_GET_VAR: COULD NOT FIND ONE &
          &OF THE FILE DIMENSION OBJECTS FOR THE VARIABLE: "//TRIM(VAR&
          &%VARNAME), "IN THE FILE: "//trim(NCF%FNAME))
     END IF
 
     if(dbg_set(dbg_io)) write(ipt,*) "      "//trim(DIM%DIMNAME)
     VAR => ADD(VAR,DIM)
  end do
  DEALLOCATE(DIMIDS)

  if(nDims /= count_dim_list(VAR)) then
     if(dbg_set(dbg_log)) call print_dim_list(NCF)
     if(dbg_set(dbg_log)) call print_dim_list(VAR)
     call fatal_error("The number of dimensions in the variable does not m&
          &atch the number loaded in the variable object list.")
  end if


  if(dbg_set(dbg_io)) write(ipt,*) "====== ADDING VARIABLE ATTRIBUTES:"
  do i=1,nAtts
     ATT => NC_GET_ATT(VAR,i)
     if(dbg_set(dbg_io)) write(ipt,*) "      "//trim(ATT%ATTNAME)
     VAR => ADD(VAR,ATT)
  end do

  if(nAtts /= count_att_list(VAR)) then
     if(dbg_set(dbg_log)) call print_att_list(VAR)
     call fatal_error("The number of attributes in the file does not m&
          &atch the number loaded in the file object.")
  end if

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_GET_VAR"
     
END FUNCTION NC_GET_VAR
!====================================================================
!====================================================================
SUBROUTINE NC_DEF_VAR(NCF,varid)
  implicit none
  TYPE(NCFILE), target,intent(in) :: NCF
  integer, intent(in) :: varid
  TYPE(NCVAR),pointer :: VAR
  TYPE(NCATT),pointer :: ATT
  TYPE(NCFILE), pointer :: NCFP

  integer, POINTER :: dimids(:)
  integer :: status,i
  CHARACTER(LEN=120)     :: errmsg
  LOGICAL FOUND

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NC_DEF_VAR"

  VAR => FIND_VAR(NCF,varid,FOUND)
  IF (.NOT. FOUND) THEN
     if (DBG_SET(dbg_log)) call print_var_list(NCF)
     CALL FATAL_ERROR&
          &("NC_DEF_VAR: COULD NOT FIND THE FILE VARIABLE WITH CORRECT VARID W&
          &HILE DEFINING THE VARIABLE IN THE FILE")
  END IF

  ! MAKE AN ARRAY TO HOLD THE DIMIDS  
  DIMIDS => VAR_DIMIDS(VAR)

#if defined DOUBLE_PRECISION && SINGLE_OUTPUT

  IF(VAR%xtype == NF90_DOUBLE) THEN
     
     status = nf90_def_var(NCF%ncid, trim(VAR%varname), NF90_FLOAT,&
          & dimids, i)

  ELSE

     status = nf90_def_var(NCF%ncid, trim(VAR%varname), VAR%xtype,&
          & dimids, i)
  END IF
     
#else
     status = nf90_def_var(NCF%ncid, trim(VAR%varname), VAR%xtype,&
          & dimids, i)
# endif

  errmsg="NF90_DEF_VAR: ERROR"
  CALL HANDLE_NCERR(status,trim(errmsg))

# if defined (USE_NETCDF4) && (USE_COMPRESSION)
  status = nf90_def_var_deflate(NCF%ncid, i,          &
                                shuffle = 1,          &
                                deflate = 1,          &
                                deflate_level = 9  )

  errmsg="NF90_DEF_VAR_DEFLATE: ERROR"
  CALL HANDLE_NCERR(status,trim(errmsg))
# endif

  DEALLOCATE(DIMIDS)
  NULLIFY(DIMIDS)

  IF (VARID .NE. I) THEN
     CALL PRINT_VAR(VAR)
     CALL FATAL_ERROR&
       &("NC_DEF_VAR: THE VARID RETURNED BY NF90_DEF_VAR DOES NOT MATC&
       &H THE VARID FOR THE VARIABLE OBJECT")
  END IF


  Do i = 1,count_att_list(VAR)
     ATT => FIND_ATT(VAR,i,FOUND)
     IF (.NOT. FOUND) THEN
        if (DBG_SET(dbg_log)) call print_att_list(VAR)
        CALL FATAL_ERROR&
             &("NC_DEF_VAR: COULD NOT FIND THE VARIABLE ATTRIBUTE WITH CORRECT ATTID W&
             &HILE PUTTING THE ATTRIBUTE IN THE FILE")
     END IF
     
     CALL WRITE_ATT_TYPE(NCF%NCID,VARID,ATT)
  End Do

  if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NC_DEF_VAR"

END SUBROUTINE NC_DEF_VAR
!====================================================================
!====================================================================
SUBROUTINE NC_WRITE_FILE(NCF,LOCAL_ONLY,STKCNT,STKRNG)
  USE CONTROL
  implicit none
  
  LOGICAL, OPTIONAL :: LOCAL_ONLY
  INTEGER, OPTIONAL :: STKCNT   ! TO SPECIFY THE STACK TO WRITE
  INTEGER, OPTIONAL :: STKRNG(2)! TO SPECIFY THE STACK RANGE TO WRITE

  LOGICAL :: MY_LOCAL
  TYPE(NCFILE), pointer :: NCF

  CHARACTER(LEN=80) :: FNAME,PATH,EXTENSION
  INTEGER :: FCNT,I

  TYPE(NCVAR), pointer :: VAR
  TYPE(NCDIM), pointer :: DIM
  type(NCVARP), POINTER       :: CURRENT
  TYPE(NCFTIME), POINTER :: FTM

  LOGICAL :: COLLECTDATA = .false.
  INTEGER :: COLLECTOR
  INTEGER :: NEXT_STKCNT
  LOGICAL :: DUMP = .false.
  LOGICAL :: FOUND = .false.

 if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "STARTING NC_WRITE_FILE"

 if(DBG_SET(dbg_sbrio)) CALL PRINT_FILE(NCF)


 MY_LOCAL = .false.
 IF(PRESENT(LOCAL_ONLY)) MY_LOCAL = LOCAL_ONLY

  IF(SERIAL .or. MY_LOCAL) THEN
     COLLECTDATA = .false.
     DUMP    = .true.
     COLLECTOR = MSRID
     
  ELSE ! (IF PARALLEL)
     
     COLLECTDATA = .true.
     
     IF (USE_MPI_IO_MODE) then
        COLLECTOR = IOPROCID
        DUMP = .false.
     ELSE
        COLLECTOR = MSRID ! MASTER
        IF(MSR)THEN
           DUMP = .true.
        ELSE
           DUMP = .false.
        END IF
     END IF
     
  END IF

  IF(PRESENT(STKCNT) .AND. PRESENT(STKRNG)) CALL FATAL_ERROR&
       &("NC_WRITE_FILE: CAN NOT CALL WITH BOTH IDX AND RNG!")


  if(DBG_SET(dbg_sbrio)) then
     write(IPT,*)"==========================================="
     write(IPT,*)"= NC_WRITE_FILE:         CONTROL STATE"
     write(IPT,*)"= NC_WRITE_FILE:   COLLECT DATA:",COLLECTDATA
     write(IPT,*)"= NC_WRITE_FILE:           DUMP:",DUMP
     write(IPT,*)"= NC_WRITE_FILE:      COLLECTOR:",COLLECTOR
     write(IPT,*)"= NC_WRITE_FILE:           MYID:",MYID
     write(IPT,*)"= NC_WRITE_FILE:         IOPROC:",IOPROC
     write(IPT,*)"= NC_WRITE_FILE:USE_MPI_IO_MODE:",USE_MPI_IO_MODE
     IF(PRESENT(LOCAL_ONLY))THEN
        write(IPT,*)"= NC_WRITE_FILE:     LOCAL_ONLY:",LOCAL_ONLY
     END IF
     IF(PRESENT(STKCNT))THEN
        write(IPT,*)"= NC_WRITE_FILE:         STKCNT:",STKCNT
     END IF
     IF(PRESENT(STKRNG))THEN
        write(IPT,*)"= NC_WRITE_FILE:         STKRNG:",STKRNG
     END IF

     write(IPT,*)"==========================================="
  end if

  IF(LEN_TRIM(NCF%FNAME) == 0) THEN
     CALL PRINT_FILE(NCF)
     CALL FATAL_ERROR("NC_WRITE_FILE: CALLED WRITE WITH NO FILE NAME!")
  END IF
  
  CALL PATH_SPLIT(NCF%FNAME,PATH,FNAME,EXTENSION)

!!$  IF (COLLECTOR .EQ. MYID) THEN
  IF (COLLECTOR .EQ. MYID_iogroup) THEN
     
     IF(NCF%CONNECTED) THEN
        
        CALL NC_OPEN(NCF)

        ! ADD CHECK FOR STK_LEN VS FILE UNLIMDIM?
     ELSE        
       
        IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! CREATING FILE: "//trim(NCF%FNAME)

        CALL NC_CREATE(NCF)
        CALL NC_SAVE(NCF)
     END IF
  ELSE

     NCF%CONNECTED = .TRUE.

  END IF

  IF(ASSOCIATED(NCF%FTIME)) THEN
     FTM => NCF%FTIME
     
     IF(PRESENT(STKCNT) )THEN

        NEXT_STKCNT = STKCNT
        FTM%NEXT_STKCNT = STKCNT+1
        FTM%PREV_STKCNT =STKCNT

     ELSEIF(PRESENT(STKRNG))THEN
        NEXT_STKCNT = STKRNG(2)
        FTM%NEXT_STKCNT = STKRNG(2)
        FTM%PREV_STKCNT = STKRNG(1)

     ELSE ! DEFAULT CASE - NO SPECIFED FILE STACK, USE FTM%NEXT_STKCNT
        
        NEXT_STKCNT = FTM%NEXT_STKCNT

     END IF

     IF(FTM%NEXT_STKCNT .LT. 0) CALL FATAL_ERROR &
          &("NC_WRITE_FILE: FILE OBJECT STKCNT LESS THAN ZERO",&
          &"FILE NAME: "//trim(NCF%FNAME))

        
     ! IF WE ARE MAKING THE FILE LONGER INCRIMENT THE STK_LEN
     FTM%STK_LEN = MAX(FTM%NEXT_STKCNT,FTM%STK_LEN)
  ELSE

     IF(PRESENT(STKCNT) .OR. PRESENT(STKRNG)) CALL FATAL_ERROR&
          &("IT IS NONSENSE TO PASS A STK OR RNG TO NC_WRITE_FILE",&
          & "WITH A FILE THAT DOES NOT HAVE AN ASSOCIATED FILETIME!")

     ! NOT A TIME VARRYING FILE... WRITE IT ANYWAY!
     NEXT_STKCNT = 0


!     CALL PRINT_FILE(NCF)
!     CALL FATAL_ERROR ("NC_WRITE_FILE: FILE OBJECT FTIME IS NOT ALLOCATED ")

  END IF


  IF (NEXT_STKCNT == 0) THEN
     IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! DUMPING STATIC DATA TO FILE: "&
          &//trim(FNAME)
  ELSE
     
     IF(PRESENT(STKRNG)) THEN
        
        IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! DUMPING DATA TO FILE: "&
             &//trim(FNAME)//"; Records#",STKRNG
     ELSE
        IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! DUMPING DATA TO FILE: "&
             &//trim(FNAME)//"; Record#",NEXT_STKCNT
     END IF
  END IF
  

  CURRENT => NCF%VARS%NEXT

  IF(.NOT. ASSOCIATED(CURRENT)) &
       & CALL FATAL_ERROR("NC_WRITE_FILE: FILE OBJECT HAS NO VARIABLES",&
       &"FILE NAME: "//trim(NCF%FNAME))


! LOOP THROUGH VARIABLES AND WRITE THE DATA
! ===================================================================== 
  DO
     IF(.NOT. ASSOCIATED(CURRENT)) THEN

        IF (IOPROC .AND. COLLECTDATA) THEN
           ! JUST FINISHED COLLECTING - NOW DUMP
           COLLECTDATA = .false.
           DUMP    = .true.
           !START LOOP OVER AND DUMP THE DATA
           CURRENT => NCF%VARS%NEXT        
        ELSE
!            YOUR DONE!
!                FVCOM GROUP ALWAYS GOES THROUGH THE LIST ONCE
!                IOPROC GOES THROUGH TWICE - COLLECT ONCE
!                                          - DUMP ONCE
!!$           IF (COLLECTOR .EQ. MYID) CALL NC_CLOSE(NCF)
           IF (COLLECTOR .EQ. MYID_iogroup) CALL NC_CLOSE(NCF)

           if(DBG_SET(dbg_sbr)) &
                & write(IPT,*) "END NC_WRITE_FILE"

           RETURN
        END IF
        
     END IF

     IF(.NOT. ASSOCIATED(CURRENT%VAR)) CALL FATAL_ERROR &
          &("NC_WRITE_FILE:",&
          & "FILE OBJECT HAS UNASSOCIATED VARIBLE OBJECT IN LINK LIST",&
          &"FILE NAME: "//trim(NCF%FNAME))

     VAR => CURRENT%VAR
        

     ! DECIDE WHETHER THIS IS A DYNAMIC OR STATIC VARIABLE     
     FOUND = HAS_UNLIMITED(VAR)

     IF (FOUND .AND. (NEXT_STKCNT .GT. 0) ) THEN 
        ! THIS IS A TIME STEP OUT PUT - WRITE THE VARIABLE
        IF(PRESENT(STKRNG)) THEN
           CALL NC_WRITE_VAR(VAR,DUMP,COLLECTDATA,COLLECTOR,STKRNG=STKRNG)
        ELSE
           CALL NC_WRITE_VAR(VAR,DUMP,COLLECTDATA,COLLECTOR,STKCNT=NEXT_STKCNT)
        END IF
     ELSEIF (.NOT. FOUND .AND. (NEXT_STKCNT .EQ. 0) ) THEN 
        ! STKCNT == 0 - JUST WROTE NEW FILE ADD STATIC VARIABLES
        CALL NC_WRITE_VAR(VAR,DUMP,COLLECTDATA,COLLECTOR)
     END IF

     CURRENT => CURRENT%NEXT
  END DO

  CALL FATAL_ERROR("NC_WRITE_FILE REACHED AN IMPOSSIBLE STATE",&
       &"PLEASE SET YOUR COMPUTER ON FIRE AND EXIT THE BUILDING QUICKLY")

END SUBROUTINE NC_WRITE_FILE
!====================================================================
!====================================================================
! Watch out - pointing to part of an array is tricky. The index may be reset
SUBROUTINE NC_WRITE_VAR(VAR,DUMP,COLLECTDATA,COLLECTOR,STKCNT,STKRNG,IOSTART,IOCOUNT,IOSTRIDE)
! DUMP - If T, write the data to the disk
! COLLECTDATA - If T, Collect data to the processor specified by COLLECTOR
! COLLECTOR - The Processor which does write/collect

  USE CONTROL, only:msr,ioproc,use_mpi_io_mode,SERIAL
  implicit none
  TYPE(NCVAR), POINTER :: VAR
  INTEGER, OPTIONAL :: STKCNT
  INTEGER, OPTIONAL :: STKRNG(2)
  INTEGER, ALLOCATABLE,TARGET, OPTIONAL :: IOSTART(:), IOCOUNT(:), IOSTRIDE(:)
  LOGICAL, INTENT(IN) ::DUMP
  LOGICAL, INTENT(IN) :: COLLECTDATA
  INTEGER, INTENT(IN) :: COLLECTOR

  INTEGER :: CODE
  INTEGER :: XTYPE
  INTEGER :: CNT,DIMCNT
  INTEGER :: NSIZE
!  integer :: dim1,dim2,dim3,dim4

 ! NF90 SUBSET VARIABLES
  INTEGER, POINTER :: NSTART(:), NCOUNT(:), NSTRIDE(:)
#if defined (PIO)
 ! NF90 SUBSET VARIABLES FOR PARALLEL IO  
  INTEGER, ALLOCATABLE :: NSTART_PIO(:), NCOUNT_PIO(:), NSTRIDE_PIO(:)
  INTEGER PARTITION_PIO,rest_PIO
#endif  
  ! The size of the data returned by nf90_get/put_var
  INTEGER, POINTER :: RDIMS(:)
  ! TEMPORARY VARIABLES ONLY USED FOR VEC_CHR
  INTEGER, POINTER :: NSTRT(:), NCNT(:), NSTRD(:)

  ! POINTERS AND SUCH
  TYPE(NCDIM), POINTER :: DIM
  TYPE(NCDIMP), POINTER :: DIMLINK
  LOGICAL :: FOUND

# if defined(MULTIPROCESSOR)
  !MPI COMM STUFF
  INTEGER, PARAMETER :: WVD_TAG = 40001
  INTEGER :: DEST, SOURCE, IERR,TMPID
  INTEGER :: STAT(MPI_STATUS_SIZE)
  TYPE(MAP), POINTER :: GMAP(:)
  INTEGER, POINTER :: LSizes(:), NPsize(:) 
# endif

  INTEGER, PARAMETER :: case_scl_int = 1
  INTEGER, PARAMETER :: case_vec_int = 2
  INTEGER, PARAMETER :: case_arr_int = 3
  INTEGER, PARAMETER :: case_cub_int = 4
  INTEGER, PARAMETER :: case_fda_int = 5

  INTEGER, PARAMETER :: case_scl_flt = 6
  INTEGER, PARAMETER :: case_vec_flt = 7
  INTEGER, PARAMETER :: case_arr_flt = 8
  INTEGER, PARAMETER :: case_cub_flt = 9
  INTEGER, PARAMETER :: case_fda_flt = 10

  INTEGER, PARAMETER :: case_scl_dbl = 11
  INTEGER, PARAMETER :: case_vec_dbl = 12
  INTEGER, PARAMETER :: case_arr_dbl = 13
  INTEGER, PARAMETER :: case_cub_dbl = 14
  INTEGER, PARAMETER :: case_fda_dbl = 15

  INTEGER, PARAMETER :: case_scl_chr = 16
  INTEGER, PARAMETER :: case_vec_chr = 17

  ! TEMPORARY STORAGE FOR DATA IF COLLECTED TO MASTER PROC
  INTEGER, POINTER                :: SCL_INT
  INTEGER, POINTER,DIMENSION(:)   :: GVEC_INT
  INTEGER, POINTER,DIMENSION(:,:) :: GARR_INT
  INTEGER, POINTER,DIMENSION(:,:,:) :: GCUB_INT
  INTEGER, POINTER,DIMENSION(:,:,:,:) :: GFDA_INT
  INTEGER, ALLOCATABLE :: GVEC_INT_temp(:),GARR_INT_temp(:,:),GCUB_INT_temp(:,:,:),GFDA_INT_temp(:,:,:,:)

  INTEGER, POINTER,DIMENSION(:)   :: LVEC_INT
  INTEGER, POINTER,DIMENSION(:,:) :: LARR_INT
  INTEGER, POINTER,DIMENSION(:,:,:) :: LCUB_INT
  INTEGER, POINTER,DIMENSION(:,:,:,:) :: LFDA_INT

  
  REAL(SPA), POINTER                :: SCL_FLT
  REAL(SPA), POINTER,DIMENSION(:)   :: LVEC_FLT
  REAL(SPA), POINTER,DIMENSION(:,:) :: LARR_FLT
  REAL(SPA), POINTER,DIMENSION(:,:,:) :: LCUB_FLT
  REAL(SPA), POINTER,DIMENSION(:,:,:,:) :: LFDA_FLT

  REAL(SPA), POINTER,DIMENSION(:)   :: GVEC_FLT
  REAL(SPA), POINTER,DIMENSION(:,:) :: GARR_FLT
  REAL(SPA), POINTER,DIMENSION(:,:,:) :: GCUB_FLT
  REAL(SPA), POINTER,DIMENSION(:,:,:,:) :: GFDA_FLT
  REAL(SPA), ALLOCATABLE :: GVEC_FLT_temp(:),GARR_FLT_temp(:,:),GCUB_FLT_temp(:,:,:),GFDA_FLT_temp(:,:,:,:)
  
  REAL(DP), POINTER                :: SCL_DBL
  REAL(DP), POINTER,DIMENSION(:)   :: GVEC_DBL
  REAL(DP), POINTER,DIMENSION(:,:) :: GARR_DBL
  REAL(DP), POINTER,DIMENSION(:,:,:) :: GCUB_DBL
  REAL(DP), POINTER,DIMENSION(:,:,:,:) :: GFDA_DBL

  REAL(DP), POINTER,DIMENSION(:)   :: LVEC_DBL
  REAL(DP), POINTER,DIMENSION(:,:) :: LARR_DBL
  REAL(DP), POINTER,DIMENSION(:,:,:) :: LCUB_DBL
  REAL(DP), POINTER,DIMENSION(:,:,:,:) :: LFDA_DBL
  REAL(DP), ALLOCATABLE :: GVEC_DBL_temp(:),GARR_DBL_temp(:,:),GCUB_DBL_temp(:,:,:),GFDA_DBL_temp(:,:,:,:)

  CHARACTER(LEN=80), POINTER              :: SCL_CHR
  CHARACTER(LEN=80), POINTER,DIMENSION(:) :: VEC_CHR

  CHARACTER(len=3) :: char1,char2,char3
  ! DATA FOR PUT VAR COMMANDS:
  INTEGER :: STATUS, I
  CHARACTER(LEN=120)          :: errmsg

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START NC_WRITE_VAR:"

  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_WRITE_VAR: Variable object argument is not assocaited!")

  ! INITIALIZE SOME MEMORY
  STATUS=0
  I=0
  
  CODE=0
  XTYPE=0
  CNT=0
  DIMCNT=0
  NSIZE=0
  FOUND=.FALSE.
  
# if defined(MULTIPROCESSOR)
  NULLIFY(GMAP)
  IERR=0
  SOURCE=0
  DEST=0
# endif
  !NULLIFY POINTERS
  NULLIFY(NSTART,NCOUNT,NSTRIDE)
  NULLIFY(RDIMS,NSTRT,NCNT,NSTRD)
  NULLIFY(DIM,DIMLINK)
  NULLIFY(SCL_INT,GVEC_INT,GARR_INT,GCUB_INT,GFDA_INT,LVEC_INT,LARR_INT,LCUB_INT,LFDA_INT)
  
  NULLIFY(SCL_FLT,GVEC_FLT,GARR_FLT,GCUB_FLT,GFDA_FLT,LVEC_FLT,LARR_FLT,LCUB_FLT,LFDA_FLT)
  NULLIFY(SCL_DBL,GVEC_DBL,GARR_DBL,GCUB_DBL,GFDA_DBL,LVEC_DBL,LARR_DBL,LCUB_DBL,LFDA_DBL)
  NULLIFY(SCL_CHR,VEC_CHR)
  

  ! COUNT THE NUMBER OF DIMENSIONS IN THE VARIABLE
  DIMCNT = count_dim_list(VAR)


  IF (SERIAL .and. (COLLECTOR .NE. MYID .OR. COLLECTDATA)) THEN
     CALL PRINT_VAR(VAR)
     CALL FATAL_ERROR("NC_WRITE_VAR: SERIAL JOB CALLED A PARALLEL WRITE?")
  END IF

  IF (.NOT. DUMP .and. .NOT. COLLECTDATA) THEN
     CALL PRINT_VAR(VAR)
     CALL FATAL_ERROR("NC_WRITE_VAR: CALLED WITH BAD ARGUMENTS;",&
          & "DUMP or COLLECTDATA or both must be true?")
  END IF

  IF(DBG_SET(DBG_SBRIO)) THEN

        write(char2,'(I3.3)')collector
        write(char3,'(I3.3)')myid
        WRITE(IPT,*)"NC_WRITE_VAR Arguments:"
        call print_var(var)
        WRITE(IPT,*)"DUMP=",DUMP,"; COLLECTDATA=",COLLECTDATA,"; COLLECTOR="//char2//"; MYID="//char3
        IF(PRESENT(STKCNT)) THEN
           WRITE(IPT,*) "STKCNT=",STKCNT
        ELSE
           WRITE(IPT,*) "STKCNT= NONE"
        END IF

        IF(PRESENT(STKRNG)) THEN
           WRITE(IPT,*) "STKRNG=",STKRNG
        ELSE
           WRITE(IPT,*) "STKRNG= NONE"
        END IF

        IF(PRESENT(IOSTART)) THEN
           WRITE(IPT,*) "IOSTART=",IOSTART
        ELSE
           WRITE(IPT,*) "IOSTART= NONE"
        END IF

        IF(PRESENT(IOCOUNT)) THEN
           WRITE(IPT,*) "IOCOUNT=",IOCOUNT
        ELSE
           WRITE(IPT,*) "IOCOUNT= NONE"
        END IF

        IF(PRESENT(IOSTRIDE)) THEN
           WRITE(IPT,*) "IOSTRIDE=",IOSTRIDE
        ELSE
           WRITE(IPT,*) "IOSTRIDE= NONE"
        END IF
        
  END IF


  IF(VAR%NCID == -1 .and. DUMP) THEN
     CALL PRINT_VAR(VAR)
     CALL FATAL_ERROR("NC_WRITE_VAR: CAN NOT WRITE TO FILE, IT IS NOT OPEN!")
  END IF


  IF ( PRESENT(STKCNT) ) THEN

     IF ( PRESENT(STKRNG) .or. PRESENT(IOSTART) .or. PRESENT(IOCOUNT) .or. PRESENT(IOSTRIDE))THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_WRITE_VAR: You can not pass both STKCNT and STKRNG/START/COUNT/STRIDE !",&
             &"Set STKCNT to write a time slice filling all other dimensions. OR",&
             &"Set IOSTART/IOCOUNT/(IOSTRIDE) to read a specific range.")
     END IF

     DIM => FIND_UNLIMITED(VAR,FOUND)
     IF(.NOT.FOUND) THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR&
             &("NC_WRITE_VAR: CALLED WITH STKCNT ARGUMENT BUT VARIABLE IS NOT UNLIMITED?")
     END IF
     NULLIFY(DIM)
     
     ALLOCATE(NSTART(DIMCNT),NCOUNT(DIMCNT),NSTRIDE(DIMCNT))
     
     
     ! SET THE VARIABLES CURRENT STACK COUNT
     VAR%CURR_STKCNT = STKCNT
     
     ! SET THE NF90_PUT_VAR DIMENSIONS
     NSTART=1 ! START AT ONE, EXCEPT FOR THE TIME VARIABLE
     NSTART(DIMCNT) = STKCNT

     
     DIMLINK => VAR%DIMS%NEXT
     DO I = 1,DIMCNT ! GET THE OUTPUT VARIABLE DIMENSIONS
        DIM => DIMLINK%DIM  
        NCOUNT(I)= DIM%DIM
        DIMLINK => DIMLINK%NEXT
     END DO
     NCOUNT(DIMCNT)=1 ! SET THE TIME OUTPUT DIMENSION TO 1
     
     NSTRIDE=1 ! ALWAYS USE STRIDE 1 FOR STKCNT INPUT

     ! FOR TIME DEPENDANT DATA THE RANK OF THE ALLOCATED MEMORY IS
     ! ONE LESS THAN THE RANK OF THE FILE's VARIALBE!
     DIMCNT = DIMCNT -1

     ! THE DIMENSIONS WRITTEN WILL BE THE VALUES OF NCOUNT, NOT
     ! INCLUDING TIME
     IF(DIMCNT > 0)THEN
        ALLOCATE(RDIMS(DIMCNT))
        RDIMS(1:DIMCNT)=NCOUNT(1:DIMCNT)
     ELSE IF (DIMCNT == 0) THEN
         ALLOCATE(RDIMS(1))
        RDIMS(1)=NCOUNT(1)
     ELSE
        nullify(RDIMS)
     END IF

     ELSEIF ( PRESENT(STKRNG) ) THEN

     IF ( PRESENT(STKCNT) .or. PRESENT(IOSTART) .or. PRESENT(IOCOUNT) .or. PRESENT(IOSTRIDE))THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_WRITE_VAR: You can not pass both STKRNG and STKCNT/START/COUNT/STRIDE !",&
             &"Set STKRNG to write a time range filling all other dimensions. OR",&
             &"Set IOSTART/IOCOUNT/(IOSTRIDE) to read a specific range.")
     END IF

     DIM => FIND_UNLIMITED(VAR,FOUND)
     IF(.NOT.FOUND) THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR&
             &("NC_WRITE_VAR: CALLED WITH STKRNG ARGUMENT BUT VARIABLE IS NOT UNLIMITED?")
     END IF
     NULLIFY(DIM)
     
     ALLOCATE(NSTART(DIMCNT),NCOUNT(DIMCNT),NSTRIDE(DIMCNT))
     
     
     ! SET THE VARIABLES CURRENT STACK COUNT
     VAR%CURR_STKCNT = -1
     
     ! SET THE NF90_PUT_VAR DIMENSIONS
     NSTART=1 ! START AT ONE, EXCEPT FOR THE TIME VARIABLE
     NSTART(DIMCNT) = STKRNG(1)

     
     DIMLINK => VAR%DIMS%NEXT
     DO I = 1,DIMCNT ! GET THE OUTPUT VARIABLE DIMENSIONS
        DIM => DIMLINK%DIM  
        NCOUNT(I)= DIM%DIM
        DIMLINK => DIMLINK%NEXT
     END DO
     NCOUNT(DIMCNT)=STKRNG(2) - STKRNG(1)+1 ! SET THE TIME OUTPUT
     ! DIMENSION TO NUMBER OF STATES IN THE RANGE
     
     NSTRIDE=1 ! ALWAYS USE STRIDE 1 FOR STKCNT INPUT

     ! THE DIMENSIONS WRITTEN WILL BE THE VALUES OF NCOUNT, NOT
     ! INCLUDING TIME
     RDIMS=>NCOUNT

  ELSE IF( PRESENT(IOSTART) .and. PRESENT(IOCOUNT)) THEN
     
     NSTART=>IOSTART
     NCOUNT=>IOCOUNT
     
     IF(.not. PRESENT(IOSTRIDE)) THEN
        ALLOCATE(NSTRIDE(DIMCNT))
        NSTRIDE=1
     ELSE
        NSTRIDE=>IOSTRIDE
     END IF

     IF(DIMCNT /= size(NSTART) .or. &
          & DIMCNT /= size(NCOUNT) .or. &
          & DIMCNT /= size(NSTRIDE) ) THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR&
        & ("NC_WRITE_VAR: Variable's # of file dimensions does not matach size(NSTART/NCOUNT/NSTRIDE) arugments?")
     END IF

     ! SET THE VARIABLES CURRENT STACK COUNT: not defined for this
     ! kind of read/write
     VAR%CURR_STKCNT = -1

     ! ONLY COUNT THE NONE SINGLETON DIMENSIONS OF A VARIABLE. 
     CNT = 0
     DO I = 1,DIMCNT
        IF(NCOUNT(I)>1) CNT=CNT+1
     END DO     
     ! NOW RECORD THE DIMENSIONS OF THE DATA THAT WILL BE READ INTO MEMORY
     IF (CNT > 0) THEN
        
        ALLOCATE(RDIMS(CNT))
        CNT = 0
        DO I = 1,DIMCNT
           IF(NCOUNT(I)>1) THEN
              CNT=CNT+1
              RDIMS(CNT)=NCOUNT(I)
           END IF
        END DO

     ELSE
        ALLOCATE(RDIMS(1)) 
        RDIMS(1) = 1
     END IF

     
     ! NOW SET THE DIMENSION OF THE DATA VARIABLE IN MEMORY     
     DIMCNT=CNT

  ELSE IF( .not. (PRESENT(IOSTART) .or. PRESENT(IOCOUNT) .or.&
       & PRESENT(STKCNT) .or. PRESENT(STKRNG) .or. PRESENT(IOSTRIDE))) THEN
     
     ALLOCATE(NSTART(DIMCNT),NCOUNT(DIMCNT),NSTRIDE(DIMCNT))
     
     ! SET THE VARIABLES CURRENT STACK COUNT
     VAR%CURR_STKCNT = 0
     
     ! SET THE NF90_PUT_VAR DIMENSIONS
     NSTART=1 ! START AT ONE, EXCEPT FOR THE TIME VARIABLE
     
     DIMLINK => VAR%DIMS%NEXT
     DO I = 1,DIMCNT ! GET THE OUTPUT VARIABLE DIMENSIONS
        DIM => DIMLINK%DIM  
        NCOUNT(I)= DIM%DIM
        DIMLINK => DIMLINK%NEXT
     END DO
     
     NSTRIDE=1 ! ALWAYS USE STRIDE 1 IF NO ARGUMENTS ARE PASSED

     ! ONLY COUNT THE NONE SINGLETON DIMENSIONS OF A VARIABLE. 
     CNT = 0
     DO I = 1,DIMCNT
        IF(NCOUNT(I)>1) CNT=CNT+1
     END DO     
     ! NOW RECORD THE DIMENSIONS OF THE DATA THAT WILL BE READ INTO MEMORY
     IF (CNT > 0) THEN
        
        ALLOCATE(RDIMS(CNT))
        CNT = 0
        DO I = 1,DIMCNT
           IF(NCOUNT(I)>1) THEN
              CNT=CNT+1
              RDIMS(CNT)=NCOUNT(I)
           END IF
        END DO

     ELSE
        ALLOCATE(RDIMS(1)) 
        RDIMS(1) = 1
     END IF

     
     ! NOW SET THE DIMENSION OF THE DATA VARIABLE IN MEMORY     
     DIMCNT=CNT

  ELSE
     
     IF(DBG_SET(DBG_LOG)) THEN
        write(ipt,*) "# IOSTART  ::",PRESENT(IOSTART)
        write(ipt,*) "# IOCOUNT  ::",PRESENT(IOCOUNT)
        write(ipt,*) "# IOSTRIDE ::",PRESENT(IOSTRIDE)
        write(ipt,*) "# STKCNT   ::",PRESENT(STKCNT)
        write(ipt,*) "# STKRNG   ::",PRESENT(STKRNG)
     END IF
     
     CALL FATAL_ERROR("NC_WRITE_VAR: YOU SPECIFIED AN ILLEGAL COMBINATION OF AGUMENTS?",&
          & "Valid choices are STKCNT or STKRNG or NSTART,NCOUNT,(NSTRIDE) or none")
  END IF

  IF(DBG_SET(DBG_SBRIO)) THEN
     write(IPT,*) "FILE DIMENSION COUNT IS::", count_dim_list(VAR)
     write(IPT,*) "MEMORY DIMENSION COUNT IS::",DIMCNT
  END IF

  !DETERMIN WHICH CASE WE ARE WRITING DATA FOR
  code = -1

  select case(VAR%XTYPE)
  case(NF90_BYTE)
     call Fatal_error("NC_WRITE_VAR: NOT SET UP TO WRITE BYTE DATA")
  case(NF90_SHORT)
     call Fatal_error("NC_WRITE_VAR: NOT SET UP TO WRITE SHORT DATA")
  case(NF90_INT)
     if (DIMCNT == 0) CODE = case_scl_int
     if (DIMCNT == 1) CODE = case_vec_int
     if (DIMCNT == 2) CODE = case_arr_int
     if (DIMCNT == 3) CODE = case_cub_int
     if (DIMCNT == 4) CODE = case_fda_int
  case(NF90_FLOAT)
     if (DIMCNT == 0) CODE = case_scl_flt
     if (DIMCNT == 1) CODE = case_vec_flt
     if (DIMCNT == 2) CODE = case_arr_flt
     if (DIMCNT == 3) CODE = case_cub_flt
     if (DIMCNT == 4) CODE = case_fda_flt
  case(NF90_DOUBLE)
     if (DIMCNT == 0) CODE = case_scl_dbl
     if (DIMCNT == 1) CODE = case_vec_dbl
     if (DIMCNT == 2) CODE = case_arr_dbl
     if (DIMCNT == 3) CODE = case_cub_dbl
     if (DIMCNT == 4) CODE = case_fda_dbl
  case(NF90_CHAR)
     
     IF(NCOUNT(1) == 1) THEN
        WRITE(IPT,*) "SINGLETON CHARACTER DATA!"
        IF(.not. ASSOCIATED(RDIMS,NCOUNT)) THEN
           DEALLOCATE(RDIMS)
        ELSE
           NULLIFY(RDIMS)
        END IF

        DIMCNT = DIMCNT+1

        ALLOCATE(RDIMS(DIMCNT))
        CNT = 1
        RDIMS(1) = NCOUNT(1)
        DO I = 2,size(ncount)
           IF(NCOUNT(I)>1) THEN
              CNT=CNT+1
              RDIMS(CNT)=NCOUNT(I)
           END IF
        END DO
     END IF

     if (DIMCNT == 1) CODE = case_scl_chr
     if (DIMCNT == 2) CODE = case_vec_chr
     ! First dim is length of string
     ! Second dim is time
  case default
     call Fatal_error("NC_WRITE_VAR: Unkown data type?")
  end select

  ! BASED ON CODE WRITE THE DATA
  errmsg="NC_WRITE_VAR: VARIABLE: "//VAR%varname//"; Can not be writen by nf90_put_var!"

  SELECT CASE(CODE)
!********************************************************************
! =====   SCALAR INTEGER DATA
!********************************************************************
  CASE(case_scl_int)

     IF(.NOT. ASSOCIATED(VAR%SCL_INT))THEN 

        IF(ASSOCIATED(VAR%VEC_INT))THEN
           IF(size(VAR%VEC_INT)==1) VAR%SCL_INT=>VAR%VEC_INT(1)
        ELSE IF(ASSOCIATED(VAR%ARR_INT))THEN
           IF(size(VAR%ARR_INT)==1) VAR%SCL_INT=>VAR%ARR_INT(1,1)
        ELSE IF(ASSOCIATED(VAR%CUB_INT))THEN
           IF(size(VAR%CUB_INT)==1) VAR%SCL_INT=>VAR%CUB_INT(1,1,1)
        ELSE IF(ASSOCIATED(VAR%FDA_INT))THEN
           IF(size(VAR%FDA_INT)==1) VAR%SCL_INT=>VAR%FDA_INT(1,1,1,1)
        ELSE        

           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects SCL_INT data is NOT assocaited!")
        END IF
     END IF

     SCL_INT => VAR%SCL_INT

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA .AND. USE_MPI_IO_MODE) THEN
# if defined(MULTIPROCESSOR)
        DEST = IOPROCID - 1
!!$        SOURCE = MSRID - 1
        SOURCE = 0
        NSIZE = 1
!!$        IF (MSR) CALL MPI_SEND&
!!$             & (scl_int,NSIZE,MPI_INTEGER,DEST,WVD_TAG,MPI_FVCOM_GROUP,IERR)
!!$        IF (IOPROC) CALL MPI_RECV &
!!$             & (scl_int,NSIZE,MPI_INTEGER,SOURCE,WVD_TAG,MPI_FVCOM_GROUP,STAT,IERR)
        IF (MYID_iogroup == 1) CALL MPI_SEND&
             & (scl_int,NSIZE,MPI_INTEGER,DEST,WVD_TAG,MPI_IO_GROUP,IERR)
        IF (IOPROC) CALL MPI_RECV &
             & (scl_int,NSIZE,MPI_INTEGER,SOURCE,WVD_TAG,MPI_IO_GROUP,STAT,IERR)
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

     IF (DUMP) THEN

        IF (SIZE(NSTART) .GT. 0) THEN
            
           if (product(nCOUNT) .NE. 1) CALL FATAL_ERROR&
                & ("NC_WRITE_VAR: NCOUNT dimension invalid while reading scl_int?")

           ! ARGUMENT TO NF90_PUT_VAR MUST BE A VECTOR
           allocate(GVEC_INT(1)); GVEC_INT(1) = SCL_INT
           IF(MYID_fgroup == 1) THEN
             STATUS = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_INT,NSTART,NCOUNT,NSTRIDE)
             CALL HANDLE_NCERR(status,trim(errmsg))
           ENDIF
           
           deallocate(GVEC_INT)
        ELSE
           IF(MYID_fgroup == 1) THEN
             STATUS = NF90_PUT_VAR(VAR%NCID,VAR%VARID,SCL_INT)
             CALL HANDLE_NCERR(status,trim(errmsg))
           ENDIF
        END IF

     END IF
     
     NULLIFY(SCL_INT)

!********************************************************************
! =====   VECTOR INTEGER DATA
!********************************************************************
  CASE(case_vec_int)

    IF(.NOT. ASSOCIATED(VAR%VEC_INT))THEN 

        IF(ASSOCIATED(VAR%ARR_INT))THEN
           IF(size(VAR%ARR_INT,1)==1) VAR%VEC_INT=>VAR%ARR_INT(1,1:)
           IF(size(VAR%ARR_INT,2)==1) VAR%VEC_INT=>VAR%ARR_INT(1:,1)
        ELSE IF(ASSOCIATED(VAR%CUB_INT))THEN
           IF(size(VAR%CUB_INT,1)==1) THEN
              IF(size(VAR%CUB_INT,2)==1) VAR%VEC_INT=>VAR%CUB_INT(1,1,1:)
              IF(size(VAR%CUB_INT,3)==1) VAR%VEC_INT=>VAR%CUB_INT(1,1:,1)
           END IF
           IF(size(VAR%CUB_INT,1)==2) THEN
              IF(size(VAR%CUB_INT,3)==1) VAR%VEC_INT=>VAR%CUB_INT(1:,1,1)
           END IF           
!        ELSE IF(ASSOCIATED(VAR%FDA_INT))THEN
!           IF(size(VAR%FDA_INT,1)==1) THEN
!              IF(size(VAR%FDA_INT,2)==1) VAR%VEC_INT=>VAR%FDA_INT(1,1,1,1:)
!              IF(size(VAR%FDA_INT,3)==1) VAR%VEC_INT=>VAR%FDA_INT(1,1,1:,1)
!              IF(size(VAR%FDA_INT,4)==1) VAR%VEC_INT=>VAR%FDA_INT(1,1:,1,1)
!           END IF
!           IF(size(VAR%FDA_INT,1)==2) THEN
!              IF(size(VAR%FDA_INT,4)==1) VAR%VEC_INT=>VAR%FDA_INT(1:,1,1,1)
!           END IF           
        ELSE        

           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects VEC_FLT data is NOT assocaited!")
        END IF
     END IF

     nsize=ubound(VAR%VEC_INT,1)

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)

        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_IO_GROUP
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,ierr)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE

        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)

        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           ALLOCATE(GVEC_INT(RDIMS(1)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR: Allocate VEC_INT failed!")                   
 
 
           LVEC_INT => VAR%VEC_INT(1:nsize)
           
        ELSE

           IF (IOPROC) THEN
              IF (UBOUND(VAR%VEC_INT,1) .LT. RDIMS(1)) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              GVEC_INT => VAR%VEC_INT(1:RDIMS(1))
	      GVEC_INT = -HUGE(GVEC_INT) 
           ELSE

              LVEC_INT => VAR%VEC_INT
           END IF

        END IF

!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LVEC_INT,GVEC_INT)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LVEC_INT,GVEC_INT)
        
 
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN

        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
!!$        If(.NOT. COLLECTDATA) GVEC_INT => VAR%VEC_INT(1:RDIMS(1)) 

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_INT,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        If(.NOT. COLLECTDATA) GVEC_INT => VAR%VEC_INT(1:RDIMS(1))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GVEC_INT_temp(RDIMS(1)))
            CALL MPI_ALLREDUCE(GVEC_INT,GVEC_INT_temp,RDIMS(1),MPI_INTEGER,MPI_MAX,MPI_fvcom_group,ierr)
            GVEC_INT = GVEC_INT_temp
            DEALLOCATE(GVEC_INT_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GVEC_INT)
            GVEC_INT => VAR%VEC_INT(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1)
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_INT,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_INT,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_INT,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_INT,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GVEC_INT)
        

     END IF

     nullify(Gvec_INT)
     nullify(Lvec_INT)
     
 
!********************************************************************
! =====   ARRAY INTEGER DATA
!********************************************************************
  CASE(case_arr_int)

    IF(.NOT. ASSOCIATED(VAR%ARR_INT))THEN 

        IF(ASSOCIATED(VAR%CUB_INT))THEN
           IF(size(VAR%CUB_INT,1)==1) VAR%ARR_INT=>VAR%CUB_INT(1,1:,1:)
           IF(size(VAR%CUB_INT,2)==1) VAR%ARR_INT=>VAR%CUB_INT(1:,1,1:)
           IF(size(VAR%CUB_INT,3)==1) VAR%ARR_INT=>VAR%CUB_INT(1:,1:,1)
        ELSE        
           
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects ARR_INT data is NOT assocaited!")
        END IF
     END IF
     nsize=ubound(VAR%ARR_INT,1)

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)

        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE

        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)


        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           ALLOCATE(GARR_INT(RDIMS(1),RDIMS(2)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR: Allocate VEC_INT failed!")                   
 
           IF (UBOUND(VAR%ARR_INT,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
                & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE WRITE COUNT",& 
                & "varname: "//TRIM(var%VARNAME))

           LARR_INT => VAR%ARR_INT(1:nsize,1:RDIMS(2))
           
        ELSE

           IF (IOPROC) THEN
              IF (   UBOUND(VAR%ARR_INT,1) .LT. RDIMS(1) .or. &
                   & UBOUND(VAR%ARR_INT,2) .LT. RDIMS(2) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              GARR_INT => VAR%ARR_INT(1:RDIMS(1),1:RDIMS(2))
              GARR_INT = -HUGE(GARR_INT)
           ELSE

              IF (UBOUND(VAR%ARR_INT,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              LARR_INT => VAR%ARR_INT(1:nsize,1:RDIMS(2))
           END IF

        END IF

!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LARR_INT,GARR_INT)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LARR_INT,GARR_INT)
        
 
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN

        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GARR_INT => VAR%ARR_INT(1:RDIMS(1),1:RDIMS(2))
!        write(ipt,*)VAR%ARR_INT(1:DIM1,1:DIM2)

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_INT,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GARR_INT_temp(RDIMS(1),RDIMS(2)))
            CALL MPI_ALLREDUCE(GARR_INT,GARR_INT_temp,RDIMS(1)*RDIMS(2),MPI_INTEGER,MPI_MAX,MPI_fvcom_group,ierr)
            GARR_INT = GARR_INT_temp
            DEALLOCATE(GARR_INT_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GARR_INT)
            GARR_INT => VAR%ARR_INT(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1,1:RDIMS(2))
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_INT,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_INT,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_INT,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_INT,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GARR_INT)
        
     END IF

     nullify(GARR_INT)
     nullify(LARR_INT)


!********************************************************************
! =====   CUBE INTEGER DATA
!********************************************************************
  CASE(case_cub_int)

    IF(.NOT. ASSOCIATED(VAR%cub_INT))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects CUB_INT data is NOT assocaited!")
     END IF

     nsize=ubound(VAR%CUB_INT,1)

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)

        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)

        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           ALLOCATE(GCUB_INT(RDIMS(1),RDIMS(2),RDIMS(3)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR: Allocate CUB_INT failed!")                   

           IF (   UBOUND(VAR%CUB_INT,2) .LT. RDIMS(2) .or. &
                & UBOUND(VAR%CUB_INT,3) .LT. RDIMS(3) ) CALL FATAL_ERROR &
                & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                & "varname: "//TRIM(var%VARNAME))
           
           LCUB_INT => VAR%CUB_INT(1:nsize,1:RDIMS(2),1:RDIMS(3))
           
        ELSE

           IF (IOPROC) THEN
              IF (   UBOUND(VAR%CUB_INT,1) .LT. RDIMS(1) .or. &
                   & UBOUND(VAR%CUB_INT,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%CUB_INT,3) .LT. RDIMS(3) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              
              GCUB_INT => VAR%CUB_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))
              GCUB_INT = -HUGE(GCUB_INT)
           ELSE

              IF (   UBOUND(VAR%CUB_INT,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%CUB_INT,3) .LT. RDIMS(3) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              
              LCUB_INT => VAR%CUB_INT(1:nsize,1:RDIMS(2),1:RDIMS(3))
           END IF

        END IF
        
!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LCUB_INT,GCUB_INT)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LCUB_INT,GCUB_INT)
        
 
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN

       ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GCUB_INT => VAR%CUB_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_INT,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GCUB_INT_temp(RDIMS(1),RDIMS(2),RDIMS(3)))
            CALL MPI_ALLREDUCE(GCUB_INT,GCUB_INT_temp,RDIMS(1)*RDIMS(2)*RDIMS(3),MPI_INTEGER,MPI_MAX,MPI_fvcom_group,ierr)
            GCUB_INT = GCUB_INT_temp
            DEALLOCATE(GCUB_INT_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GCUB_INT)
            GCUB_INT => VAR%CUB_INT(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1,1:RDIMS(2),1:RDIMS(3))
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_INT,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_INT,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_INT,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_INT,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF        

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GCUB_INT)
        

     END IF

     
     nullify(GCUB_INT)
     nullify(LCUB_INT)

!********************************************************************
! =====   FOUR DIMENSION ARRAY INTEGER DATA
!********************************************************************
  CASE(case_fda_int)

    IF(.NOT. ASSOCIATED(VAR%fda_INT))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects FDA_INT data is NOT assocaited!")
     END IF

     nsize=ubound(VAR%FDA_INT,1)

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)

        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)

        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           ALLOCATE(GFDA_INT(RDIMS(1),RDIMS(2),RDIMS(3),RDIMS(4)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR: Allocate FDA_INT failed!")                   

           IF (   UBOUND(VAR%FDA_INT,2) .LT. RDIMS(2) .or. &
                & UBOUND(VAR%FDA_INT,3) .LT. RDIMS(3) .or. &
                & UBOUND(VAR%FDA_INT,3) .LT. RDIMS(4) ) CALL FATAL_ERROR &
                & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                & "varname: "//TRIM(var%VARNAME))
           
           LFDA_INT => VAR%FDA_INT(1:nsize,1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
           
        ELSE

           IF (IOPROC) THEN
              IF (   UBOUND(VAR%FDA_INT,1) .LT. RDIMS(1) .or. &
                   & UBOUND(VAR%FDA_INT,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%FDA_INT,3) .LT. RDIMS(3) .or. &
                   & UBOUND(VAR%FDA_INT,4) .LT. RDIMS(4) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              
              GFDA_INT => VAR%FDA_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
              GFDA_INT = -HUGE(GFDA_INT)
           ELSE

              IF (   UBOUND(VAR%FDA_INT,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%FDA_INT,3) .LT. RDIMS(3) .or. &
                   & UBOUND(VAR%FDA_INT,4) .LT. RDIMS(4) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              
              LFDA_INT => VAR%FDA_INT(1:nsize,1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
           END IF

        END IF
        
!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LFDA_INT,GFDA_INT)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LFDA_INT,GFDA_INT)
        
 
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN

       ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GFDA_INT => VAR%FDA_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_INT,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GFDA_INT_temp(RDIMS(1),RDIMS(2),RDIMS(3),RDIMS(4)))
            CALL MPI_ALLREDUCE(GFDA_INT,GFDA_INT_temp,RDIMS(1)*RDIMS(2)*RDIMS(3)*RDIMS(4),MPI_INTEGER,MPI_MAX,MPI_fvcom_group,ierr)
            GFDA_INT = GFDA_INT_temp
            DEALLOCATE(GFDA_INT_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GFDA_INT)
            GFDA_INT => VAR%FDA_INT(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1,1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_INT,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_INT,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_INT,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_INT,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GFDA_INT)
        

     END IF

     
     nullify(GFDA_INT)
     nullify(LFDA_INT)

!********************************************************************
! =====   SCALAR FLOATING POINT DATA
!********************************************************************
  CASE(case_scl_flt)

     IF(.NOT. ASSOCIATED(VAR%SCL_FLT))THEN 

        IF(ASSOCIATED(VAR%VEC_FLT))THEN
           IF(size(VAR%VEC_FLT)==1) VAR%SCL_FLT=>VAR%VEC_FLT(1)
        ELSE IF(ASSOCIATED(VAR%ARR_FLT))THEN
           IF(size(VAR%ARR_FLT)==1) VAR%SCL_FLT=>VAR%ARR_FLT(1,1)
        ELSE IF(ASSOCIATED(VAR%CUB_FLT))THEN
           IF(size(VAR%CUB_FLT)==1) VAR%SCL_FLT=>VAR%CUB_FLT(1,1,1)
        ELSE        

           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects SCL_FLT data is NOT assocaited!")
        END IF
     END IF

     SCL_FLT => VAR%SCL_FLT

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA .AND. USE_MPI_IO_MODE) THEN
# if defined(MULTIPROCESSOR)
        DEST = IOPROCID - 1
!!$        SOURCE = MSRID - 1
        SOURCE = 0
        NSIZE = 1
!!$        IF (MSR) CALL MPI_SEND(scl_flt,NSIZE,MPI_REAL,DEST,WVD_TAG&
!!$             &,MPI_FVCOM_GROUP,IERR)
!!$        IF (IOPROC) CALL MPI_RECV(scl_flt,NSIZE,MPI_REAL,SOURCE,WVD_TAG&
!!$             & ,MPI_FVCOM_GROUP,STAT,IERR)
        IF (MYID_iogroup == 1) CALL MPI_SEND(scl_flt,NSIZE,MPI_REAL,DEST,WVD_TAG&
             &,MPI_IO_GROUP,IERR)
        IF (IOPROC) CALL MPI_RECV(scl_flt,NSIZE,MPI_REAL,SOURCE,WVD_TAG&
             & ,MPI_IO_GROUP,STAT,IERR)
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

     IF (DUMP) THEN

        IF (SIZE(NSTART) .GT. 0) THEN

           if (product(nCOUNT) .NE. 1) CALL FATAL_ERROR&
                & ("NC_WRITE_VAR: NCOUNT dimension invalid while reading scl_flt?")

           ! ARGUMENT TO NF90_PUT_VAR MUST BE A VECTOR
           allocate(GVEC_FLT(1)); GVEC_FLT(1) = SCL_FLT
           IF(MYID_fgroup == 1 ) THEN        
             status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_FLT,NSTART,NCOUNT,NSTRIDE)
             CALL HANDLE_NCERR(status,trim(errmsg))
           ENDIF
           
           deallocate(GVEC_FLT)
        ELSE
           IF(MYID_fgroup == 1 ) THEN
             status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,SCL_FLT)
             CALL HANDLE_NCERR(status,trim(errmsg))
           ENDIF
        END IF
     END IF
     
     NULLIFY(SCL_FLT)

!********************************************************************
! =====   VECTOR FLOATING POINT DATA
!********************************************************************
  CASE(case_vec_flt)

    IF(.NOT. ASSOCIATED(VAR%VEC_FLT))THEN 
       
       IF(ASSOCIATED(VAR%ARR_FLT))THEN
          IF(size(VAR%ARR_FLT,1)==1) VAR%VEC_FLT=>VAR%ARR_FLT(1,1:)
          IF(size(VAR%ARR_FLT,2)==1) VAR%VEC_FLT=>VAR%ARR_FLT(1:,1)
       ELSE IF(ASSOCIATED(VAR%CUB_FLT))THEN
          IF(size(VAR%CUB_FLT,1)==1) THEN
             IF(size(VAR%CUB_FLT,2)==1) VAR%VEC_FLT=>VAR%CUB_FLT(1,1,1:)
             IF(size(VAR%CUB_FLT,3)==1) VAR%VEC_FLT=>VAR%CUB_FLT(1,1:,1)
          END IF
          IF(size(VAR%CUB_FLT,1)==2) THEN
             IF(size(VAR%CUB_FLT,3)==1) VAR%VEC_FLT=>VAR%CUB_FLT(1:,1,1)
          END IF
       ELSE        
          
          CALL PRINT_VAR(VAR)
          CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects VEC_FLT data is NOT assocaited!")
       END IF
    END IF

    nsize=ubound(VAR%VEC_FLT,1)

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)

        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
!!$        IF (IOPROC) LSIZES(IOPROCID) = 0
        IF (IOPROC) LSIZES(MYID) = 0
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)
 
        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           ALLOCATE(GVEC_FLT(RDIMS(1)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR:&
                & Allocate VEC_FLT failed!")                   
 
           LVEC_FLT => VAR%VEC_FLT
           
        ELSE

           IF (IOPROC) THEN
              IF (UBOUND(VAR%VEC_FLT,1) .LT. RDIMS(1)) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              GVEC_FLT => VAR%VEC_FLT(1:RDIMS(1))
              GVEC_FLT = -HUGE(GVEC_FLT)
           ELSE
              LVEC_FLT => VAR%VEC_FLT
           END IF

        END IF

!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LVEC_FLT,GVEC_FLT)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LVEC_FLT,GVEC_FLT)
        
 
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN

        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GVEC_FLT => VAR%VEC_FLT(1:RDIMS(1))


!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_FLT,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GVEC_FLT_temp(RDIMS(1)))
            CALL MPI_ALLREDUCE(GVEC_FLT,GVEC_FLT_temp,RDIMS(1),MPI_REAL,MPI_MAX,MPI_fvcom_group,ierr)
            GVEC_FLT = GVEC_FLT_temp
            DEALLOCATE(GVEC_FLT_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GVEC_FLT)
            GVEC_FLT => VAR%VEC_FLT(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1)
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_FLT,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_FLT,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_FLT,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_FLT,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GVEC_FLT)
        
     END IF

     nullify(Gvec_flt)
     nullify(Lvec_flt)
     

!********************************************************************
! =====   ARRAY FLOATING POINT DATA
!********************************************************************
  CASE(case_arr_flt)

    IF(.NOT. ASSOCIATED(VAR%ARR_FLT))THEN 

        IF(ASSOCIATED(VAR%CUB_FLT))THEN
           IF(size(VAR%CUB_FLT,1)==1) VAR%ARR_FLT=>VAR%CUB_FLT(1,1:,1:)
           IF(size(VAR%CUB_FLT,2)==1) VAR%ARR_FLT=>VAR%CUB_FLT(1:,1,1:)
           IF(size(VAR%CUB_FLT,3)==1) VAR%ARR_FLT=>VAR%CUB_FLT(1:,1:,1)
        ELSE        

           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects ARR_FLT data is NOT assocaited!")
        END IF
     END IF

     nsize=ubound(VAR%ARR_FLT,1)

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)
 
         ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)

        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           
           ALLOCATE(GARR_FLT(RDIMS(1),RDIMS(2)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR:&
                & Allocate VEC_FLT failed!")                   

           IF (UBOUND(VAR%ARR_FLT,2) .LT. RDIMS(2) ) CALL FATAL_ERROR &
                & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                & "varname: "//TRIM(var%VARNAME))
           
           LARR_FLT => VAR%ARR_FLT(1:nsize,1:RDIMS(2))
                         
        ELSE

           IF (IOPROC) THEN
              IF (   UBOUND(VAR%ARR_FLT,1) .LT. RDIMS(1) .or. &
                   & UBOUND(VAR%ARR_FLT,2) .LT. RDIMS(2) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              GARR_FLT => VAR%ARR_FLT(1:RDIMS(1),1:RDIMS(2))
!              GARR_FLT => VAR%ARR_FLT
              GARR_FLT = -HUGE(GARR_FLT)
           ELSE

              IF (UBOUND(VAR%ARR_FLT,2) .LT. RDIMS(2) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              LARR_FLT => VAR%ARR_FLT(1:nsize,1:RDIMS(2))
           END IF

        END IF
           

!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LARR_FLT,GARR_FLT)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LARR_FLT,GARR_FLT)

# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN

        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GARR_FLT => VAR%ARR_FLT(1:RDIMS(1),1:RDIMS(2))

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_FLT,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GARR_FLT_temp(RDIMS(1),RDIMS(2)))
            CALL MPI_ALLREDUCE(GARR_FLT,GARR_FLT_temp,RDIMS(1)*RDIMS(2),MPI_REAL,MPI_MAX,MPI_fvcom_group,ierr)
            GARR_FLT = GARR_FLT_temp
            DEALLOCATE(GARR_FLT_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GARR_FLT)
            GARR_FLT => VAR%ARR_FLT(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1,1:RDIMS(2))
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_FLT,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_FLT,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif 
#endif 
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_FLT,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_FLT,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF        

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GARR_FLT)
        

     END IF
     
     nullify(LARR_flt)
     nullify(GARR_flt)
     

!********************************************************************
! =====   CUBE FLOATING POINT DATA
!********************************************************************
  CASE(case_cub_flt)

    IF(.NOT. ASSOCIATED(VAR%cub_FLT))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects CUB_FLT data is NOT assocaited!")
     END IF

     nsize=ubound(VAR%CUB_FLT,1)

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)

        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)
        
        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           
           ALLOCATE(GCUB_FLT(RDIMS(1),RDIMS(2),RDIMS(3)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR: Allocate CUB_FLT failed!")                   

           IF (   UBOUND(VAR%CUB_FLT,2) .LT. RDIMS(2) .or. &
                & UBOUND(VAR%CUB_FLT,3) .LT. RDIMS(3) ) CALL FATAL_ERROR &
                & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                & "varname: "//TRIM(var%VARNAME))
           
           LCUB_FLT => VAR%CUB_FLT(1:nsize,1:RDIMS(2),1:RDIMS(3))
                         
        ELSE

           IF (IOPROC) THEN
              IF (   UBOUND(VAR%CUB_FLT,1) .LT. RDIMS(1) .or. &
                   & UBOUND(VAR%CUB_FLT,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%CUB_FLT,3) .LT. RDIMS(3) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              GCUB_FLT => VAR%CUB_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))
              GCUB_FLT = -HUGE(GCUB_FLT)
           ELSE

              IF (   UBOUND(VAR%CUB_FLT,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%CUB_FLT,3) .LT. RDIMS(3) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              LCUB_FLT => VAR%CUB_FLT(1:nsize,1:RDIMS(2),1:RDIMS(3))
           END IF

        END IF
        

!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LCUB_FLT,GCUB_FLT)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LCUB_FLT,GCUB_FLT)

# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN
        
        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GCUB_FLT => VAR%CUB_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_FLT,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GCUB_FLT_temp(RDIMS(1),RDIMS(2),RDIMS(3)))
            CALL MPI_ALLREDUCE(GCUB_FLT,GCUB_FLT_temp,RDIMS(1)*RDIMS(2)*RDIMS(3),MPI_REAL,MPI_MAX,MPI_fvcom_group,ierr)
            GCUB_FLT = GCUB_FLT_temp
            DEALLOCATE(GCUB_FLT_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GCUB_FLT)
            GCUB_FLT => VAR%CUB_FLT(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1,1:RDIMS(2),1:RDIMS(3))
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_FLT,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_FLT,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_FLT,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_FLT,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GCUB_FLT)
        
     END IF

     nullify(LCUB_flt)
     nullify(GCUB_flt)
     

!********************************************************************
! =====   FOUR DIMENSION ARRAY FLOATING POINT DATA
!********************************************************************
  CASE(case_fda_flt)

    IF(.NOT. ASSOCIATED(VAR%fda_FLT))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects FDA_FLT data is NOT assocaited!")
     END IF

     nsize=ubound(VAR%FDA_FLT,1)

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)

        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)
        
        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           
           ALLOCATE(GFDA_FLT(RDIMS(1),RDIMS(2),RDIMS(3),RDIMS(4)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR: Allocate FDA_FLT failed!")                   

           IF (   UBOUND(VAR%FDA_FLT,2) .LT. RDIMS(2) .or. &
                & UBOUND(VAR%FDA_FLT,3) .LT. RDIMS(3) .or. &
                & UBOUND(VAR%FDA_FLT,4) .LT. RDIMS(4) ) CALL FATAL_ERROR &
                & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                & "varname: "//TRIM(var%VARNAME))
           
           LFDA_FLT => VAR%FDA_FLT(1:nsize,1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
                         
        ELSE

           IF (IOPROC) THEN
              IF (   UBOUND(VAR%FDA_FLT,1) .LT. RDIMS(1) .or. &
                   & UBOUND(VAR%FDA_FLT,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%FDA_FLT,3) .LT. RDIMS(3) .or. &
                   & UBOUND(VAR%FDA_FLT,4) .LT. RDIMS(4) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              GFDA_FLT => VAR%FDA_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
              GFDA_FLT = -HUGE(GFDA_FLT)
           ELSE

              IF (   UBOUND(VAR%FDA_FLT,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%FDA_FLT,3) .LT. RDIMS(3) .or. &
                   & UBOUND(VAR%FDA_FLT,4) .LT. RDIMS(4) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              LFDA_FLT => VAR%FDA_FLT(1:nsize,1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
           END IF

        END IF
        

!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LFDA_FLT,GFDA_FLT)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LFDA_FLT,GFDA_FLT)

# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN
        
        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GFDA_FLT => VAR%FDA_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_FLT,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GFDA_FLT_temp(RDIMS(1),RDIMS(2),RDIMS(3),RDIMS(4)))
            CALL MPI_ALLREDUCE(GFDA_FLT,GFDA_FLT_temp,RDIMS(1)*RDIMS(2)*RDIMS(3)*RDIMS(4),MPI_REAL,MPI_MAX,MPI_fvcom_group,ierr)
            GFDA_FLT = GFDA_FLT_temp
            DEALLOCATE(GFDA_FLT_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GFDA_FLT)
            GFDA_FLT => VAR%FDA_FLT(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1,1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_FLT,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_FLT,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_FLT,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_FLT,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GFDA_FLT)
        
     END IF

     nullify(LFDA_flt)
     nullify(GFDA_flt)
     

!********************************************************************
! =====   SCALAR DOUBLE DATA
!********************************************************************
  CASE(case_scl_dbl)

     IF(.NOT. ASSOCIATED(VAR%SCL_DBL))THEN 

        IF(ASSOCIATED(VAR%VEC_DBL))THEN
           IF(size(VAR%VEC_DBL)==1) VAR%SCL_DBL=>VAR%VEC_DBL(1)
        ELSE IF(ASSOCIATED(VAR%ARR_DBL))THEN
           IF(size(VAR%ARR_DBL)==1) VAR%SCL_DBL=>VAR%ARR_DBL(1,1)
        ELSE IF(ASSOCIATED(VAR%CUB_DBL))THEN
           IF(size(VAR%CUB_DBL)==1) VAR%SCL_DBL=>VAR%CUB_DBL(1,1,1)
        ELSE        
           
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects SCL_DBL data is NOT assocaited!")
        END IF
     END IF

     SCL_DBL => VAR%SCL_DBL

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA .AND. USE_MPI_IO_MODE) THEN
# if defined(MULTIPROCESSOR)
        DEST = IOPROCID - 1
!!$        SOURCE = MSRID - 1
        SOURCE = 0
        NSIZE = 1
!!$        IF (MSR) CALL MPI_SEND(scl_dbl,NSIZE,MPI_DP,DEST,WVD_TAG&
!!$             &,MPI_FVCOM_GROUP,IERR)
!!$        IF (IOPROC) CALL MPI_RECV(scl_dbl,NSIZE,MPI_DP,SOURCE,WVD_TAG&
!!$             & ,MPI_FVCOM_GROUP,STAT,IERR)
        IF (MYID_iogroup == 1) CALL MPI_SEND(scl_dbl,NSIZE,MPI_DP,DEST,WVD_TAG&
             &,MPI_IO_GROUP,IERR)
        IF (IOPROC) CALL MPI_RECV(scl_dbl,NSIZE,MPI_DP,SOURCE,WVD_TAG&
             & ,MPI_IO_GROUP,STAT,IERR)
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

     IF (DUMP) THEN
        IF (SIZE(NSTART) .GT. 0) THEN

            if (product(nCOUNT) .NE. 1) CALL FATAL_ERROR&
                & ("NC_WRITE_VAR: NCOUNT dimension invalid while reading scl_dbl?")

           ! ARGUMENT TO NF90_PUT_VAR MUST BE A VECTOR
           allocate(GVEC_DBL(1)); GVEC_DBL(1) = SCL_DBL
           IF(MYID_fgroup == 1) THEN
             STATUS = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_DBL,NSTART,NCOUNT,NSTRIDE)
             CALL HANDLE_NCERR(status,trim(errmsg))
           ENDIF
           deallocate(GVEC_DBL)
        ELSE
           IF(MYID_fgroup == 1 ) THEN
             STATUS = NF90_PUT_VAR(VAR%NCID,VAR%VARID,SCL_DBL)
             CALL HANDLE_NCERR(status,trim(errmsg))
           ENDIF
        END IF
     END IF
     
     NULLIFY(SCL_DBL)

!********************************************************************
! =====   VECTOR DOUBLE DATA
!********************************************************************
  CASE(case_vec_dbl)


    IF(.NOT. ASSOCIATED(VAR%VEC_DBL))THEN 

       IF(ASSOCIATED(VAR%ARR_DBL))THEN
          IF(size(VAR%ARR_DBL,1)==1) VAR%VEC_DBL=>VAR%ARR_DBL(1,1:)
          IF(size(VAR%ARR_DBL,2)==1) VAR%VEC_DBL=>VAR%ARR_DBL(1:,1)
       ELSE IF(ASSOCIATED(VAR%CUB_DBL))THEN
          IF(size(VAR%CUB_DBL,1)==1) THEN
             IF(size(VAR%CUB_DBL,2)==1) VAR%VEC_DBL=>VAR%CUB_DBL(1,1,1:)
             IF(size(VAR%CUB_DBL,3)==1) VAR%VEC_DBL=>VAR%CUB_DBL(1,1:,1)
          END IF
          IF(size(VAR%CUB_DBL,1)==2) THEN
             IF(size(VAR%CUB_DBL,3)==1) VAR%VEC_DBL=>VAR%CUB_DBL(1:,1,1)
          END IF
       ELSE        
          
          CALL PRINT_VAR(VAR)
          CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects VEC_DBL data is NOT assocaited!")
       END IF
    END IF

    nsize=ubound(VAR%VEC_DBL,1)

     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)

        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)
 
        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           ALLOCATE(GVEC_DBL(RDIMS(1)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR:&
                & Allocate VEC_FLT failed!")                   
 
           LVEC_DBL => VAR%VEC_DBL
           
        ELSE

           IF (IOPROC) THEN
              IF (UBOUND(VAR%VEC_DBL,1) .LT. RDIMS(1)) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              GVEC_DBL => VAR%VEC_DBL(1:RDIMS(1))
              GVEC_DBL = -HUGE(GVEC_DBL)
           ELSE
              LVEC_DBL => VAR%VEC_DBL
           END IF

        END IF

!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LVEC_DBL,GVEC_DBL)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LVEC_DBL,GVEC_DBL)
        
 
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN


        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GVEC_DBL => VAR%VEC_DBL(1:RDIMS(1))

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_DBL,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GVEC_DBL_temp(RDIMS(1)))
            CALL MPI_ALLREDUCE(GVEC_DBL,GVEC_DBL_temp,RDIMS(1),MPI_DP,MPI_MAX,MPI_fvcom_group,ierr)
            GVEC_DBL = GVEC_DBL_temp
            DEALLOCATE(GVEC_DBL_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GVEC_DBL)
            GVEC_DBL => VAR%VEC_DBL(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1)
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_DBL,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_DBL,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_DBL,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GVEC_DBL,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF        

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GVEC_DBL)
        
     END IF

     nullify(Gvec_DBL)
     nullify(Lvec_DBL)
     

!********************************************************************
! =====   ARRAY DOUBLE DATA
!********************************************************************
  CASE(case_arr_dbl)

     IF(.NOT. ASSOCIATED(VAR%ARR_DBL))THEN 
        
        IF(ASSOCIATED(VAR%CUB_DBL))THEN
           IF(size(VAR%CUB_DBL,1)==1) VAR%ARR_DBL=>VAR%CUB_DBL(1,1:,1:)
           IF(size(VAR%CUB_DBL,2)==1) VAR%ARR_DBL=>VAR%CUB_DBL(1:,1,1:)
           IF(size(VAR%CUB_DBL,3)==1) VAR%ARR_DBL=>VAR%CUB_DBL(1:,1:,1)
        ELSE        
           
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects ARR_DBL data is NOT assocaited!")
        END IF
     END IF
     
     nsize=ubound(VAR%ARR_DBL,1)

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)

        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)
 
        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           
           ALLOCATE(GARR_DBL(RDIMS(1),RDIMS(2)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR:&
                & Allocate ARR_DBL failed!")                   

           IF (UBOUND(VAR%ARR_DBL,2) .LT. RDIMS(2) ) CALL FATAL_ERROR &
                & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                & "varname: "//TRIM(var%VARNAME))
           
           LARR_DBL => VAR%ARR_DBL(1:nsize,1:RDIMS(2))
                         
        ELSE

           IF (IOPROC) THEN
              IF (   UBOUND(VAR%ARR_DBL,1) .LT. RDIMS(1) .or. &
                   & UBOUND(VAR%ARR_DBL,2) .LT. RDIMS(2) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              GARR_DBL => VAR%ARR_DBL(1:RDIMS(1),1:RDIMS(2))
              GARR_DBL = -HUGE(GARR_DBL)
           ELSE

              IF (UBOUND(VAR%ARR_DBL,2) .LT. RDIMS(2) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))

              LARR_DBL => VAR%ARR_DBL(1:nsize,1:RDIMS(2))
           END IF

        END IF
           

!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LARR_DBL,GARR_DBL)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LARR_DBL,GARR_DBL)

# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN

        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GARR_DBL => VAR%ARR_DBL(1:RDIMS(1),1:RDIMS(2))

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_DBL,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GARR_DBL_temp(RDIMS(1),RDIMS(2)))
            CALL MPI_ALLREDUCE(GARR_DBL,GARR_DBL_temp,RDIMS(1)*RDIMS(2),MPI_DP,MPI_MAX,MPI_fvcom_group,ierr)
            GARR_DBL = GARR_DBL_temp
            DEALLOCATE(GARR_DBL_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GARR_DBL)
            GARR_DBL => VAR%ARR_DBL(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1,1:RDIMS(2))
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_DBL,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_DBL,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_DBL,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GARR_DBL,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF 

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GARR_DBL)
        
     END IF

     nullify(LARR_DBL)
     nullify(GARR_DBL)

!********************************************************************
! =====   CUBE DOUBLE DATA
!********************************************************************
  CASE(case_cub_dbl)

    IF(.NOT. ASSOCIATED(VAR%CUB_DBL))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects CUB_DBL data is NOT assocaited!")
     END IF

     nsize=ubound(VAR%CUB_DBL,1)

     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)
        
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)
 
        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           
           ALLOCATE(GCUB_DBL(RDIMS(1),RDIMS(2),RDIMS(3)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR:&
                & Allocate ARR_DBL failed!")                   
           
           IF (   UBOUND(VAR%CUB_DBL,2) .LT. RDIMS(2) .or. &
                & UBOUND(VAR%CUB_DBL,3) .LT. RDIMS(3) ) CALL FATAL_ERROR &
                & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                & "varname: "//TRIM(var%VARNAME))

           LCUB_DBL => VAR%CUB_DBL(1:nsize,1:RDIMS(2),1:RDIMS(3))
                         
        ELSE

           IF (IOPROC) THEN
              IF (   UBOUND(VAR%CUB_DBL,1) .LT. RDIMS(1) .or. &
                   & UBOUND(VAR%CUB_DBL,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%CUB_DBL,3) .LT. RDIMS(3) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              GCUB_DBL => VAR%CUB_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))
              GCUB_DBL = -HUGE(GCUB_DBL)
           ELSE

              IF (   UBOUND(VAR%CUB_DBL,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%CUB_DBL,3) .LT. RDIMS(3) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              
              LCUB_DBL => VAR%CUB_DBL(1:nSIZE,1:RDIMS(2),1:RDIMS(3))
           END IF

        END IF
           
!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LCUB_DBL,GCUB_DBL)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LCUB_DBL,GCUB_DBL)

# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN
        
        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GCUB_DBL => VAR%CUB_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_DBL,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GCUB_DBL_temp(RDIMS(1),RDIMS(2),RDIMS(3)))
            CALL MPI_ALLREDUCE(GCUB_DBL,GCUB_DBL_temp,RDIMS(1)*RDIMS(2)*RDIMS(3),MPI_DP,MPI_MAX,MPI_fvcom_group,ierr)
            GCUB_DBL = GCUB_DBL_temp
            DEALLOCATE(GCUB_DBL_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GCUB_DBL)
            GCUB_DBL => VAR%CUB_DBL(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1,1:RDIMS(2),1:RDIMS(3))
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_DBL,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_DBL,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_DBL,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GCUB_DBL,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF        

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GCUB_DBL)
        
     END IF

     nullify(LCUB_DBL)
     nullify(GCUB_DBL)
 
!********************************************************************
! =====   FOUR DIMENSION ARRAY DOUBLE DATA
!********************************************************************
  CASE(case_fda_dbl)

    IF(.NOT. ASSOCIATED(VAR%FDA_DBL))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects FDA_DBL data is NOT assocaited!")
     END IF

     nsize=ubound(VAR%FDA_DBL,1)

     IF(COLLECTDATA) THEN
# if defined(MULTIPROCESSOR)
        
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_WRITE_VAR INLCUDE THE IOPROC - MPI_COMM_FVCOM
        CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_COMM_FVCOM,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => Lsizes(1:nprocs)
        ! COLLECT OPERATION USE THE INTERNAL MAP - IT IS SMALLER
        GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPSize,FOUND)
        IF(.NOT. FOUND) CALL FATAL_ERROR &
             &("NC_WRITE_VAR: DON'T KNOW HOW TO COLLECT DATA WITHOUT A MAP THAT FITS!",&
             & "varname: "//TRIM(var%VARNAME))
        NULLIFY(NPSIZE)
        DEALLOCATE(LSIZES)
 
        IF (DUMP) THEN ! DUMP IS ONLY TRUE DURING COLLECT IF THE DATA
           ! IS TO BE DUMPED IMMEDIATLY AND SPACE NEEDS TO BE
           ! ALLOCATED TO STORE THE DATA ON A SINGLE PROC
           
           ALLOCATE(GFDA_DBL(RDIMS(1),RDIMS(2),RDIMS(3),RDIMS(4)),stat=status)
           IF (status /= 0 ) Call fatal_error("NC_WRITE_VAR:&
                & Allocate FDA_DBL failed!")                   
           
           IF (   UBOUND(VAR%FDA_DBL,2) .LT. RDIMS(2) .or. &
                & UBOUND(VAR%FDA_DBL,3) .LT. RDIMS(3) .or. &
                & UBOUND(VAR%FDA_DBL,4) .LT. RDIMS(4) ) CALL FATAL_ERROR &
                & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                & "varname: "//TRIM(var%VARNAME))

           LFDA_DBL => VAR%FDA_DBL(1:nsize,1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
                         
        ELSE

           IF (IOPROC) THEN
              IF (   UBOUND(VAR%FDA_DBL,1) .LT. RDIMS(1) .or. &
                   & UBOUND(VAR%FDA_DBL,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%FDA_DBL,3) .LT. RDIMS(3) .or. &
                   & UBOUND(VAR%FDA_DBL,4) .LT. RDIMS(4) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              GFDA_DBL => VAR%FDA_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
              GFDA_DBL = -HUGE(GFDA_DBL)
           ELSE

              IF (   UBOUND(VAR%FDA_DBL,2) .LT. RDIMS(2) .or. &
                   & UBOUND(VAR%FDA_DBL,3) .LT. RDIMS(3) .or. &
                   & UBOUND(VAR%FDA_DBL,4) .LT. RDIMS(4) ) CALL FATAL_ERROR &
                   & ("NC_WRITE_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH THE DIMENSIONS OF THE WRITE COUNT",& 
                   & "varname: "//TRIM(var%VARNAME))
              
              LFDA_DBL => VAR%FDA_DBL(1:nSIZE,1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
           END IF

        END IF
           
!!$        CALL PCOLLECT(MYID,COLLECTOR,NPROCS,GMAP,LFDA_DBL,GFDA_DBL)
        CALL PCOLLECT_IO(MYID,COLLECTOR,NPROCS_TOTAL,GMAP,LFDA_DBL,GFDA_DBL)

# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

 
     IF (DUMP) THEN
        
        ! IF YOU ARE THE IOPROC REASSIGN THE POINTER 
        If(.NOT. COLLECTDATA) GFDA_DBL => VAR%FDA_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))

!!$        status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_DBL,NSTART,NCOUNT,NSTRIDE)
!!$        CALL HANDLE_NCERR(status,trim(errmsg))
        IF ( USE_MPI_IO_MODE ) THEN
          IF(NPROCS_fvcom > 1) THEN
# if defined(MULTIPROCESSOR)
            ALLOCATE(GFDA_DBL_temp(RDIMS(1),RDIMS(2),RDIMS(3),RDIMS(4)))
            CALL MPI_ALLREDUCE(GFDA_DBL,GFDA_DBL_temp,RDIMS(1)*RDIMS(2)*RDIMS(3)*RDIMS(4),MPI_DP,MPI_MAX,MPI_fvcom_group,ierr)
            GFDA_DBL = GFDA_DBL_temp
            DEALLOCATE(GFDA_DBL_temp)
#if defined (PIO)
            allocate(NSTART_PIO(size(NSTART)))
            allocate(NCOUNT_PIO(size(NCOUNT)))
            allocate(NSTRIDE_PIO(size(NSTRIDE)))
            NSTART_PIO=NSTART
            NCOUNT_PIO=NCOUNT
            NSTRIDE_PIO=NSTRIDE
            PARTITION_PIO=NCOUNT_PIO(1)/NPROCS_fvcom
            rest_PIO= MOD(NCOUNT_PIO(1),NPROCS_fvcom)
            IF (MYID_fgroup <= rest_PIO ) THEN
              NSTART_PIO(1)=(PARTITION_PIO+1)*(MYID_fgroup-1)+1
              NCOUNT_PIO(1)=PARTITION_PIO+1
            ELSE
              NSTART_PIO(1)=PARTITION_PIO*(MYID_fgroup-1)+rest_PIO+1
              NCOUNT_PIO(1)=PARTITION_PIO
            ENDIF
            nullify(GFDA_DBL)
            GFDA_DBL => VAR%FDA_DBL(NSTART_PIO(1):NSTART_PIO(1)+NCOUNT_PIO(1)-1,1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_DBL,NSTART_PIO,NCOUNT_PIO,NSTRIDE_PIO)
            CALL HANDLE_NCERR(status,trim(errmsg))
            deallocate(NSTART_PIO)
            deallocate(NCOUNT_PIO)
            deallocate(NSTRIDE_PIO)
#else
            IF(MYID_fgroup == 1) THEN
              status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_DBL,NSTART,NCOUNT,NSTRIDE)
              CALL HANDLE_NCERR(status,trim(errmsg))
            ENDIF
#endif
#endif
          ELSE
            status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_DBL,NSTART,NCOUNT,NSTRIDE)
            CALL HANDLE_NCERR(status,trim(errmsg))
          ENDIF
        ELSE
          status = NF90_PUT_VAR(VAR%NCID,VAR%VARID,GFDA_DBL,NSTART,NCOUNT,NSTRIDE)
          CALL HANDLE_NCERR(status,trim(errmsg))
        ENDIF        

        ! IF YOU ARE COLLECTING AN DUMPING IN THE SAME PASS
        ! DEALLOCATE THE MEMORY YOU COLLECTED INTO
        IF (COLLECTDATA) deallocate(GFDA_DBL)
        
     END IF

     nullify(LFDA_DBL)
     nullify(GFDA_DBL)
 
!********************************************************************
! =====   SCALAR CHARACTER STRING DATA
!********************************************************************
  CASE(case_scl_chr)

     IF(.NOT. ASSOCIATED(VAR%SCL_CHR))THEN 

        IF (ASSOCIATED(VAR%vec_chr))THEN 
           IF(SIZE(VAR%vec_chr)==1) VAR%scl_chr => VAR%vec_chr(1)
        ELSE
           
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects SCL_CHR data is NOT assocaited!")
        END IF
     END IF
     
     SCL_CHR => VAR%SCL_CHR

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA .AND. USE_MPI_IO_MODE) THEN
# if defined(MULTIPROCESSOR)
        DEST = IOPROCID - 1
!!$        SOURCE = MSRID - 1
        SOURCE = 0
!!$        IF (MSR) THEN
        IF (MYID_iogroup == 1) THEN
           NSIZE = LEN(SCL_CHR)

           CALL MPI_SEND(NSIZE,1,MPI_INTEGER,DEST,WVD_TAG&
             &,MPI_IO_GROUP,IERR)

           CALL MPI_SEND(SCL_CHR,NSIZE,MPI_CHARACTER,DEST,WVD_TAG&
                &,MPI_IO_GROUP,IERR)
        END IF

        IF (IOPROC) THEN

           CALL MPI_RECV(NSIZE,1,MPI_INTEGER,SOURCE,WVD_TAG&
             & ,MPI_IO_GROUP,STAT,IERR)

           CALL MPI_RECV(scl_chr,NSIZE,MPI_CHARACTER,SOURCE,WVD_TAG&
             & ,MPI_IO_GROUP,STAT,IERR)

        END IF
# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

!!$     IF (DUMP) THEN
     IF (DUMP .and. MYID_fgroup == 1) THEN
        ! IF THE STRING IS SMALLER THAN THE DIMENSION
        ! AND NCOUNT IS TOO BIG
        NSIZE = LEN_TRIM(SCL_CHR)
        IF(NSIZE < RDIMS(1)) NCOUNT(1)=NSIZE

        STATUS = NF90_PUT_VAR(VAR%NCID,VAR%VARID,SCL_CHR,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))

     END IF
     
     NULLIFY(SCL_CHR)
!********************************************************************
! =====   VECTOR CHARACTER STRING DATA
!********************************************************************
  CASE(case_vec_chr)

     IF(.NOT. ASSOCIATED(VAR%VEC_CHR))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_WRITE_VAR: Variable objects VEC_CH&
             &R data is NOT assocaited!")
     END IF

     VEC_CHR => VAR%VEC_CHR

     ! ONLY COLLECT SCLs IF USING IOPROC: MSR SHOULD ALREADY HAVE THE DATA...
     IF(COLLECTDATA .AND. USE_MPI_IO_MODE) THEN
# if defined(MULTIPROCESSOR)
        DEST = IOPROCID - 1
!!$        SOURCE = MSRID - 1
        SOURCE = 0
        DO I = 1, size(Vec_chr)
           SCL_CHR => Vec_Chr(I)
           
!!$           IF (MSR) THEN
           IF (MYID_iogroup == 1) THEN
              NSIZE = LEN(SCL_CHR)
              
              CALL MPI_SEND(NSIZE,1,MPI_INTEGER,DEST,WVD_TAG&
                   &,MPI_IO_GROUP,IERR)
              
              CALL MPI_SEND(SCL_CHR,NSIZE,MPI_CHARACTER,DEST,WVD_TAG&
                   &,MPI_IO_GROUP,IERR)
           END IF
           
           IF (IOPROC) THEN
              
              CALL MPI_RECV(NSIZE,1,MPI_INTEGER,SOURCE,WVD_TAG&
                   & ,MPI_IO_GROUP,STAT,IERR)
              
              CALL MPI_RECV(scl_chr,NSIZE,MPI_CHARACTER,SOURCE,WVD_TAG&
                   & ,MPI_IO_GROUP,STAT,IERR)
              
           END IF

           Nullify(scl_chr)

        END DO

# else
        CALL FATAL_ERROR("NC_WRITE_VAR: COLLECT DATA FLAG IS T&
             &RUE BUT THIS IS NOT A MULTIPROCESSOR CASE?")
# endif
     END IF

     ! WRITING CHARACTER VECTORS IS A PROBLEM
     ! THE LENGTH HAS TO BE SET FOR EACH WRITE IN A DO LOOP
     ! POINT THE SCL_CHR TO EACH STRING AND DO A SEPERATE WRITE...
!!$     IF (DUMP) THEN
     IF (DUMP .and. MYID_fgroup == 1) THEN

        CNT=SIZE(NCOUNT)
        allocate(nSTRT(CNT),nCNT(CNT),nSTRD(CNT))
        nSTRT=NSTART
        nCNT=NCOUNT
        nSTRD=NSTRIDE
        
        
        DO I = 1,RDIMS(2)
           
           scl_chr => vec_chr(I)

           ! IF THE STRING IS SMALLER THAN THE DIMENSION
           ! AND NCOUNT IS TOO BIG 
           NSIZE = LEN_TRIM(SCL_CHR)
           IF(NSIZE < RDIMS(1))THEN
              nCNT(1)=NSIZE
           ELSE
              nCNT(1)= RDIMS(1)
           END IF
           
           ! WRITE THE I'th Entry (one entry)
           nSTRT(2)=I
           nCNT(2)=1

           STATUS = NF90_PUT_VAR(VAR%NCID,VAR%VARID,SCL_CHR,nSTRT,nCNT,nSTRD)
           CALL HANDLE_NCERR(status,trim(errmsg))

           nullify(scl_chr)

        END DO
        deallocate(nSTRT,nCNT,nSTRD)

     END IF
     
     NULLIFY(VEC_CHR)
     
  CASE default
     call print_var(VAR)
     call Fatal_error("NC_WRITE_VAR: UNKNOWN CASE")
  END SELECT

  IF(DUMP) THEN
     status = NF90_SYNC(VAR%ncid)
     CALL HANDLE_NCERR(status,trim(errmsg))
  END IF

  !Only deallocate if it is not pointing to Ncount
  IF(.not.ASSOCIATED(RDIMS,NCOUNT)) THEN
     IF(ASSOCIATED(RDIMS)) DEALLOCATE(RDIMS)
  ELSE
     NULLIFY(RDIMS)
  END IF


  IF(PRESENT(IOSTART)) THEN
     nullify(nstart)
  ELSE
     deallocate(nstart)
  END IF

  IF(PRESENT(IOCOUNT)) THEN
     nullify(nCOUNT)
  ELSE
     deallocate(nCOUNT)
  END IF

  IF(PRESENT(IOSTRIDE)) THEN
     nullify(nstride)
  ELSE
     deallocate(nstride)
  END IF



  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END NC_WRITE_VAR:"

END SUBROUTINE NC_WRITE_VAR
!====================================================================
!====================================================================
! Watch out - pointing to part of an array is tricky. The index may be reset
SUBROUTINE NC_READ_VAR(VAR,STKCNT,STKRNG,IOSTART,IOCOUNT,IOSTRIDE,DEALERID,PARALLEL)
  USE CONTROL
  implicit none
  TYPE(NCVAR), POINTER :: VAR
  INTEGER, INTENT(IN), OPTIONAL :: STKCNT
  INTEGER, INTENT(IN), OPTIONAL :: STKRNG(2)
  INTEGER, INTENT(IN), OPTIONAL :: DEALERID
  LOGICAL, INTENT(IN), OPTIONAL :: PARALLEL
  INTEGER, ALLOCATABLE,TARGET, OPTIONAL :: IOSTART(:), IOCOUNT(:), IOSTRIDE(:)

  ! LOCAL CONTROL VARIALBES
  INTEGER :: DEALER
  LOGICAL :: PAR_READ
  LOGICAL :: SER_READ

  ! NF90 SUBSET VARIABLES
  INTEGER, POINTER :: NSTART(:), NCOUNT(:), NSTRIDE(:)

  ! The size of the data returned by nf90_get/put_var
  INTEGER, POINTER :: RDIMS(:)
  ! For use in put/get_var for vector character data
  INTEGER, POINTER :: NSTRT(:), NCNT(:), NSTRD(:)

  INTEGER :: CODE
  INTEGER :: XTYPE
  INTEGER :: CNT,DIMCNT
  INTEGER :: NSIZE, MYSIZE
!  integer :: dim1,dim2,dim3,dim4

  ! POINTERS AND SUCH
  TYPE(NCDIM), POINTER :: DIM
  TYPE(NCDIMP), POINTER :: DIMLINK
  TYPE(NCFILE), POINTER :: NCF
  LOGICAL FOUND

# if defined(MULTIPROCESSOR)
  !MPI COMM STUFF
  INTEGER, PARAMETER :: WVD_TAG = 40002
  INTEGER :: DEST, SOURCE, IERR
!  INTEGER :: STAT(MPI_STATUS_SIZE)
  TYPE(MAP), POINTER :: GMAP(:)
  INTEGER, POINTER :: Lsizes(:), NPsize(:)
# endif

  INTEGER, PARAMETER :: case_scl_int = 1
  INTEGER, PARAMETER :: case_vec_int = 2
  INTEGER, PARAMETER :: case_arr_int = 3
  INTEGER, PARAMETER :: case_cub_int = 4
  INTEGER, PARAMETER :: case_fda_int = 5

  INTEGER, PARAMETER :: case_scl_flt = 6
  INTEGER, PARAMETER :: case_vec_flt = 7
  INTEGER, PARAMETER :: case_arr_flt = 8
  INTEGER, PARAMETER :: case_cub_flt = 9
  INTEGER, PARAMETER :: case_fda_flt = 10

  INTEGER, PARAMETER :: case_scl_dbl = 11
  INTEGER, PARAMETER :: case_vec_dbl = 12
  INTEGER, PARAMETER :: case_arr_dbl = 13
  INTEGER, PARAMETER :: case_cub_dbl = 14
  INTEGER, PARAMETER :: case_fda_dbl = 15

  INTEGER, PARAMETER :: case_scl_chr = 16
  INTEGER, PARAMETER :: case_vec_chr = 17

  ! TEMPORARY STORAGE FOR DATA IF COLLECTED TO MASTER PROC
  INTEGER, POINTER                :: SCL_INT
  INTEGER, POINTER,DIMENSION(:)   :: GVEC_INT
  INTEGER, POINTER,DIMENSION(:,:) :: GARR_INT
  INTEGER, POINTER,DIMENSION(:,:,:) :: GCUB_INT
  INTEGER, POINTER,DIMENSION(:,:,:,:) :: GFDA_INT

  INTEGER, POINTER,DIMENSION(:)   :: LVEC_INT
  INTEGER, POINTER,DIMENSION(:,:) :: LARR_INT
  INTEGER, POINTER,DIMENSION(:,:,:) :: LCUB_INT
  INTEGER, POINTER,DIMENSION(:,:,:,:) :: LFDA_INT

  REAL(SPA), POINTER                :: SCL_FLT
  REAL(SPA), POINTER,DIMENSION(:)   :: LVEC_FLT
  REAL(SPA), POINTER,DIMENSION(:,:) :: LARR_FLT
  REAL(SPA), POINTER,DIMENSION(:,:,:) :: LCUB_FLT
  REAL(SPA), POINTER,DIMENSION(:,:,:,:) :: LFDA_FLT

  REAL(SPA), POINTER,DIMENSION(:)   :: GVEC_FLT
  REAL(SPA), POINTER,DIMENSION(:,:) :: GARR_FLT
  REAL(SPA), POINTER,DIMENSION(:,:,:) :: GCUB_FLT
  REAL(SPA), POINTER,DIMENSION(:,:,:,:) :: GFDA_FLT
  
  REAL(DP), POINTER                :: SCL_DBL
  REAL(DP), POINTER,DIMENSION(:)   :: GVEC_DBL
  REAL(DP), POINTER,DIMENSION(:,:) :: GARR_DBL
  REAL(DP), POINTER,DIMENSION(:,:,:) :: GCUB_DBL
  REAL(DP), POINTER,DIMENSION(:,:,:,:) :: GFDA_DBL

  REAL(DP), POINTER,DIMENSION(:)   :: LVEC_DBL
  REAL(DP), POINTER,DIMENSION(:,:) :: LARR_DBL
  REAL(DP), POINTER,DIMENSION(:,:,:) :: LCUB_DBL
  REAL(DP), POINTER,DIMENSION(:,:,:,:) :: LFDA_DBL

  CHARACTER(LEN=80), POINTER        :: SCL_CHR
  CHARACTER(LEN=80), POINTER        :: VEC_CHR(:)
  INTEGER :: nlen

  CHARACTER(len=3) :: char1,char2,char3
  ! DATA FOR PUT VAR COMMANDS:
  INTEGER :: STATUS, I
  CHARACTER(LEN=120)          :: errmsg

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START NC_READ_VAR:"

  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("NC_READ_VAR: Variable object argument is not assocaited!")

  ! INITIALIZE SOME MEMORY
  nlen=0
  STATUS=0
  I=0
  CODE=0
  XTYPE=0
  CNT=0
  DIMCNT=0
  NSIZE=0
  MYSIZE=0
  DEALER=0
  PAR_READ=.false.
  SER_READ=.false.
  FOUND=.FALSE.


# if defined(MULTIPROCESSOR)
  NULLIFY(GMAP)
  IERR=0
  SOURCE=0
  DEST=0
# endif

  !NULLIFY POINTERS
  NULLIFY(NSTART,NCOUNT,NSTRIDE)
  NULLIFY(RDIMS,NSTRT,NCNT,NSTRD)
  NULLIFY(DIM,DIMLINK,NCF)
  NULLIFY(SCL_INT,GVEC_INT,GARR_INT,GCUB_INT,GFDA_INT,LVEC_INT,LARR_INT,LCUB_INT,LFDA_INT)

  NULLIFY(SCL_FLT,GVEC_FLT,GARR_FLT,GCUB_FLT,GFDA_FLT,LVEC_FLT,LARR_FLT,LCUB_FLT,LFDA_FLT)
  NULLIFY(SCL_DBL,GVEC_DBL,GARR_DBL,GCUB_DBL,GFDA_DBL,LVEC_DBL,LARR_DBL,LCUB_DBL,LFDA_DBL)
  NULLIFY(SCL_CHR,VEC_CHR)

  ! CHECK TO SEE IF A DEALER ID WAS SPECIFIED?
  IF (PRESENT(DEALERID)) THEN
     DEALER = DEALERID
  ELSE
     DEALER = MSRID ! FOR SERIAL CASE DEALER MUST BE MASTER, 
  END IF


! THIS IS HARD TO DO - YOU CAN'T RELIABLE GET THE NCF POINTER?
! NOT ALL FILES ARE ADDED TO THE FILEHEAD BEFORE CALLING READ/WRITE
!!$  ! MAKE SURE THE FILE IS OPEN 
!!$  IF(DEALER == MYID) THEN
!!$     IF(.NOT. Associated(VAR%NCID)) THEN
!!$        Call Print_var(VAR)
!!$        Call Fatal_error&
!!$          &("NC_READ_VAR: VARIABLE NCID NOT ASSOCIATED?")
!!$     END IF
!!$
!!$     NCF => FIND_FILE_BYNCID(FILEHEAD ,VAR%NCID,FOUND)
!!$     IF (.NOT. FOUND)THEN
!!$        CALL PRINT_FILE_LIST(FILEHEAD)
!!$        Call Print_var(VAR)
!!$        Call Fatal_error&
!!$          &(" NC_READ_VAR: FILE NOT FOUND IN FILE THE FILEHEAD?")
!!$     END IF
!!$
!!$     IF(.NOT. NCF%OPEN) CALL NC_OPEN(NCF)
!!$
!!$     IF(.NOT. NCF%CONNECTED) THEN
!!$        CALL PRINT_FILE(NCF)
!!$        CALL FATAL_ERROR("NC_READ_VAR: Attempt to read variable from file which is not connected?",& 
!!$             & "You must call 'nc_save' or 'nc_load' before trying to read data!")
!!$     END IF
!!$
!!$  END IF


  ! CHECK TO SEE IF GLOBAL PAR/SERIAL IS OVER-RIDDEN FOR THIS READ?
  IF (PRESENT(PARALLEL)) THEN
     IF(SERIAL .and. PARALLEL)THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: PARALLEL ERROR!",&
             & "The model is running on a single processor,",&
             & "but asked to do a parallel NC_READ_VAR operation!")
     END IF
     PAR_READ= PARALLEL ! SET THE LOCAL PASSED OPTION
  ELSE
     PAR_READ= PAR ! SET THE GLOBAL
  END IF
  SER_READ=.NOT. PAR_READ

  ! COUNT THE NUMBER OF DIMENSIONS IN THE VARIABLE
  DIMCNT = count_dim_list(VAR)

  IF(DBG_SET(DBG_SBRIO)) THEN

     write(char2,'(I3.3)')DEALER
     write(char3,'(I3.3)')myid
     WRITE(IPT,*)"NC_READ_VAR Arguments:"
     call print_var(var)
     WRITE(IPT,*)"; DEALER="//char2//"; MYID="//char3//"; PAR_READ=",PAR_READ
     
     
     IF(PRESENT(STKCNT)) THEN
        WRITE(IPT,*) "STKCNT=",STKCNT
     ELSE
        WRITE(IPT,*) "STKCNT= NONE"
     END IF

     IF(PRESENT(STKRNG)) THEN
        WRITE(IPT,*) "STKRNG=",STKRNG
     ELSE
        WRITE(IPT,*) "STKRNG= NONE"
     END IF
     
     IF(PRESENT(IOSTART)) THEN
        WRITE(IPT,*) "IOSTART=",IOSTART
     ELSE
        WRITE(IPT,*) "IOSTART= NONE"
     END IF
     
     IF(PRESENT(IOCOUNT)) THEN
        WRITE(IPT,*) "IOCOUNT=",IOCOUNT
     ELSE
        WRITE(IPT,*) "IOCOUNT= NONE"
     END IF
     
     IF(PRESENT(IOSTRIDE)) THEN
        WRITE(IPT,*) "IOSTRIDE=",IOSTRIDE
     ELSE
        WRITE(IPT,*) "IOSTRIDE= NONE"
     END IF
     
  END IF

  IF(VAR%NCID == -1 .and. DEALER == MYID) THEN
     CALL PRINT_VAR(VAR)
     CALL FATAL_ERROR("NC_READ_VAR: CAN NOT READ FILE, IT IS NOT OPEN!")
  END IF
  


  IF ( PRESENT(STKCNT) ) THEN

     IF (PRESENT(STKRNG) .or. PRESENT(IOSTART) .or. PRESENT(IOCOUNT) .or. PRESENT(IOSTRIDE))THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: You can not pass both STKCNT and STKRNG/START/COUNT/STRIDE !",&
             &"Set STKCNT to read a time slice filling all other dimensions. OR",&
             &"Set IOSTART/IOCOUNT/(IOSTRIDE) to read a specific range.")
     END IF

     DIM => FIND_UNLIMITED(VAR,FOUND)
     IF(.NOT.FOUND) THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR&
             &("NC_READ_VAR: CALLED WITH STKCNT ARGUMENT BUT VARIABLE IS NOT UNLIMITED?")
     END IF
     NULLIFY(DIM)


     ALLOCATE(NSTART(DIMCNT),NCOUNT(DIMCNT),NSTRIDE(DIMCNT))


     ! SET THE VARIABLES CURRENT STACK COUNT
     VAR%CURR_STKCNT = STKCNT
     
     ! SET THE NF90_PUT_VAR DIMENSIONS
     NSTART=1 ! START AT ONE, EXCEPT FOR THE TIME VARIABLE
     NSTART(DIMCNT) = STKCNT

     
     DIMLINK => VAR%DIMS%NEXT
     DO I = 1,DIMCNT ! GET THE OUTPUT VARIABLE DIMENSIONS
        DIM => DIMLINK%DIM  
        NCOUNT(I)= DIM%DIM
        DIMLINK => DIMLINK%NEXT
     END DO
     NCOUNT(DIMCNT)=1 ! SET THE TIME OUTPUT DIMENSION TO 1
     
     NSTRIDE=1 ! ALWAYS USE STRIDE 1 FOR STKCNT INPUT

     ! FOR TIME DEPENDANT DATA THE RANK OF THE ALLOCATED MEMORY IS
     ! ONE LESS THAN THE RANK OF THE FILE's VARIALBE!
     DIMCNT = DIMCNT -1

     ! THE DIMENSIONS READ WILL BE THE VALUES OF NCOUNT, NOT
     ! INCLUDING TIME
     IF(DIMCNT > 0)THEN
        ALLOCATE(RDIMS(DIMCNT))
        RDIMS(1:DIMCNT)=NCOUNT(1:DIMCNT)
     ELSE IF (DIMCNT == 0) THEN
         ALLOCATE(RDIMS(1))
        RDIMS(1)=NCOUNT(1)
     ELSE
        nullify(RDIMS)
     END IF

  ELSEIF ( PRESENT(STKRNG) ) THEN

     IF (PRESENT(STKCNT) .or. PRESENT(IOSTART) .or. PRESENT(IOCOUNT) .or. PRESENT(IOSTRIDE))THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: You can not pass both STKRNG and STKCNT/START/COUNT/STRIDE !",&
             &"Set STKRNG to read a time range filling all other dimensions. OR",&
             &"Set IOSTART/IOCOUNT/(IOSTRIDE) to read a specific range.")
     END IF

     DIM => FIND_UNLIMITED(VAR,FOUND)
     IF(.NOT.FOUND) THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR&
             &("NC_READ_VAR: CALLED WITH STKRNG ARGUMENT BUT VARIABLE IS NOT UNLIMITED?")
     END IF
     NULLIFY(DIM)


     ALLOCATE(NSTART(DIMCNT),NCOUNT(DIMCNT),NSTRIDE(DIMCNT))


     ! SET THE VARIABLES CURRENT STACK COUNT
     VAR%CURR_STKCNT = -1
     
     ! SET THE NF90_PUT_VAR DIMENSIONS
     NSTART=1 ! START AT ONE, EXCEPT FOR THE TIME VARIABLE
     NSTART(DIMCNT) = STKRNG(1)

     
     DIMLINK => VAR%DIMS%NEXT
     DO I = 1,DIMCNT ! GET THE OUTPUT VARIABLE DIMENSIONS
        DIM => DIMLINK%DIM  
        NCOUNT(I)= DIM%DIM
        DIMLINK => DIMLINK%NEXT
     END DO
     NCOUNT(DIMCNT)=STKRNG(2)-STKRNG(1)+1 ! SET THE TIME RANGE
     
     NSTRIDE=1 ! ALWAYS USE STRIDE 1 FOR STKRNG INPUT

     ! FOR A TIME RANGE THE RANK OF THE ALLOCATED MEMORY IS
     ! THE SAME AS THE RANK OF THE FILE's VARIALBE!
!     DIMCNT = DIMCNT

     ! THE DIMENSIONS READ WILL BE THE VALUES OF NCOUNT, NOT
     ! INCLUDING TIME
     RDIMS=>NCOUNT

  ELSE IF( PRESENT(IOSTART) .and. PRESENT(IOCOUNT)) THEN

     NSTART=>IOSTART
     NCOUNT=>IOCOUNT

     IF(.not. PRESENT(IOSTRIDE)) THEN
        ALLOCATE(NSTRIDE(DIMCNT))
        NSTRIDE=1
     ELSE
        NSTRIDE=>IOSTRIDE
     END IF

     IF(DIMCNT /= size(NSTART) .or. &
          & DIMCNT /= size(NCOUNT) .or. &
          & DIMCNT /= size(NSTRIDE) ) THEN
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR&
        & ("NC_READ_VAR: Variable's # of file dimensions does not matach size(NSTART/NCOUNT/NSTRIDE) arugments?")
     END IF

     ! SET THE VARIABLES CURRENT STACK COUNT: not defined for this
     ! kind of read/write
     VAR%CURR_STKCNT = -1

     ! ONLY COUNT THE NONE SINGLETON DIMENSIONS OF A VARIABLE. 
     CNT = 0
     DO I = 1,DIMCNT
        IF(NCOUNT(I)>1) CNT=CNT+1
     END DO

     ! NOW RECORD THE DIMENSIONS OF THE DATA THAT WILL BE READ INTO MEMORY
     IF (CNT > 0) THEN
        
        ALLOCATE(RDIMS(CNT))
        CNT = 0
        DO I = 1,DIMCNT
           IF(NCOUNT(I)>1) THEN
              CNT=CNT+1
              RDIMS(CNT)=NCOUNT(I)
           END IF
        END DO

     ELSE
        ALLOCATE(RDIMS(1)) 
        RDIMS(1) = 1
     END IF

     
     ! NOW SET THE DIMENSION OF THE DATA VARIABLE IN MEMORY
     DIMCNT=CNT

  ELSE IF( .not. (PRESENT(IOSTART) .or. PRESENT(IOCOUNT) .or.&
       & PRESENT(STKCNT) .or. PRESENT(IOSTRIDE))) THEN

     ALLOCATE(NSTART(DIMCNT),NCOUNT(DIMCNT),NSTRIDE(DIMCNT))

     ! SET THE VARIABLES CURRENT STACK COUNT
     VAR%CURR_STKCNT = -1
     
     ! SET THE NF90_PUT_VAR DIMENSIONS
     NSTART=1 ! START AT ONE, EXCEPT FOR THE TIME VARIABLE
     
     DIMLINK => VAR%DIMS%NEXT
     DO I = 1,DIMCNT ! GET THE OUTPUT VARIABLE DIMENSIONS
        DIM => DIMLINK%DIM  
        NCOUNT(I)= DIM%DIM
        DIMLINK => DIMLINK%NEXT
     END DO
     
     NSTRIDE=1 ! ALWAYS USE STRIDE 1 IF NO ARGUMENTS ARE PASSED

     ! ONLY COUNT THE NONE SINGLETON DIMENSIONS OF A VARIABLE. 
     CNT = 0
     DO I = 1,DIMCNT
        IF(NCOUNT(I)>1) CNT=CNT+1
     END DO

     ! NOW RECORD THE DIMENSIONS OF THE DATA THAT WILL BE READ INTO MEMORY
     IF (CNT > 0) THEN
        
        ALLOCATE(RDIMS(CNT))
        CNT = 0
        DO I = 1,DIMCNT
           IF(NCOUNT(I)>1) THEN
              CNT=CNT+1
              RDIMS(CNT)=NCOUNT(I)
           END IF
        END DO

     ELSE
        ALLOCATE(RDIMS(1)) 
        RDIMS(1) = 1
     END IF

     
     ! NOW SET THE DIMENSION OF THE DATA VARIABLE IN MEMORY
     DIMCNT=CNT

  ELSE

     IF(DBG_SET(DBG_LOG)) THEN
        write(ipt,*) "# IOSTART  ::",PRESENT(IOSTART)
        write(ipt,*) "# IOCOUNT  ::",PRESENT(IOCOUNT)
        write(ipt,*) "# IOSTRIDE ::",PRESENT(IOSTRIDE)
        write(ipt,*) "# STKCNT   ::",PRESENT(STKCNT)
        write(ipt,*) "# STKRNG   ::",PRESENT(STKRNG)
     END IF
     
     CALL FATAL_ERROR("NC_READ_VAR: YOU SPECIFIED AN ILLEGAL COMBINATION OF AGUMENTS?",&
          & "Valid choices are STKCNT or STKRNG or NSTART,NCOUNT,(NSTRIDE) or none")
  END IF

  IF(DBG_SET(DBG_SBRIO)) THEN
     write(IPT,*) "FILE DIMENSION COUNT IS   ::", count_dim_list(VAR)
     write(IPT,*) "MEMORY DIMENSION COUNT IS ::",DIMCNT
  END IF



  !DETERMIN WHICH CASE WE ARE WRITING DATA FOR
  code = -1

  select case(VAR%XTYPE)
  case(NF90_BYTE)
     call Fatal_error("NC_READ_VAR: NOT SET UP TO WRITE BYTE DATA")
  case(NF90_SHORT)
     call Fatal_error("NC_READ_VAR: NOT SET UP TO WRITE SHORT DATA")
  case(NF90_INT)
     if (DIMCNT == 0) CODE = case_scl_int
     if (DIMCNT == 1) CODE = case_vec_int
     if (DIMCNT == 2) CODE = case_arr_int
     if (DIMCNT == 3) CODE = case_cub_int
     if (DIMCNT == 4) CODE = case_fda_int
  case(NF90_FLOAT)
!--Single Precision Coding------------------------------------------------------!
# if !defined(DOUBLE_PRECISION)
     ! READ SINGLE VARIABLES INTO SINGLE MEMORY POINTERS
     if (DIMCNT == 0) CODE = case_scl_flt
     if (DIMCNT == 1) CODE = case_vec_flt
     if (DIMCNT == 2) CODE = case_arr_flt
     if (DIMCNT == 3) CODE = case_cub_flt
     if (DIMCNT == 4) CODE = case_fda_flt
#  else 
!--Double Precision Coding------------------------------------------------------!
     ! READ SINGLE VARIABLES INTO DOUBLE MEMORY POINTERS
     if (DIMCNT == 0) CODE = case_scl_dbl
     if (DIMCNT == 1) CODE = case_vec_dbl
     if (DIMCNT == 2) CODE = case_arr_dbl
     if (DIMCNT == 3) CODE = case_cub_dbl
     if (DIMCNT == 4) CODE = case_fda_dbl
# endif

  case(NF90_DOUBLE)
     ! ALWAYS READ DOUBLE VARIABLES INTO DOUBLE MEMORY POINTERS
     if (DIMCNT == 0) CODE = case_scl_dbl
     if (DIMCNT == 1) CODE = case_vec_dbl
     if (DIMCNT == 2) CODE = case_arr_dbl
     if (DIMCNT == 3) CODE = case_cub_dbl
     if (DIMCNT == 4) CODE = case_fda_dbl
  case(NF90_CHAR)

     IF(NCOUNT(1) == 1) THEN
        WRITE(IPT,*) "SINGLETON CHARACTER DATA!"
        IF(.not. ASSOCIATED(RDIMS,NCOUNT)) THEN
           DEALLOCATE(RDIMS)
        ELSE
           NULLIFY(RDIMS)
        END IF

        DIMCNT = DIMCNT+1

        ALLOCATE(RDIMS(DIMCNT))
        CNT = 1
        RDIMS(1) = NCOUNT(1)
        DO I = 2,size(ncount)
           IF(NCOUNT(I)>1) THEN
              CNT=CNT+1
              RDIMS(CNT)=NCOUNT(I)
           END IF
        END DO
     END IF


     if (DIMCNT == 1) CODE = case_scl_chr
     if (DIMCNT == 2) CODE = case_vec_chr
     ! First dim is length of string
     ! Second dim is time
  case default
     call Fatal_error("NC_READ_VAR: Unkown data type?")
  end select

  ! BASED ON CODE WRITE THE DATA
  errmsg="NC_READ_VAR: VARIABLE: "//VAR%varname//"; Can not be read by nf90_get_var!"

  SELECT CASE(CODE)
!********************************************************************
! =====   SCALAR FLOATING POINT DATA
!********************************************************************
  CASE(case_scl_FLT)
     
     IF(.NOT. ASSOCIATED(VAR%SCL_FLT))THEN 

        IF(ASSOCIATED(VAR%VEC_FLT))THEN
           IF(size(VAR%VEC_FLT)==1) VAR%SCL_FLT=>VAR%VEC_FLT(1)
        ELSE IF(ASSOCIATED(VAR%ARR_FLT))THEN
           IF(size(VAR%ARR_FLT)==1) VAR%SCL_FLT=>VAR%ARR_FLT(1,1)
        ELSE IF(ASSOCIATED(VAR%CUB_FLT))THEN
           IF(size(VAR%CUB_FLT)==1) VAR%SCL_FLT=>VAR%CUB_FLT(1,1,1)
        ELSE IF(ASSOCIATED(VAR%FDA_FLT))THEN
           IF(size(VAR%FDA_FLT)==1) VAR%SCL_FLT=>VAR%FDA_FLT(1,1,1,1)
        ELSE        
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_READ_VAR: Variable objects SCL_FLT data is NOT assocaited!")
        END IF

     END IF
  
     IF (SER_READ .OR. DEALER .EQ. MYID) THEN
        
        
        IF (SIZE(NSTART) .GT. 0) THEN
           
           if (product(nCOUNT) .NE. 1) CALL FATAL_ERROR&
                & ("NC_READ_VAR: NCOUNT dimension invalid while reading scl_flt?")
           
           allocate(gvec_FLT(1))
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GVEC_FLT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
           VAR%SCL_FLT = gvec_flt(1)
           deallocate(gvec_flt)
        ELSE
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,VAR%SCL_FLT)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF
     END IF
     
     IF (PAR_READ) THEN
# if defined (MULTIPROCESSOR)
        
        SOURCE = DEALER -1 
        CALL MPI_BCAST(VAR%SCL_FLT,1,MPI_REAL,SOURCE&
             &,MPI_FVCOM_GROUP,IERR)
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
       
     END IF     

     nullify(gvec_flt)

!********************************************************************
! =====   VECTOR FLOATING POINT DATA
!********************************************************************
  CASE(case_vec_flt)

    IF(.NOT. ASSOCIATED(VAR%VEC_FLT))THEN 

        IF(ASSOCIATED(VAR%ARR_FLT))THEN
           IF(size(VAR%ARR_FLT,1)==1) VAR%VEC_FLT=>VAR%ARR_FLT(1,1:)
           IF(size(VAR%ARR_FLT,2)==1) VAR%VEC_FLT=>VAR%ARR_FLT(1:,1)
        ELSE IF(ASSOCIATED(VAR%CUB_FLT))THEN
           IF(size(VAR%CUB_FLT,1)==1) THEN
              IF(size(VAR%CUB_FLT,2)==1) VAR%VEC_FLT=>VAR%CUB_FLT(1,1,1:)
              IF(size(VAR%CUB_FLT,3)==1) VAR%VEC_FLT=>VAR%CUB_FLT(1,1:,1)
           END IF
           IF(size(VAR%CUB_FLT,1)==2) THEN
              IF(size(VAR%CUB_FLT,3)==1) VAR%VEC_FLT=>VAR%CUB_FLT(1:,1,1)
           END IF           
        ELSE        

           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_READ_VAR: Variable objects VEC_FLT data is NOT assocaited!")
        END IF
     END IF


     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%VEC_FLT,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME))
        
        GVEC_FLT => VAR%VEC_FLT(1:RDIMS(1))
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GVEC_FLT,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))
        
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%VEC_FLT,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================


        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GVEC_FLT(1:RDIMS(1)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%VEC_FLT,1) .EQ. RDIMS(1)) THEN
              GVEC_FLT => VAR%VEC_FLT(1:RDIMS(1))
              
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GVEC_FLT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF ! YOU ARE THE DEALER
           
        
        LVEC_FLT => VAR%VEC_FLT(1:UBOUND(VAR%VEC_FLT,1))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GVEC_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GVEC_FLT,LVEC_FLT)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GVEC_FLT)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LVEC_FLT,RDIMS(1),MPI_REAL,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: VEC_FLT")

           NULLIFY(GVEC_FLT)
           ALLOCATE(GVEC_FLT(RDIMS(1)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GVEC_FLT=0.0_SPA

           IF (MYID==DEALER) GVEC_FLT=LVEC_FLT

           CALL MPI_BCAST(GVEC_FLT,RDIMS(1),MPI_REAL,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: VEC_FLTT")

           IF(MYID/=DEALER) LVEC_FLT=GVEC_FLT
           DEALLOCATE(GVEC_FLT)



        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(gvec_flt)
     nullify(lvec_flt)

!********************************************************************
! =====   ARRAY FLOATING POINT DATA
!********************************************************************
  CASE(case_ARR_FLT)

    IF(.NOT. ASSOCIATED(VAR%ARR_FLT))THEN 

        IF(ASSOCIATED(VAR%CUB_FLT))THEN
           IF(size(VAR%CUB_FLT,1)==1) VAR%ARR_FLT=>VAR%CUB_FLT(1,1:,1:)
           IF(size(VAR%CUB_FLT,2)==1) VAR%ARR_FLT=>VAR%CUB_FLT(1:,1,1:)
           IF(size(VAR%CUB_FLT,3)==1) VAR%ARR_FLT=>VAR%CUB_FLT(1:,1:,1)
        ELSE        

           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_READ_VAR: Variable objects ARR_FLT data is NOT assocaited!")
        END IF
     END IF
     

     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%ARR_FLT,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM1 does not match!")

      
        IF (UBOUND(VAR%ARR_FLT,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM2 does not match!")
        
        GARR_FLT => VAR%ARR_FLT(1:RDIMS(1),1:RDIMS(2))
        
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GARR_FLT,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))

        
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)


        MYSIZE  = UBOUND(VAR%ARR_FLT,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================

        
        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GARR_FLT(1:RDIMS(1),1:RDIMS(2)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%ARR_FLT,1) .EQ. RDIMS(1) .and.&
                & UBOUND(VAR%ARR_FLT,2) .LE. RDIMS(2) ) THEN

              GARR_FLT => VAR%ARR_FLT(1:RDIMS(1),1:RDIMS(2))
              
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA",&
                   & "Dimensions do not match allocated space")
              
           END IF
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GARR_FLT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF ! YOU ARE THE DEALER
           
        
        LARR_FLT => VAR%ARR_FLT(1:ubound(VAR%ARR_FLT,1),1:RDIMS(2))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GVEC_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GARR_FLT,LARR_FLT)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GARR_FLT)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LARR_FLT,RDIMS(1)*RDIMS(2),MPI_REAL,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: ARR_FLT")

           NULLIFY(GARR_FLT)
           ALLOCATE(GARR_FLT(RDIMS(1),RDIMS(2)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GARR_FLT=0.0_SPA
           
           IF (MYID==DEALER) GARR_FLT=LARR_FLT
           
           CALL MPI_BCAST(GARR_FLT,RDIMS(1)*RDIMS(2),MPI_REAL,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: ARR_FLT")
           
           IF(MYID/=DEALER) LARR_FLT=GARR_FLT
           DEALLOCATE(GARR_FLT)

        END IF

# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(gARR_flt)
     nullify(lARR_flt)

!********************************************************************
! =====   CUBE FLOATING POINT DATA
!********************************************************************
  CASE(case_CUB_FLT)

    IF(.NOT. ASSOCIATED(VAR%CUB_FLT))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: Variable objects CUB_FLT data is NOT assocaited!")
     END IF
     

     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%CUB_FLT,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM1 does not match!")

        IF (UBOUND(VAR%CUB_FLT,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM2 does not match!")

        IF (UBOUND(VAR%CUB_FLT,3) .LT. RDIMS(3)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM3 does not match!")
        
        GCUB_FLT => VAR%CUB_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))
        
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GCUB_FLT,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))

        
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%CUB_FLT,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================
        
        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GCUB_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")

           ELSEIF (UBOUND(VAR%CUB_FLT,1) .EQ. RDIMS(1) .and.&
                &  UBOUND(VAR%CUB_FLT,2) .LE. RDIMS(2) .and.&
                &  UBOUND(VAR%CUB_FLT,3) .LE. RDIMS(3)) THEN

              GCUB_FLT => VAR%CUB_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))              
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GCUB_FLT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF ! YOU ARE THE DEALER
           
        
        LCUB_FLT => VAR%CUB_FLT(1:ubound(VAR%CUB_FLT,1),1:RDIMS(2),1:RDIMS(3))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GVEC_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GCUB_FLT,LCUB_FLT)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GCUB_FLT)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LCUB_FLT,RDIMS(1)*RDIMS(2)*RDIMS(3),MPI_REAL,SOURCE,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: CUB_FLT")

           NULLIFY(GCUB_FLT)
           ALLOCATE(GCUB_FLT(RDIMS(1),RDIMS(2),RDIMS(3)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GCUB_FLT=0.0_SPA

           IF (MYID==DEALER) GCUB_FLT=LCUB_FLT

           CALL MPI_BCAST(GCUB_FLT,RDIMS(1)*RDIMS(2)*RDIMS(3),MPI_REAL,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: CUB_FLT")

           IF(MYID/=DEALER) LCUB_FLT=GCUB_FLT
           DEALLOCATE(GCUB_FLT)

        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(gCUB_flt)
     nullify(lCUB_flt)

!********************************************************************
! =====   FOUR DIMENSION ARRAY FLOATING POINT DATA
!********************************************************************
  CASE(case_FDA_FLT)

    IF(.NOT. ASSOCIATED(VAR%FDA_FLT))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: Variable objects FDA_FLT data is NOT assocaited!")
     END IF
     

     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%FDA_FLT,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM1 does not match!")

        IF (UBOUND(VAR%FDA_FLT,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM2 does not match!")

        IF (UBOUND(VAR%FDA_FLT,3) .LT. RDIMS(3)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM3 does not match!")
        
        IF (UBOUND(VAR%FDA_FLT,4) .LT. RDIMS(4)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM4 does not match!")
        
        GFDA_FLT => VAR%FDA_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
        
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GFDA_FLT,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))

        
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%FDA_FLT,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================
        
        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GFDA_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")

           ELSEIF (UBOUND(VAR%FDA_FLT,1) .EQ. RDIMS(1) .and.&
                &  UBOUND(VAR%FDA_FLT,2) .LE. RDIMS(2) .and.&
                &  UBOUND(VAR%FDA_FLT,3) .LE. RDIMS(3) .and.&
                &  UBOUND(VAR%FDA_FLT,4) .LE. RDIMS(4)) THEN

              GFDA_FLT => VAR%FDA_FLT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))              
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GFDA_FLT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF ! YOU ARE THE DEALER
           
        
        LFDA_FLT => VAR%FDA_FLT(1:ubound(VAR%FDA_FLT,1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GVEC_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GFDA_FLT,LFDA_FLT)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GFDA_FLT)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LFDA_FLT,RDIMS(1)*RDIMS(2)*RDIMS(3)*RDIMS(4),MPI_REAL,SOURCE,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: FDA_FLT")

           NULLIFY(GFDA_FLT)
           ALLOCATE(GFDA_FLT(RDIMS(1),RDIMS(2),RDIMS(3),RDIMS(4)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GFDA_FLT=0.0_SPA

           IF (MYID==DEALER) GFDA_FLT=LFDA_FLT

           CALL MPI_BCAST(GFDA_FLT,RDIMS(1)*RDIMS(2)*RDIMS(3)*RDIMS(4),MPI_REAL,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: FDA_FLT")

           IF(MYID/=DEALER) LFDA_FLT=GFDA_FLT
           DEALLOCATE(GFDA_FLT)

        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(gFDA_flt)
     nullify(lFDA_flt)

!********************************************************************
! =====   SCALAR DOUBLE DATA
!********************************************************************
  CASE(case_scl_DBL)
     
     IF(.NOT. ASSOCIATED(VAR%SCL_DBL))THEN 

        IF(ASSOCIATED(VAR%VEC_DBL))THEN
           IF(size(VAR%VEC_DBL)==1) VAR%SCL_DBL=>VAR%VEC_DBL(1)
        ELSE IF(ASSOCIATED(VAR%ARR_DBL))THEN
           IF(size(VAR%ARR_DBL)==1) VAR%SCL_DBL=>VAR%ARR_DBL(1,1)
        ELSE IF(ASSOCIATED(VAR%CUB_DBL))THEN
           IF(size(VAR%CUB_DBL)==1) VAR%SCL_DBL=>VAR%CUB_DBL(1,1,1)
        ELSE        
           
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_READ_VAR: Variable objects SCL_DBL data is NOT assocaited!")
        END IF

     END IF
  
     IF (SER_READ .OR. DEALER .EQ. MYID) THEN
        
        
        IF (SIZE(NSTART).GT.0) THEN

            if (product(nCOUNT) .NE. 1) CALL FATAL_ERROR&
                & ("NC_READ_VAR: NCOUNT dimension size invalid while reading scl_dbl?")

           allocate(gvec_DBL(1))
           !gvec_flt = VAR%SCL_FLT ! Was this here for a reason?
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GVEC_DBL,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
           VAR%SCL_DBL = gvec_DBL(1)
           deallocate(gvec_DBL)
           
        ELSE
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,VAR%SCL_DBL)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF
     END IF
     
     IF (PAR_READ) THEN
# if defined (MULTIPROCESSOR)
        
        SOURCE = DEALER -1 
        CALL MPI_BCAST(VAR%SCL_DBL,1,MPI_DP,SOURCE,MPI_FVCOM_GROUP,IERR)
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
        
     END IF     

     nullify(gvec_dbl)

!********************************************************************
! =====   VECTOR DOUBLE DATA
!********************************************************************
  CASE(case_vec_dbl)

    IF(.NOT. ASSOCIATED(VAR%VEC_DBL))THEN 

        IF(ASSOCIATED(VAR%ARR_DBL))THEN
           IF(size(VAR%ARR_DBL,1)==1) VAR%VEC_DBL=>VAR%ARR_DBL(1,1:)
           IF(size(VAR%ARR_DBL,2)==1) VAR%VEC_DBL=>VAR%ARR_DBL(1:,1)
        ELSE IF(ASSOCIATED(VAR%CUB_DBL))THEN
           IF(size(VAR%CUB_DBL,1)==1) THEN
              IF(size(VAR%CUB_DBL,2)==1) VAR%VEC_DBL=>VAR%CUB_DBL(1,1,1:)
              IF(size(VAR%CUB_DBL,3)==1) VAR%VEC_DBL=>VAR%CUB_DBL(1,1:,1)
           END IF
           IF(size(VAR%CUB_DBL,1)==2) THEN
              IF(size(VAR%CUB_DBL,3)==1) VAR%VEC_DBL=>VAR%CUB_DBL(1:,1,1)
           END IF           
        ELSE        
           
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_READ_VAR: Variable objects VEC_DBL data is NOT assocaited!")
        END IF
     END IF

     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%VEC_DBL,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME))
        
        GVEC_DBL => VAR%VEC_DBL(1:RDIMS(1))
          
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GVEC_DBL,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))
          
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%VEC_DBL,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================

        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GVEC_DBL(1:RDIMS(1)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%VEC_DBL,1) .EQ. RDIMS(1)) THEN
              GVEC_DBL => VAR%VEC_DBL(1:RDIMS(1))
              
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
             
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GVEC_DBL,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
             
        END IF ! YOU ARE THE DEALER
           
        
        LVEC_DBL => VAR%VEC_DBL(1:ubound(VAR%VEC_DBL,1))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GVEC_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GVEC_DBL,LVEC_DBL)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GVEC_DBL)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LVEC_DBL,RDIMS(1),MPI_DP,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: VEC_DBL")

           NULLIFY(GVEC_DBL)
           ALLOCATE(GVEC_DBL(RDIMS(1)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GVEC_DBL=0.0_DP

           IF (MYID==DEALER) GVEC_DBL=LVEC_DBL

           CALL MPI_BCAST(GVEC_DBL,RDIMS(1),MPI_DP,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT:VEC_DBL")

           IF(MYID/=DEALER) LVEC_DBL=GVEC_DBL
           DEALLOCATE(GVEC_DBL)

        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(gvec_dbl)
     nullify(lvec_dbl)


!********************************************************************
! =====   ARRAY DOUBLE DATA
!********************************************************************
  CASE(case_ARR_DBL)

    IF(.NOT. ASSOCIATED(VAR%ARR_DBL))THEN 

        IF(ASSOCIATED(VAR%CUB_DBL))THEN
           IF(size(VAR%CUB_DBL,1)==1) VAR%ARR_DBL=>VAR%CUB_DBL(1,1:,1:)
           IF(size(VAR%CUB_DBL,2)==1) VAR%ARR_DBL=>VAR%CUB_DBL(1:,1,1:)
           IF(size(VAR%CUB_DBL,3)==1) VAR%ARR_DBL=>VAR%CUB_DBL(1:,1:,1)
        ELSE        

           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_READ_VAR: Variable objects ARR_DBL data is NOT assocaited!")
        END IF
     END IF
     
     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%ARR_DBL,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM1 does not match!")

        IF (UBOUND(VAR%ARR_DBL,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM2 does not match!")

        GARR_DBL => VAR%ARR_DBL(1:RDIMS(1),1:RDIMS(2))
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GARR_DBL,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))
                
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)


        MYSIZE  = UBOUND(VAR%ARR_DBL,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================
        
        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GARR_DBL(RDIMS(1),RDIMS(2)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%ARR_DBL,1) .EQ. RDIMS(1) .and. &
                &  UBOUND(VAR%ARR_DBL,2) .LE. RDIMS(2)) THEN
         
              GARR_DBL => VAR%ARR_DBL(1:RDIMS(1),1:RDIMS(2))
              
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
                   
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GARR_DBL,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
                   
        END IF ! YOU ARE THE DEALER
           
        
        LARR_DBL => VAR%ARR_DBL(1:ubound(VAR%ARR_DBL,1),1:RDIMS(2))
        
        IF (FOUND) THEN ! THE DATA IS IN GVEC_FLT

           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GARR_DBL,LARR_DBL)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GARR_DBL)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LARR_DBL,RDIMS(1)*RDIMS(2),MPI_DP,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: ARR_DBL")

           NULLIFY(GARR_DBL)
           ALLOCATE(GARR_DBL(RDIMS(1),RDIMS(2)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GARR_DBL=0.0_DP

           IF (MYID==DEALER) GARR_DBL=LARR_DBL

           CALL MPI_BCAST(GARR_DBL,RDIMS(1)*RDIMS(2),MPI_DP,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT:ARR_DBL")

           IF(MYID/=DEALER) LARR_DBL=GARR_DBL
           DEALLOCATE(GARR_DBL)


        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

  NULLIFY(GARR_DBL)
  NULLIFY(LARR_DBL)


!********************************************************************
! =====   CUBE DOUBLE DATA
!********************************************************************
  CASE(case_CUB_DBL)

    IF(.NOT. ASSOCIATED(VAR%CUB_DBL))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: Variable objects CUB_DBL data is NOT assocaited!")
     END IF
     
     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%CUB_DBL,1) .NE.  RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM1 does not match!")
        
        IF (UBOUND(VAR%CUB_DBL,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM2 does not match!")

        IF (UBOUND(VAR%CUB_DBL,3) .LT. RDIMS(3)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM3 does not match!")
   
        GCUB_DBL => VAR%CUB_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GCUB_DBL,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))
                
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%CUB_DBL,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================

        
        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GCUB_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%CUB_DBL,1) .EQ. RDIMS(1) .and.&
                &  UBOUND(VAR%CUB_DBL,2) .LE. RDIMS(2) .and.&
                &  UBOUND(VAR%CUB_DBL,3) .LE. RDIMS(3)) THEN

              GCUB_DBL => VAR%CUB_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))    
              
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
                   
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GCUB_DBL,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
                   
        END IF ! YOU ARE THE DEALER
           
        
        LCUB_DBL => VAR%CUB_DBL(1:ubound(VAR%CUB_DBL,1),1:RDIMS(2),1:RDIMS(3))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GCUB_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GCUB_DBL,LCUB_DBL)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GCUB_DBL)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LCUB_DBL,RDIMS(1)*RDIMS(2)*RDIMS(3),MPI_DP,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: CUB_DBL")

           NULLIFY(GCUB_DBL)
           ALLOCATE(GCUB_DBL(RDIMS(1),RDIMS(2),RDIMS(3)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GCUB_DBL=0.0_DP

           IF (MYID==DEALER) GCUB_DBL=LCUB_DBL

           CALL MPI_BCAST(GCUB_DBL,RDIMS(1)*RDIMS(2)*RDIMS(3),MPI_DP,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT:CUB_DBL")

           IF(MYID/=DEALER) LCUB_DBL=GCUB_DBL
           DEALLOCATE(GCUB_DBL)


        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(gcub_dbl)
     nullify(lcub_dbl)

!********************************************************************
! =====   FOUR DIMENSION ARRAY DOUBLE DATA
!********************************************************************
  CASE(case_FDA_DBL)

    IF(.NOT. ASSOCIATED(VAR%FDA_DBL))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: Variable objects FDA_DBL data is NOT assocaited!")
     END IF
     
     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%FDA_DBL,1) .NE.  RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM1 does not match!")
        
        IF (UBOUND(VAR%FDA_DBL,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM2 does not match!")

        IF (UBOUND(VAR%FDA_DBL,3) .LT. RDIMS(3)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM3 does not match!")
   
        IF (UBOUND(VAR%FDA_DBL,4) .LT. RDIMS(4)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM4 does not match!")
   
        GFDA_DBL => VAR%FDA_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GFDA_DBL,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))
                
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%FDA_DBL,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================

        
        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GFDA_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%FDA_DBL,1) .EQ. RDIMS(1) .and.&
                &  UBOUND(VAR%FDA_DBL,2) .LE. RDIMS(2) .and.&
                &  UBOUND(VAR%FDA_DBL,3) .LE. RDIMS(3) .and.&
                &  UBOUND(VAR%FDA_DBL,4) .LE. RDIMS(4)) THEN

              GFDA_DBL => VAR%FDA_DBL(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))    
              
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
                   
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GFDA_DBL,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
                   
        END IF ! YOU ARE THE DEALER
           
        
        LFDA_DBL => VAR%FDA_DBL(1:ubound(VAR%FDA_DBL,1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GFDA_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GFDA_DBL,LFDA_DBL)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GFDA_DBL)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LFDA_DBL,RDIMS(1)*RDIMS(2)*RDIMS(3)*RDIMS(4),MPI_DP,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: FDA_DBL")

           NULLIFY(GFDA_DBL)
           ALLOCATE(GFDA_DBL(RDIMS(1),RDIMS(2),RDIMS(3),RDIMS(4)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GFDA_DBL=0.0_DP

           IF (MYID==DEALER) GFDA_DBL=LFDA_DBL

           CALL MPI_BCAST(GFDA_DBL,RDIMS(1)*RDIMS(2)*RDIMS(3)*RDIMS(4),MPI_DP,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT:FDA_DBL")

           IF(MYID/=DEALER) LFDA_DBL=GFDA_DBL
           DEALLOCATE(GFDA_DBL)


        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(gfda_dbl)
     nullify(lfda_dbl)

!********************************************************************
! =====   SCALAR INTEGER DATA
!********************************************************************
  CASE(case_scl_int)
     
     IF(.NOT. ASSOCIATED(VAR%SCL_INT))THEN 

        IF(ASSOCIATED(VAR%VEC_INT))THEN
           IF(size(VAR%VEC_INT)==1) VAR%SCL_INT=>VAR%VEC_INT(1)
        ELSE IF(ASSOCIATED(VAR%ARR_INT))THEN
           IF(size(VAR%ARR_INT)==1) VAR%SCL_INT=>VAR%ARR_INT(1,1)
        ELSE IF(ASSOCIATED(VAR%CUB_INT))THEN
           IF(size(VAR%CUB_INT)==1) VAR%SCL_INT=>VAR%CUB_INT(1,1,1)
        ELSE        
           
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_READ_VAR: Variable objects SCL_INT data is NOT assocaited!")
        END IF
     END IF
     
     IF (SER_READ .OR. DEALER .EQ. MYID) THEN
        
        
        IF (SIZE(NSTART).GT.0) THEN
           
           if (product(nCOUNT) .NE. 1) CALL FATAL_ERROR&
                & ("NC_READ_VAR: NCOUNT dimension invalid while reading scl_int?")
 
           allocate(gvec_int(1))
           gvec_int = VAR%SCL_INT
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GVEC_INT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))

           VAR%SCL_INT = gvec_int(1)
           deallocate(gvec_int)
           
        ELSE
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,VAR%SCL_INT)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF
     END IF
     
     IF (PAR_READ) THEN
# if defined (MULTIPROCESSOR)
        
        SOURCE = DEALER -1 
        CALL MPI_BCAST(VAR%SCL_INT,1,MPI_INTEGER,SOURCE&
             &,MPI_FVCOM_GROUP,IERR)
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
        
     END IF    
     nullify(gvec_int)
 
!********************************************************************
! =====   VECTOR INTEGER DATA
!********************************************************************
  CASE(case_vec_int)

    IF(.NOT. ASSOCIATED(VAR%VEC_INT))THEN 

        IF(ASSOCIATED(VAR%ARR_INT))THEN
           IF(size(VAR%ARR_INT,1)==1) VAR%VEC_INT=>VAR%ARR_INT(1,1:)
           IF(size(VAR%ARR_INT,2)==1) VAR%VEC_INT=>VAR%ARR_INT(1:,1)
        ELSE IF(ASSOCIATED(VAR%CUB_INT))THEN
           IF(size(VAR%CUB_INT,1)==1) THEN
              IF(size(VAR%CUB_INT,2)==1) VAR%VEC_INT=>VAR%CUB_INT(1,1,1:)
              IF(size(VAR%CUB_INT,3)==1) VAR%VEC_INT=>VAR%CUB_INT(1,1:,1)
           END IF
           IF(size(VAR%CUB_INT,1)==2) THEN
              IF(size(VAR%CUB_INT,3)==1) VAR%VEC_INT=>VAR%CUB_INT(1:,1,1)
           END IF           
        ELSE        
           
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_READ_VAR: Variable objects VEC_INT data is NOT assocaited!")
        END IF
     END IF

     IF (SER_READ) THEN

        IF (UBOUND(VAR%VEC_INT,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME))
        
        GVEC_INT => VAR%VEC_INT(1:RDIMS(1))
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GVEC_INT,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))
        
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%VEC_INT,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================

        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           
           IF (FOUND) THEN
              allocate(GVEC_INT(1:RDIMS(1)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%VEC_INT,1) .EQ. RDIMS(1)) THEN
              GVEC_INT => VAR%VEC_INT(1:RDIMS(1))
              
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
           END IF
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GVEC_INT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF ! YOU ARE THE DEALER
           
        
        LVEC_INT => VAR%VEC_INT(1:UBOUND(VAR%VEC_INT,1))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GVEC_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GVEC_INT,LVEC_INT)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GVEC_INT)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LVEC_INT,RDIMS(1),MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: VEC_INT")

           NULLIFY(GVEC_INT)
           ALLOCATE(GVEC_INT(RDIMS(1)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GVEC_INT=0

           IF (MYID==DEALER) GVEC_INT=LVEC_INT

           CALL MPI_BCAST(GVEC_INT,RDIMS(1),MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT:VEC_INT")

           IF(MYID/=DEALER) LVEC_INT=GVEC_INT
           DEALLOCATE(GVEC_INT)

        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(lvec_int)
     nullify(gvec_int)
!********************************************************************
! =====   ARRAY INTEGER DATA
!********************************************************************
  CASE(case_ARR_INT)

    IF(.NOT. ASSOCIATED(VAR%ARR_INT))THEN 

        IF(ASSOCIATED(VAR%CUB_INT))THEN
           IF(size(VAR%CUB_INT,1)==1) VAR%ARR_INT=>VAR%CUB_INT(1,1:,1:)
           IF(size(VAR%CUB_INT,2)==1) VAR%ARR_INT=>VAR%CUB_INT(1:,1,1:)
           IF(size(VAR%CUB_INT,3)==1) VAR%ARR_INT=>VAR%CUB_INT(1:,1:,1)
        ELSE        
           
           CALL PRINT_VAR(VAR)
           CALL FATAL_ERROR("NC_READ_VAR: Variable objects ARR_INT data is NOT assocaited!")
        END IF
     END IF

     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%ARR_INT,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM1 does not match!")
        
        IF (UBOUND(VAR%ARR_INT,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM2 does not match!")

        GARR_INT => VAR%ARR_INT(1:RDIMS(1),1:RDIMS(2))
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GARR_INT,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))
        
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%ARR_INT,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE

        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================

        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GARR_INT(RDIMS(1),RDIMS(2)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%ARR_INT,1) .EQ. RDIMS(1) .and. &
                &  UBOUND(VAR%ARR_INT,2) .LE. RDIMS(2)) THEN

              GARR_INT => VAR%ARR_INT(1:RDIMS(1),1:RDIMS(2))
                        
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GARR_INT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF ! YOU ARE THE DEALER
           
        
        LARR_INT => VAR%ARR_INT(1:ubound(VAR%ARR_INT,1),1:RDIMS(2))
     
        
        IF (FOUND) THEN ! THE DATA IS IN GVEC_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GARR_INT,LARR_INT)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GARR_INT)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%ARR_INT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LARR_INT,RDIMS(1)*RDIMS(2),MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: ARR_INT")

           NULLIFY(GARR_INT)
           ALLOCATE(GARR_INT(RDIMS(1),RDIMS(2)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GARR_INT=0

           IF (MYID==DEALER) GARR_INT=LARR_INT

           CALL MPI_BCAST(GARR_INT,RDIMS(1)*RDIMS(2),MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: ARR_INT")

           IF(MYID/=DEALER) LARR_INT=GARR_INT
           DEALLOCATE(GARR_INT)

        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(lARR_int)
     nullify(gARR_int)

!********************************************************************
! =====   CUBE INTEGER DATA
!********************************************************************
  CASE(case_CUB_INT)

    IF(.NOT. ASSOCIATED(VAR%CUB_INT))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: Variable objects CUB_INT data is NOT assocaited!")
     END IF
     
     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%CUB_INT,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM1 does not match!")
        
        IF (UBOUND(VAR%CUB_INT,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM2 does not match!")

        IF (UBOUND(VAR%CUB_INT,3) .LT. RDIMS(3)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM3 does not match!")

        GCUB_INT => VAR%CUB_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GCUB_INT,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))
        
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%CUB_INT,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================


        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GCUB_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%CUB_INT,1) .EQ. RDIMS(1) .and.&
                &  UBOUND(VAR%CUB_INT,2) .LE. RDIMS(2) .and.&
                &  UBOUND(VAR%CUB_INT,3) .LE. RDIMS(3)) THEN

              GCUB_INT => VAR%CUB_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3))
          
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GCUB_INT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF ! YOU ARE THE DEALER
           
        
        LCUB_INT => VAR%CUB_INT(1:ubound(VAR%CUB_INT,1),1:RDIMS(2),1:RDIMS(3))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GCUB_FLT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GCUB_INT,LCUB_INT)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GCUB_INT)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LCUB_INT,RDIMS(1)*RDIMS(2)*RDIMS(3),MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: CUB_INT")

           NULLIFY(GCUB_INT)
           ALLOCATE(GCUB_INT(RDIMS(1),RDIMS(2),RDIMS(3)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GCUB_INT=0

           IF (MYID==DEALER) GCUB_INT=LCUB_INT

           CALL MPI_BCAST(GCUB_INT,RDIMS(1)*RDIMS(2)*RDIMS(3),MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT:CUB_INT")

           IF(MYID/=DEALER) LCUB_INT=GCUB_INT
           DEALLOCATE(GCUB_INT)


        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(lCUB_int)
     nullify(gCUB_int)

!********************************************************************
! =====   FOUR DIMENSION ARRAY INTEGER DATA
!********************************************************************
  CASE(case_FDA_INT)

    IF(.NOT. ASSOCIATED(VAR%FDA_INT))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: Variable objects FDA_INT data is NOT assocaited!")
     END IF
     
     IF (SER_READ) THEN
        
        IF (UBOUND(VAR%FDA_INT,1) .NE. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM1 does not match!")
        
        IF (UBOUND(VAR%FDA_INT,2) .LT. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM2 does not match!")

        IF (UBOUND(VAR%FDA_INT,3) .LT. RDIMS(3)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM3 does not match!")

        IF (UBOUND(VAR%FDA_INT,4) .LT. RDIMS(4)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE DATA IN SERIAL CASE",& 
             & "varname: "//TRIM(var%VARNAME),&
             & "DIM4 does not match!")

        GFDA_INT => VAR%FDA_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
        
        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GFDA_INT,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))
        
     ELSE ! PAR CASE 
        
# if defined (MULTIPROCESSOR)

        MYSIZE  = UBOUND(VAR%FDA_INT,1)

        ! LOOK TO SEE IF THERE IS A MAP FOR THIS DATA
        !======================================================================
        ! ALLOCATE SPACE FOR LOCAL ARRAY SIZE - KEEP NPROCS_TOTAL TO
        ! MATCH MAP ALLOCATION SIZE
        ALLOCATE(Lsizes(NPROCS_TOTAL)); LSIZES=0        
        ! FOR NC_READ_VAR DO NOT INLCUDE THE IOPROC - MPI_FVCOM_GROUP
        CALL MPI_ALLGATHER(MYSIZE,1,MPI_INTEGER,LSizes,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
        ! DO NOT USE ALLOCATED SIZE OF IOPROC DATA - THAT IS THE GLOBAL STORAGE
        
        NPsize => LSizes(1:nprocs)
        ! DEAL OPERATION USE THE HALO MAP AS DEFAULT
        GMAP => FIND_MAP(HALO_MAPS,RDIMS(1),NPsize,FOUND)
        IF(.NOT.FOUND) GMAP => FIND_MAP(INTERNAL_MAPS,RDIMS(1),NPsize,FOUND)

        nullify(NPsize)
        DEALLOCATE(LSIZES)
        !======================================================================


        ! IF YOU ARE THE DEALER
        IF(DEALER .EQ. MYID) THEN
           
           IF (FOUND) THEN
              allocate(GFDA_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4)),stat=status)
              if (0/=status) CALL FATAL_ERROR &
                   & ("NC_READ_VAR: COULD NOT ALLOCATE SPACE FOR READ")
           ELSEIF (UBOUND(VAR%FDA_INT,1) .EQ. RDIMS(1) .and.&
                &  UBOUND(VAR%FDA_INT,2) .LE. RDIMS(2) .and.&
                &  UBOUND(VAR%FDA_INT,3) .LE. RDIMS(3) .and.&
                &  UBOUND(VAR%FDA_INT,4) .LE. RDIMS(4)) THEN

              GFDA_INT => VAR%FDA_INT(1:RDIMS(1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
          
           ELSE
              CALL FATAL_ERROR &
                   & ("NC_READ_VAR: CAN NOT FIND MAP TO READ DATA")
              
           END IF
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,GFDA_INT,NSTART,NCOUNT,NSTRIDE)
           CALL HANDLE_NCERR(status,trim(errmsg))
           
        END IF ! YOU ARE THE DEALER
           
        
        LFDA_INT => VAR%FDA_INT(1:ubound(VAR%FDA_INT,1),1:RDIMS(2),1:RDIMS(3),1:RDIMS(4))
        
        
        IF (FOUND) THEN ! THE DATA IS IN GFDA_INT
           CALL PDEAL(MYID,DEALER,NPROCS,GMAP,GFDA_INT,LFDA_INT)
           
           IF (MYID .EQ. DEALER) DEALLOCATE(GFDA_INT)

        ELSE ! THE DATA IS ALREADY LOADED IN VAR%VEC_FLT
           SOURCE = DEALER -1 

           ! CAN NOT PASS A POINTER WHICH MIGHT NOT USE CONTIGUOUS MEMORY - THIS
           ! IS A BUG IN IFORT/MVAPICH
           ! CALL MPI_BCAST(LFDA_INT,RDIMS(1)*RDIMS(2)*RDIMS(3)*RDIMS(4),MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
           ! IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT: FDA_INT")

           NULLIFY(GFDA_INT)
           ALLOCATE(GFDA_INT(RDIMS(1),RDIMS(2),RDIMS(3),RDIMS(4)),stat=STATUS)
           IF(STATUS/=0) &
                &CALL FATAL_ERROR("CAN'T ALLOCATE TEMPORARY ARRAY FOR BROADCAST OF DATA FROM NETCDF FILE")
           GFDA_INT=0

           IF (MYID==DEALER) GFDA_INT=LFDA_INT

           CALL MPI_BCAST(GFDA_INT,RDIMS(1)*RDIMS(2)*RDIMS(3)*RDIMS(4),MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
           IF(IERR/=0) CALL FATAL_ERROR("NC_READ_VAR: COULD NOT BROADCAST RESULT:FDA_INT")

           IF(MYID/=DEALER) LFDA_INT=GFDA_INT
           DEALLOCATE(GFDA_INT)


        END IF
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF

     nullify(lFDA_int)
     nullify(gFDA_int)

!********************************************************************
! =====   SCALAR STRING DATA
!********************************************************************
  CASE(case_scl_chr)

    IF(.NOT. ASSOCIATED(VAR%scl_chr))THEN 

       IF (ASSOCIATED(VAR%vec_chr))THEN 
          IF(SIZE(VAR%vec_chr)==1) VAR%scl_chr => VAR%vec_chr(1)
       ELSE
          CALL PRINT_VAR(VAR)
          CALL FATAL_ERROR("NC_READ_VAR: Variable objects scl_chr&
               & data is NOT assocaited!")
       END IF

     END IF
     

     IF (SER_READ .OR. DEALER .EQ. MYID) THEN
        
        IF (LEN(VAR%SCL_CHR) .LT. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE CHARACTER DATA",& 
             & "varname: "//TRIM(var%VARNAME))

        var%scl_chr = ACHAR(0) ! INITIALIZE with null char
!        scl_chr => var%scl_chr(1:RDIMS(1))
        scl_chr => var%scl_chr

! NOTE: TO PRINT POINTER CHARACTER YOU MUST USE TRIM(chr) OR chr(:) !         
!        scl_chr = "testing"
!        write(ipt,*) "a test: '"//scl_chr(:)//"'"

        status = NF90_GET_VAR(VAR%NCID,VAR%VARID,scl_chr,NSTART,NCOUNT,NSTRIDE)
        CALL HANDLE_NCERR(status,trim(errmsg))

        nlen = index(scl_chr,ACHAR(0))
        ! If no null byte is found use the whole string
        if (nlen == 0) then
           nlen = RdimS(1) 
        else
           nlen = nlen -1
           scl_chr = scl_chr(1:nlen)
        end if

     END IF
        
     IF(PAR_READ) THEN
# if defined (MULTIPROCESSOR)
        SOURCE = DEALER -1 
        scl_chr => var%scl_chr
        CALL MPI_BCAST(nlen,1,MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
        CALL MPI_BCAST(scl_chr,nlen,MPI_CHARACTER,SOURCE,MPI_FVCOM_GROUP,IERR)
        scl_chr = scl_chr(1:nlen)
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF
     nullify(scl_chr)

!********************************************************************
! =====   VECTOR STRING DATA
!********************************************************************
  CASE(case_vec_chr)

    IF(.NOT. ASSOCIATED(VAR%vec_chr))THEN 
        CALL PRINT_VAR(VAR)
        CALL FATAL_ERROR("NC_READ_VAR: Variable objects vec_chr&
             & data is NOT assocaited!")
     END IF

     IF (SER_READ .OR. DEALER .EQ. MYID) THEN
        
        IF (SIZE(Var%VEC_CHR) .NE. RDIMS(2)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED Character array DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE CHARACTER DATA",& 
             & "varname: "//TRIM(var%VARNAME))

        IF (LEN(VAR%VEC_CHR(1)) .LT. RDIMS(1)) CALL FATAL_ERROR &
             & ("NC_READ_VAR: FILE OBJECT ALLOCATED SPACE DOES NOT MATCH&
             & THE DIMENSIONS OF THE FILE CHARACTER DATA",& 
             & "varname: "//TRIM(var%VARNAME))
        
        CNT=SIZE(NCOUNT)
        allocate(nSTRT(CNT),nCNT(CNT),nSTRD(CNT))
        nSTRT=NSTART
        nCNT=NCOUNT
        nSTRD=NSTRIDE

        Do I = 1,RDIMS(2)
           
           
           VAR%vec_chr(I) = ACHAR(0) ! INITIALIZE with null char
           scl_chr => VAR%vec_chr(I)
         
           ! READ THE I'th Entry (one entry)
           nSTRT(2)=I
           nCNT(2)=1
           
           status = NF90_GET_VAR(VAR%NCID,VAR%VARID,SCL_CHR,NSTRT,NCNT,NSTRD)
           CALL HANDLE_NCERR(status,trim(errmsg))

           nlen = index(scl_chr,ACHAR(0))
           ! If no null byte is found use the whole string
           if (nlen == 0) then
              nlen = RdimS(1) 
           else
              nlen = nlen -1
              scl_chr = scl_chr(1:nlen)
           end if
           
           nullify(scl_chr)
           
        END DO
        deallocate(nSTRT,nCNT,nSTRD)

     END IF
        
     IF(PAR) THEN
# if defined (MULTIPROCESSOR)
        DO I = 1,RDIMS(2)
           SOURCE = DEALER -1 
           scl_chr => var%vec_chr(I)
           nlen = len_Trim(scl_chr)
           CALL MPI_BCAST(nlen,1,MPI_INTEGER,SOURCE,MPI_FVCOM_GROUP,IERR)
           CALL MPI_BCAST(scl_chr,nlen,MPI_CHARACTER,SOURCE,MPI_FVCOM_GROUP,IERR)
           scl_chr = scl_chr(1:nlen)
        END DO
# else
        CALL FATAL_ERROR("NC_READ_VAR: TRIED TO ENTER PARALLEL OPTION IN SERIAL CODING!")
# endif
     END IF
     
     nullify(scl_chr)
     
  CASE default
     call print_var(VAR)
     call Fatal_error("NC_WRITE_VAR: UNKNOWN CASE")
     
  END SELECT


  !Only deallocate if it is not pointing to Ncount
  IF(.not.ASSOCIATED(RDIMS,NCOUNT)) THEN
     IF(ASSOCIATED(RDIMS)) DEALLOCATE(RDIMS)
  ELSE
     NULLIFY(RDIMS)
  END IF


  IF(PRESENT(IOSTART)) THEN
     nullify(nstart)
  ELSE
     deallocate(nstart)
  END IF

  IF(PRESENT(IOCOUNT)) THEN
     nullify(nCOUNT)
  ELSE
     deallocate(nCOUNT)
  END IF

  IF(PRESENT(IOSTRIDE)) THEN
     nullify(nstride)
  ELSE
     deallocate(nstride)
  END IF


  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END NC_READ_VAR:"


END SUBROUTINE NC_READ_VAR
!====================================================================
!====================================================================
FUNCTION IS_VALID_DATETIME(VAR,tzone) RESULT(RES)
  IMPLICIT NONE
  TYPE(NCATT), POINTER :: ATT
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM
  logical :: found, res
  character(len=80), intent(out) :: tzone
  
  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("IS_VALID_DATETIME: Variable object argument is not assocai&
       &ted!")


  RES = .FALSE.

  IF (VAR%XTYPE /= NF90_CHAR) RETURN

  ATT => FIND_ATT(VAR,'description',FOUND)
  IF (FOUND) THEN
     IF (ATT%CHR(1) == "GMT time")THEN
        RES =.TRUE.
        tzone="UTC"
        RETURN
     END IF
  END IF
  
  ATT => FIND_ATT(VAR,'time_zone',FOUND)
  IF (FOUND) THEN
     IF (.not. is_valid_timezone(ATT%chr(1))) return
     tzone=ATT%chr(1)
     RES=.TRUE.
  END IF
  

  RETURN
 
END FUNCTION IS_VALID_DATETIME
!====================================================================
!====================================================================
 FUNCTION DATETIME_OBJECT(DIMSTR,DIMTIME,timezone,size) RESULT(VAR)
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER:: DIMSTR
    TYPE(NCDIM),  POINTER, OPTIONAL :: DIMTIME
    INTEGER, OPTIONAL :: SIZE
    TYPE(NCATT), POINTER :: ATT
    CHARACTER(LEN=80),pointer :: Data_vec(:)
    CHARACTER(LEN=80),pointer :: Data_scl
    CHARACTER(LEN=*), optional :: timezone

    IF(PRESENT(SIZE)) THEN
       ALLOCATE(DATA_vec(SIZE))
    ELSE
       ALLOCATE(DATA_vec(1))
       DATA_scl =>DATA_vec(1)
    END IF

    IF (PRESENT(DIMTIME)) THEN
       VAR  => NC_MAKE_PVAR(name='Times', values=Data_vec, DIM1= DIMSTR, DIM2= DIMtime)
       VAR%SCL_CHR => VAR%VEC_CHR(1)
    ELSE
       VAR  => NC_MAKE_PVAR(name='Times', values=Data_scl, DIM1= DIMStr)
    END IF

    IF (PRESENT(TIMEZONE))THEN
       ATT  => NC_MAKE_ATT(name='time_zone',values=trim(timezone)) 
       VAR  => ADD(VAR,ATT)
    ELSE
       ATT  => NC_MAKE_ATT(name='time_zone',values='UTC') 
       VAR  => ADD(VAR,ATT)
    END IF

    
  END FUNCTION DATETIME_OBJECT
!====================================================================
!====================================================================
  SUBROUTINE UPDATE_DATETIME(VAR,NOW)
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT), POINTER :: ATT
    TYPE(TIME), INTENT(in) :: NOW
    CHARACTER(len=80), POINTER :: Data
    LOGICAL :: TEST
    CHARACTER(len=80):: TZONE
    CHARACTER(len=80):: TEMP
    

    TEST = IS_VALID_DateTime(VAR,tzone)
    IF(.not. TEST) THEN
       CALL PRINT_VAR(VAR)
       CALL FATAL_ERROR &
            ("CAN NOT UPDATE TIME FOR INVALID DATE TIME VARIABLE")
    END IF

    CALL NC_POINT_VAR(VAR,Data)

     TEMP = WRITE_DATETIME(NOW,TimePrec,tzone)
     Data = TRIM(TEMP)

  END SUBROUTINE UPDATE_DATETIME
!====================================================================
!====================================================================
FUNCTION IS_VALID_ITIME(VAR1,VAR2,tzone) RESULT(RES)
  IMPLICIT NONE
  character(len=80), intent(out) :: tzone
  TYPE(NCATT), POINTER :: ATT
  TYPE(NCVAR), POINTER :: VAR1,VAR2
  TYPE(NCDIM), POINTER :: DIM
  logical :: found, res

  IF(.NOT. ASSOCIATED(VAR1)) CALL FATAL_ERROR &
       & ("IS_VALID_INT2_MJD: Variable object argument is not assocaited!")

  IF(.NOT. ASSOCIATED(VAR2)) CALL FATAL_ERROR &
       & ("IS_VALID_INT2_MJD: Variable object argument is not assocaited!")

  RES = .false.

  IF (VAR1%XTYPE /= NF90_INT) RETURN
  IF (VAR2%XTYPE /= NF90_INT) RETURN

!!$  ! CHECK ATTS FOR THE DAYS VARIABLE  
!!$  DIM => FIND_DIM(VAR1,'time',FOUND)
!!$  IF (.not. FOUND) CALL WARNING &
!!$       & ("IS_VALID_INT2_MJD: VARIABLE PASSED DOES NOT HAVE DIMENSION&
!!$       & 'time'", "VARNAME: "//TRIM(VAR1%VARNAME))
!!$
!!$  IF (.NOT. DIM%UNLIMITED) CALL WARNING &
!!$       & ("IS_VALID_INT2_MJD: VARIABLE PASSED HAS DIMENSION&
!!$       & 'time' BUT IT IS NOT AN 'UNLIMITED' DIMENSION", "VARNAME: "//TRIM(VAR1%VARNAME))


  ATT => find_att(VAR1,'units',FOUND)
  IF(.not. FOUND) return

  IF(mjd_units .eq. ATT%chr(1)(1:len_trim(mjd_units)) .or.     &
     'days since '//trim(DATE_REFERENCE) .eq. ATT%chr(1)(1:len_trim(mjd_units))) THEN
!JQI  IF(mjd_units .eq. ATT%chr(1)(1:len_trim(mjd_units))) THEN
     
     ATT => find_att(VAR1,'format',FOUND)
     IF(.not. FOUND) return
     
     IF (ATT%chr(1)(1:len_trim(fmat)) .NE. fmat .and.       &
         ATT%chr(1)(1:len_trim(rfmat)) .NE. rfmat) return
!JQI     IF (ATT%chr(1)(1:len_trim(fmat)) .NE. fmat ) return
     
     ATT => find_att(VAR1,'time_zone',FOUND)
     IF(.not. FOUND) return
     tzone=ATT%chr(1)
     
  ELSE IF(days_units .eq. ATT%chr(1)(1:len_trim(days_units))) THEN

     ATT => find_att(VAR1,'time_zone',FOUND)
     IF(.not. FOUND) return
     tzone = TRIM(ATT%chr(1))
     IF (tzone /= 'none') return
  ELSE
     RETURN

  END IF
     


  
!!$  ! CHECK ATTS FOR THE Mili Seconds VARIABLE  
!!$  DIM => FIND_DIM(VAR2,'time',FOUND)
!!$  IF (.not. FOUND) CALL WARNING &
!!$       & ("IS_VALID_INT2_MJD: VARIABLE PASSED DOES NOT HAVE DIMENSION&
!!$       & 'time'", "VARNAME: "//TRIM(VAR2%VARNAME))
!!$
!!$  IF (.NOT. DIM%UNLIMITED) CALL WARNING &
!!$       & ("IS_VALID_INT2_MJD: VARIABLE PASSED HAS DIMENSION&
!!$       & 'time' BUT IT IS NOT AN 'UNLIMITED' DIMENSION", "VARNAME: "//TRIM(VAR2%VARNAME))

  ATT => find_att(VAR2,'units',FOUND)
  IF(.not. FOUND) return
  
  IF(ATT%chr(1)(1:len_trim(msec_units)) .NE. msec_units) RETURN


  ATT => find_att(VAR2,'time_zone',FOUND)
  IF(.not. FOUND) return


  ! TIME ZONE MUST BE THE SAME FOR BOTH ITEME'S
  IF (trim(ATT%chr(1)) /= trim(tzone)) RETURN


  res = .true.
  return
  

END FUNCTION IS_VALID_ITIME
!====================================================================
!====================================================================
 FUNCTION ITIME_OBJECT(use_mjd,DIM,size) RESULT(VAR)
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    logical, intent(in) :: use_mjd
    TYPE(NCDIM),  POINTER, OPTIONAL :: DIM
    INTEGER, OPTIONAL :: SIZE
    TYPE(NCATT), POINTER :: ATT
    INTEGER,POINTER :: Data_vec(:)
    INTEGER,POINTER :: Data_scl

    IF(PRESENT(SIZE)) THEN
       ALLOCATE(DATA_vec(SIZE))
    ELSE
       ALLOCATE(DATA_vec(1))
       DATA_scl =>DATA_vec(1)
    END IF
    
        ! Itime
    IF (PRESENT(DIM)) THEN
        VAR  => NC_MAKE_PVAR(name='Itime', values=Data_vec, DIM1= DIM)
        VAR%SCL_INT => VAR%VEC_INT(1)
    ELSE
       VAR  => NC_MAKE_PVAR(name='Itime', values=Data_scl)
    END IF

    IF (use_mjd) THEN
      IF(DATE_REFERENCE == 'default')THEN

       ATT  => NC_MAKE_ATT(name='units',values=mjd_units)
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='format',values=fmat)
       VAR  => ADD(VAR,ATT)

      ELSE

       ATT  => NC_MAKE_ATT(name='units',values='days since '//trim(DATE_REFERENCE))
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='format',values=rfmat)
       VAR  => ADD(VAR,ATT)

      END IF
       
       ATT  => NC_MAKE_ATT(name='time_zone',values='UTC') 
       VAR  => ADD(VAR,ATT)
    ELSE
       ATT  => NC_MAKE_ATT(name='units',values=days_units)
       VAR  => ADD(VAR,ATT)
              
       ATT  => NC_MAKE_ATT(name='time_zone',values='none') 
       VAR  => ADD(VAR,ATT)
    END IF

  END FUNCTION ITIME_OBJECT
!====================================================================
!====================================================================
 FUNCTION ITIME2_OBJECT(use_mjd,DIM,size) RESULT(VAR)
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    logical, intent(in) :: use_mjd
    TYPE(NCDIM),  POINTER, OPTIONAL :: DIM
    INTEGER, OPTIONAL :: SIZE
    TYPE(NCATT), POINTER :: ATT
    INTEGER,POINTER :: Data_vec(:)
    INTEGER,POINTER :: Data_scl

    IF(PRESENT(SIZE)) THEN
       ALLOCATE(DATA_vec(SIZE))
    ELSE
       ALLOCATE(DATA_vec(1))
       DATA_scl =>DATA_vec(1)
    END IF

    IF (PRESENT(DIM)) THEN
       
       VAR  => NC_MAKE_PVAR(name='Itime2', values=Data_vec, DIM1= DIM)
       VAR%SCL_INT => VAR%VEC_INT(1)
    ELSE
       VAR  => NC_MAKE_PVAR(name='Itime2', values=Data_scl)
    END IF

    IF (use_mjd) THEN
       ATT  => NC_MAKE_ATT(name='units',values=msec_units) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='time_zone',values='UTC') 
       VAR  => ADD(VAR,ATT)
    ELSE
       ATT  => NC_MAKE_ATT(name='units',values=msec_units) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='time_zone',values='none') 
       VAR  => ADD(VAR,ATT)
    END IF


  END FUNCTION ITIME2_OBJECT
!====================================================================
!====================================================================
  SUBROUTINE UPDATE_ITIME(VAR1,VAR2,NOW)
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR1
    TYPE(NCVAR),  POINTER :: VAR2
    TYPE(NCATT), POINTER :: ATT
    TYPE(TIME), INTENT(in) :: NOW
    INTEGER, POINTER :: D1,D2
    LOGICAL :: TEST2
    INTEGER :: TEST
    CHARACTER(len=80):: TZONE

    TEST2 = IS_VALID_ITIME(VAR1,VAR2,tzone)
    IF(.not. TEST2) THEN
       CALL PRINT_VAR(VAR1)
       CALL PRINT_VAR(VAR2)
       CALL FATAL_ERROR &
            ("CAN NOT UPDATE TIME FOR INVALID INTEGER TIME VARIABLES")
    END IF

    CALL NC_POINT_VAR(VAR1,D1)

    CALL NC_POINT_VAR(VAR2,D2)

    TEST = TIME2NCITIME(NOW,ReferenceDate,D1,D2)

!    if(.not. TEST) call fatal_error("That is bad times man!") 
    if(TEST==0) call fatal_error("That is bad times man!")
    ! THIS SHOULD NEVER HAPPEN?
    
  END SUBROUTINE UPDATE_ITIME
!====================================================================
!====================================================================
FUNCTION IS_VALID_FLOAT_DAYS(VAR,tzone) RESULT(RES)
  IMPLICIT NONE
  character(len=80), intent(out) :: tzone
  TYPE(NCATT), POINTER :: ATT
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM
  logical :: found, res

  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("IS_VALID_FLOAT_MJD: Variable object argument is not assocaited!")
  
  RES = .false.

!!$  DIM => FIND_DIM(VAR,'time',FOUND)
!!$  IF (.not. FOUND) CALL WARNING &
!!$       & ("IS_VALID_FLOAT_MJD: VARIABLE PASSED DOES NOT HAVE DIMENSION&
!!$       & 'time'", "VARNAME: "//TRIM(VAR%VARNAME))
!!$
!!$  IF (.NOT. DIM%UNLIMITED) CALL WARNING &
!!$       & ("IS_VALID_FLOAT_MJD: VARIABLE PASSED HAS DIMENSION&
!!$       & 'time' BUT IT IS NOT AN 'UNLIMITED' DIMENSION", "VARNAME: "//TRIM(VAR%VARNAME))


  ! CHECK ATTS FOR THE TIME VARIABLE  
  ATT => find_att(VAR,'units',FOUND)
  IF(.not. FOUND) return
  
  IF(mjd_units .eq. ATT%chr(1)(1:len_trim(mjd_units)) .or.     &
     'days since '//trim(DATE_REFERENCE) .eq. ATT%chr(1)(1:len_trim(mjd_units))) THEN
!JQI  IF(mjd_units .eq. ATT%chr(1)(1:len_trim(mjd_units)) ) THEN
       
     ATT => find_att(VAR,'format',FOUND)
     IF(.not. FOUND) return
     
     IF (ATT%chr(1)(1:len_trim(fmat)) .NE. fmat .and.      &
         ATT%chr(1)(1:len_trim(rfmat)) .NE. rfmat) return
!JQI     IF (ATT%chr(1)(1:len_trim(fmat)) .NE. fmat ) return
     
     ATT => find_att(VAR,'time_zone',FOUND)
     IF(.not. FOUND) return
     
     IF (.not. is_valid_timezone(ATT%chr(1))) return
     tzone=TRIM(ATT%chr(1))
     
  ELSE IF(days_units .eq. ATT%chr(1)(1:len_trim(days_units))) THEN
     
     ATT => find_att(VAR,'time_zone',FOUND)
     IF(.not. FOUND) return
     tzone = TRIM(ATT%chr(1))
     IF (tzone /= 'none') return
     
  ELSE
     return
  END IF

  res = .true.
  return
END FUNCTION IS_VALID_FLOAT_DAYS
!====================================================================
!====================================================================
FUNCTION IS_VALID_FLOAT_SECONDS(VAR,tzone) RESULT(RES)
  IMPLICIT NONE
  character(len=80), intent(out) :: tzone
  TYPE(NCATT), POINTER :: ATT
  TYPE(NCVAR), POINTER :: VAR
  TYPE(NCDIM), POINTER :: DIM
  logical :: found, res

  IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
       & ("IS_VALID_FLOAT_MJD: Variable object argument is not assocaited!")
  
  RES = .false.

!!$  DIM => FIND_DIM(VAR,'time',FOUND)
!!$  IF (.not. FOUND) CALL WARNING &
!!$       & ("IS_VALID_FLOAT_MJD: VARIABLE PASSED DOES NOT HAVE DIMENSION&
!!$       & 'time'", "VARNAME: "//TRIM(VAR%VARNAME))
!!$
!!$  IF (.NOT. DIM%UNLIMITED) CALL WARNING &
!!$       & ("IS_VALID_FLOAT_MJD: VARIABLE PASSED HAS DIMENSION&
!!$       & 'time' BUT IT IS NOT AN 'UNLIMITED' DIMENSION", "VARNAME: "//TRIM(VAR%VARNAME))


  ! CHECK ATTS FOR THE TIME VARIABLE  
  ATT => find_att(VAR,'units',FOUND)
  IF(.not. FOUND) return
  
  IF(seconds_units .eq. ATT%chr(1)(1:len_trim(seconds_units))) THEN
       
     ATT => find_att(VAR,'time_zone',FOUND)
     IF(.not. FOUND) THEN
        tzone = 'none'
     Else
        tzone = TRIM(ATT%chr(1))
     END IF
  ELSE
     return
  END IF

  res = .true.
  return
END FUNCTION IS_VALID_FLOAT_SECONDS
!====================================================================
!====================================================================
 FUNCTION Float_time_OBJECT(use_mjd,DIM,size) RESULT(VAR)
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    logical, intent(in) :: use_mjd
    TYPE(NCDIM),  POINTER, OPTIONAL :: DIM
    INTEGER, OPTIONAL :: SIZE
    TYPE(NCATT), POINTER :: ATT
    REAL(SP),pointer :: Data_vec(:)
    REAL(SP),pointer :: Data_scl

    IF(PRESENT(SIZE)) THEN
       ALLOCATE(DATA_vec(SIZE))
    ELSE
       ALLOCATE(DATA_vec(1))
       DATA_scl =>DATA_vec(1)
    END IF

    IF (PRESENT(DIM)) THEN
       VAR  => NC_MAKE_PVAR(name='time', values=Data_vec, DIM1= DIM)
       if(associated(var%vec_flt))then
         VAR%scl_flt => Var%vec_flt(1)
       else
         VAR%scl_dbl => Var%vec_dbl(1)
       endif
    ELSE
       VAR  => NC_MAKE_PVAR(name='time', values=Data_scl)
    END IF
    
    ATT  => NC_MAKE_ATT(name='long_name',values='time')
    VAR  => ADD(VAR,ATT)
    
    IF (use_mjd) THEN
      IF(DATE_REFERENCE == 'default')THEN

       ATT  => NC_MAKE_ATT(name='units',values=mjd_units)
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='format',values=fmat)
       VAR  => ADD(VAR,ATT)

      ELSE

       ATT  => NC_MAKE_ATT(name='units',values='days since '//trim(DATE_REFERENCE))
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='format',values=rfmat)
       VAR  => ADD(VAR,ATT)

      END IF
       
       ATT  => NC_MAKE_ATT(name='time_zone',values='UTC') 
       VAR  => ADD(VAR,ATT)
       
    ELSE
       ATT  => NC_MAKE_ATT(name='units',values=days_units)
       VAR  => ADD(VAR,ATT)       
       
       ATT  => NC_MAKE_ATT(name='time_zone',values='none') 
       VAR  => ADD(VAR,ATT)
    END IF
    

  END FUNCTION FLOAT_TIME_OBJECT
!====================================================================
!====================================================================
  SUBROUTINE UPDATE_FLOAT_TIME(VAR,NOW)
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(TIME), INTENT(in) :: NOW
    REAL(SP), POINTER :: Data
    LOGICAL :: TEST
    CHARACTER(len=80):: TZONE

    TEST = IS_VALID_FLOAT_DAYS(VAR,TZONE)
    IF(.not. TEST) THEN
       CALL PRINT_VAR(VAR)
       call print_att_list(VAR)
       CALL FATAL_ERROR &
            ("CAN NOT UPDATE TIME FOR INVALID FLOATING POINT TIME VARIABLE")
    END IF

    CALL NC_POINT_VAR(VAR,Data)
    
    Data = DAYS(NOW) - DAYS(ReferenceDate)

  END SUBROUTINE UPDATE_FLOAT_TIME
!====================================================================
!====================================================================
 FUNCTION IINT_OBJECT(DIM,size) RESULT(VAR)
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER, OPTIONAL :: DIM
    INTEGER, OPTIONAL :: SIZE
    TYPE(NCATT), POINTER :: ATT
    INTEGER,POINTER :: Data_vec(:)
    INTEGER,POINTER :: Data_scl

    IF(PRESENT(SIZE)) THEN
       ALLOCATE(DATA_vec(SIZE))
    ELSE
       ALLOCATE(DATA_vec(1))
       DATA_scl =>DATA_vec(1)
    END IF
    
    ! IINT
    IF (PRESENT(DIM)) THEN
        VAR  => NC_MAKE_PVAR(name='iint', values=Data_vec, DIM1= DIM)
        VAR%SCL_INT => VAR%VEC_INT(1)
    ELSE
       VAR  => NC_MAKE_PVAR(name='iint', values=Data_scl)
    END IF
    
    ATT  => NC_MAKE_ATT(name='long_name',values='internal mode iteration number') 
    VAR  => ADD(VAR,ATT)

  END FUNCTION IINT_OBJECT
!====================================================================
!====================================================================
 RECURSIVE SUBROUTINE RECURSIVE_FILE_BRACKET(FTM,NOW,status)
  IMPLICIT NONE
  TYPE(NCFTIME), POINTER :: FTM
  TYPE(TIME)   :: NOW
  INTEGER :: STATUS

  REAL(SP) :: DF
  INTEGER ::  FRAME
  TYPE(TIME) :: CTIME

  IF (DBG_SET(DBG_SBRIO))THEN
     WRITE(IPT,*) "%%%%%%%%%%%% RECURSIVE_FILE_BRACK INPUT %%%%%%%%%%%%%%" 
     CALL PRINT_FTIME(FTM)
     WRITE(IPT,*) "%%%%%%%%%%%%"
     CALL PRINT_REAL_TIME(NOW,IPT,"NOW")
     WRITE(IPT,*) "%%%%%%%%%%%% ========================== %%%%%%%%%%%%%%" 
     
  END IF

  ! DO SOME ERROR CHECKING FOR BAD VALUES
  IF (FTM%NEXT_STKCNT == FTM%PREV_STKCNT) THEN
     
     CALL FATAL_ERROR("PREVIOUS STKCNT  IS EQUAL TO NEXT!")

  ELSE IF (FTM%NEXT_STKCNT < FTM%PREV_STKCNT) THEN
     
     CALL FATAL_ERROR("PREVIOUS STKCNT IS GREATER THAN NEXT!")
     
  END IF

  IF ( FTM%PREV_IO == FTM%NEXT_IO ) THEN
     
     CALL FATAL_ERROR("IT SEEMS YOUR FILE HAS DUPLICATE TIME VALUES!")
  
  ELSEIF ( FTM%PREV_IO > FTM%NEXT_IO ) THEN
     
     CALL FATAL_ERROR("IT SEEMS YOUR FILE HAS NONE MONOTONIC TIME!")
     
  END IF
  

  ! DECIDE WHAT TO DO BASED ON THE CURRENT BRACKET!

  ! IF PREV AND NEXT BRACK THE CURRENT VALUE
  IF( FTM%PREV_IO < NOW .AND. NOW < FTM%NEXT_IO ) THEN

!     IF (dbg_set(dbg_log)) write(IPT,*) "PREVIOUS AND NEXT BRACKET NOW"

     IF(FTM%NEXT_STKCNT == FTM%PREV_STKCNT+1) THEN
!        IF (dbg_set(dbg_log)) write(IPT,*) "PREVIOUS AND NEXT ARE CONTIGUOUS FRAMES"

        STATUS = 0
        RETURN ! WE FOUND IT
     END IF

     ! CUT THE FRAME DIFFERENCE IN HALF AND TRY AGAIN
     DF = FTM%NEXT_STKCNT - FTM%PREV_STKCNT
     ! NOTE: DF WILL ALWAYS BE GREATER THAN OR EQUAL TO 2
     FRAME  = FTM%PREV_STKCNT + CEILING(DF/2.0_SP)
     CTIME = GET_FILE_TIME(FTM,FRAME)
     
     IF( CTIME < NOW) THEN
        FTM%PREV_STKCNT = FRAME
        FTM%PREV_IO = CTIME
        
        CALL RECURSIVE_FILE_BRACKET(FTM,NOW,status)
     ELSE IF (CTIME > NOW) THEN
        FTM%NEXT_STKCNT = FRAME
        FTM%NEXT_IO = CTIME

        CALL RECURSIVE_FILE_BRACKET(FTM,NOW,status)
     ELSE IF (CTIME == NOW) THEN

        ! FRAM CAN NEVER BE EQUAL TO ONE
        ! NO NEED TO CHECK
        
        FTM%PREV_STKCNT = FRAME-1
        FTM%PREV_IO = GET_FILE_TIME(FTM,FRAME-1)
        
        FTM%NEXT_STKCNT = FRAME
        FTM%NEXT_IO = CTIME
        
        status = 0
        RETURN
                
     ELSE
        CALL FATAL_ERROR("YOU SHOULD NOT BE HERE - I MEAN IT!")
     END IF
     
  ELSE IF( FTM%PREV_IO == NOW ) THEN

!     IF (dbg_set(dbg_log)) write(IPT,*) "PREVIOUS EQUAL TO NOW"
     ! WE HAVE THE ANSWER
     FTM%NEXT_STKCNT = FTM%PREV_STKCNT + 1
     FTM%NEXT_IO = GET_FILE_TIME(FTM,FTM%NEXT_STKCNT)
     
     status = 0
     RETURN
     
  ELSE IF( FTM%NEXT_IO == NOW ) THEN
     
!     IF (dbg_set(dbg_log)) write(IPT,*) "NEXT EQUAL TO NOW"
     ! WE HAVE THE ANSWER
     FTM%PREV_STKCNT = FTM%NEXT_STKCNT - 1
     FTM%PREV_IO = GET_FILE_TIME(FTM,FTM%PREV_STKCNT)
     
     status = 0
     RETURN
     
  ELSE IF( NOW > FTM%NEXT_IO ) THEN

!     IF (dbg_set(dbg_log)) write(IPT,*) "NOW GREATER THAN NEXT"
  
     IF (FTM%NEXT_STKCNT .GE. FTM%STK_LEN) THEN
        ! THE TIME EXCEEDS THE FILE TIME
        status = 1
        return
     END IF
 
     DF = (FTM%STK_LEN - FTM%NEXT_STKCNT)
     FRAME  = FTM%NEXT_STKCNT + CEILING(DF/2.0_SP)
     CTIME = GET_FILE_TIME(FTM,FRAME)

     ! IS THE RESULTANT TIME LESS THAN EQUAL TO OR GREATER THAN NOW?
     IF( CTIME < NOW) THEN
        FTM%PREV_STKCNT = FRAME
        FTM%PREV_IO = CTIME

        FTM%NEXT_STKCNT = FTM%STK_LEN
        FTM%NEXT_IO = GET_FILE_TIME(FTM,FTM%STK_LEN)
        CALL RECURSIVE_FILE_BRACKET(FTM,NOW,status)

     ELSE IF (CTIME > NOW) THEN

        FTM%PREV_STKCNT = FTM%NEXT_STKCNT
        FTM%PREV_IO = FTM%NEXT_IO

        FTM%NEXT_STKCNT = FRAME
        FTM%NEXT_IO = CTIME
        
        CALL RECURSIVE_FILE_BRACKET(FTM,NOW,status)

     ELSE IF (CTIME == NOW) THEN
       ! WE HAVE THE ANSWER
        
        FTM%PREV_STKCNT = FRAME-1
        FTM%PREV_IO = GET_FILE_TIME(FTM,FRAME-1)
        
        FTM%NEXT_STKCNT = FRAME
        FTM%NEXT_IO = CTIME
        
        status = 0
        RETURN
        
        
     ELSE
        CALL FATAL_ERROR("YOU SHOULD NOT BE HERE - YOU DON'T LIKE ME, DO YOU!")
     END IF

     
      
  ELSE IF ( FTM%PREV_IO > NOW) THEN

!     IF (dbg_set(dbg_log)) write(IPT,*) "PREVIOUS GREATER THAN NOW"

     IF (FTM%PREV_STKCNT .LE. 1) THEN
        ! THE TIME PRECEEDS THE FILE TIME
        status = -1
        return
     END IF

     DF = (FTM%PREV_STKCNT)
     FRAME  = CEILING(DF/2.0_SP)
     CTIME = GET_FILE_TIME(FTM,FRAME)

     ! IS THE RESULTANT TIME LESS THAN EQUAL TO OR GREATER THAN NOW?
     IF( CTIME < NOW) THEN

        FTM%NEXT_STKCNT = FTM%PREV_STKCNT
        FTM%NEXT_IO = FTM%PREV_IO

        FTM%PREV_STKCNT = FRAME
        FTM%PREV_IO = CTIME
        CALL RECURSIVE_FILE_BRACKET(FTM,NOW,status)

     ELSE IF (CTIME > NOW) THEN

        FTM%PREV_STKCNT = 1
        FTM%PREV_IO = GET_FILE_TIME(FTM,1)

        FTM%NEXT_STKCNT = FRAME
        FTM%NEXT_IO = CTIME
        
        CALL RECURSIVE_FILE_BRACKET(FTM,NOW,status)

     ELSE IF (CTIME == NOW) THEN
        ! WE HAVE THE ANSWER
        FTM%PREV_STKCNT = FRAME
        FTM%PREV_IO = CTIME
        
        FTM%NEXT_STKCNT = FRAME+1
        FTM%NEXT_IO = GET_FILE_TIME(FTM,FRAME+1)
        
        status = 0
        RETURN
     
     ELSE

        CALL FATAL_ERROR("YOU SHOULD NOT BE HERE - WHY OH WHY")
     END IF

  ELSE

     CALL FATAL_ERROR("YOU SHOULD NOT BE HERE - THIS ONE IS NOT GOOD EITHER!")
  END IF
  


END SUBROUTINE RECURSIVE_FILE_BRACKET
!====================================================================
!====================================================================
SUBROUTINE UPDATE_FILE_BRACKET(NCF,NOW,status)
  IMPLICIT NONE
  TYPE(NCFILE), POINTER:: NCF
  TYPE(NCFTIME), POINTER :: FTM
  TYPE(TIME)   :: NOW

  TYPE(TIME)   :: TIMETEST,dtime
  REAL(DP)     :: denom, numer
  INTEGER      :: STATUS

  IF(.NOT. ASSOCIATED(NCF)) CALL FATAL_ERROR &
       & ("UPDATE_FILE_BRACKET: FILE object argument is not assocaited!")

  IF (.NOT. ASSOCIATED(NCF%FTIME)) THEN
     CALL PRINT_FILE(NCF)
     CALL FATAL_ERROR("UPDATE_FILE_BRACKET: FILE object's FTIME is not assocaited!")
  END IF

  FTM => NCF%FTIME

  IF (FTM%STK_LEN == 1) THEN
     CALL PRINT_FILE(NCF)
     CALL FATAL_ERROR ("FILE BRACKET DOES NOT WORK IF THE TIME DIMENSI&
          &ON LENGTH IS ONE!")
  END IF

  IF (FTM%NEXT_STKCNT == FTM%PREV_STKCNT .or. &
       & FTM%PREV_IO == FTM%NEXT_IO ) THEN

     FTM%NEXT_STKCNT = FTM%STK_LEN
     FTM%NEXT_IO = GET_FILE_TIME(FTM,FTM%STK_LEN)


     FTM%PREV_STKCNT = 1
     FTM%PREV_IO = GET_FILE_TIME(FTM,1)

     IF(NOW >  FTM%NEXT_IO) THEN
        status = 1
        RETURN
     END IF


     IF( NOW < FTM%PREV_IO) THEN
        status = -1
        RETURN
     END IF

     CALL RECURSIVE_FILE_BRACKET(FTM,NOW,STATUS)

  END IF


  ! SET STATUS BASED ON FINDING THE CORRECT FILE TIMES 
  IF( FTM%PREV_IO < NOW .AND. NOW <= FTM%NEXT_IO .AND. &
       & FTM%NEXT_STKCNT .EQ. FTM%PREV_STKCNT+1 ) THEN

     status = 0
     ! UPDATE THE WEIGHTS AND EXIT

  ELSE IF( FTM%PREV_IO <= NOW .AND. NOW < FTM%NEXT_IO .AND. &
       & FTM%NEXT_STKCNT .EQ. FTM%PREV_STKCNT+1 ) THEN

     status = 0
     ! UPDATE THE WEIGHTS AND EXIT

  ELSE IF( NOW > FTM%NEXT_IO ) THEN

     !TRY ADVANCING STACK COUNT ONE FIRST

     IF (FTM%NEXT_STKCNT == FTM%STK_LEN) THEN
        STATUS =1
        RETURN
     END IF

     TIMETEST = GET_FILE_TIME(FTM,FTM%NEXT_STKCNT+1)
     IF (TIMETEST >= NOW) THEN
        FTM%PREV_STKCNT = FTM%NEXT_STKCNT
        FTM%PREV_IO = FTM%NEXT_IO

        FTM%NEXT_STKCNT = FTM%NEXT_STKCNT +1
        FTM%NEXT_IO = TIMETEST

        status = 0

     ELSE
        CALL RECURSIVE_FILE_BRACKET(FTM,NOW,STATUS)
        if (STATUS /= 0) return

     END IF

  ELSE IF ( FTM%PREV_IO > NOW) THEN

     IF (FTM%PREV_STKCNT == 1) THEN
        STATUS = -1
        RETURN
     END IF

     TIMETEST = GET_FILE_TIME(FTM,FTM%PREV_STKCNT-1)
     IF (TIMETEST <= NOW) THEN
        FTM%NEXT_STKCNT = FTM%PREV_STKCNT
        FTM%NEXT_IO = FTM%PREV_IO

        FTM%PREV_STKCNT = FTM%PREV_STKCNT-1
        FTM%PREV_IO = TIMETEST

        status = 0

     ELSE

        CALL RECURSIVE_FILE_BRACKET(FTM,NOW,STATUS)
        if (STATUS /= 0) return

     END IF

  ELSE
     CALL FATAL_ERROR &
          & ("And you may ask yourself", &
          & "How do I work this?", &
          & "And you may ask yourself" , &
          & "Where is that fvcom manual? - The Talking Heads")
  END IF

  ! NOW SET THE WGHT VALUES

  ! CALCULATE THE TIME DIFFERENCE IN MICROSECONDS AND CONVERT TO
  ! DOUBLE THEN DIVIDE BY ONE MILLION 
  !  dtime = NOW - FTM%PREV_IO
  !  NUMER = REAL((dtime%MuSod + MUSPD * dtime%MJD), DP)/1000000.0_DP

  !  NUMER = REAL_TIME_DIFF(NOW,FTM%PREV_IO)
  NUMER = SECONDS(NOW - FTM%PREV_IO)

  !  dtime = FTM%NEXT_IO - FTM%PREV_IO
  !  DENOM = REAL((dtime%MuSod + MUSPD * dtime%MJD), DP)/1000000.0_DP

  !  DENOM = REAL_TIME_DIFF(FTM%NEXT_IO,FTM%PREV_IO)
  DENOM = SECONDS(FTM%NEXT_IO - FTM%PREV_IO)


  ! TAKE THE RATIO IN DOUBLE PRECISION AND CONVERT IF MODEL IS NOT DOUBLE
  FTM%NEXT_WGHT = NUMER/DENOM

  FTM%PREV_WGHT = 1.0_DP - NUMER/DENOM


END SUBROUTINE UPDATE_FILE_BRACKET
  !==============================================================================|
  !==============================================================================|
  SUBROUTINE UPDATE_VAR_BRACKET(NCF,VPREV,VNEXT,NOW,status,INTERP)
    !
    !  RETURN STATUS VALUES:
    !     -1   NOW is before the first forcing time
    !      0   The data is current
    !      1   NOW is after the last forcing time
    !
    !
    IMPLICIT NONE
    TYPE(NCVAR), POINTER :: VNEXT,VPREV, VTMP
    TYPE(TIME) :: NOW
    TYPE(NCFILE), POINTER :: NCF
    INTEGER :: STATUS
    TYPE(INTERP_WEIGHTS),POINTER, OPTIONAL :: INTERP

    REAL(SP), POINTER :: VARRP(:,:),VVECP(:)
    TYPE(NCFTIME),POINTER :: FTM

    LOGICAL :: FOUND


    IF(.not. ASSOCIATED(NCF)) CALL FATAL_ERROR&
         & ("UPDATE_VAR_BRACKET: FILE OBJECT ARGUMENT IS NOT ASSOCIATED!")

    IF(.not. ASSOCIATED(VNEXT)) CALL FATAL_ERROR&
         & ("UPDATE_VAR_BRACKET: FIRST VARIABLE ARGUMENT IS NOT ASSOCIATED!")

    IF(.not. ASSOCIATED(VPREV)) CALL FATAL_ERROR&
         & ("UPDATE_VAR_BRACKET: SECOND VARIABLE ARGUMENT IS NOT ASSOCIATED!")

    IF(PRESENT(INTERP)) THEN
       IF (.not. ASSOCIATED(INTERP)) CALL FATAL_ERROR&
            & ("UPDATE_VAR_BRACKET: THE INTERP ARGUMENT IS NOT ASSOCIATED!")
    END IF

    NULLIFY(VTMP)


    IF(.not.ASSOCIATED(VPREV%ncid, target = VNEXT%ncid)) THEN
       call print_var(VPREV)
       call print_var(VNEXT)
       CALL FATAL_ERROR ("UPDATE_VAR_BRACKET: ", &
            &" VARIABLE ARGUMENTS DO NOT POINT TO THE SAME NETCDF FILE!")
    END IF

    IF(.not.ASSOCIATED(VPREV%ncid, target = NCF%ncid)) THEN
       call print_var(VPREV)
       call print_file(NCF)
       CALL FATAL_ERROR ("UPDATE_VAR_BRACKET: ", &
            &" VARIABLE ARGUMENTS DOES NOT POINT TO THE NETCDF FILE!")
    END IF




    CALL UPDATE_FILE_BRACKET(NCF,NOW,status)
    if (status /= 0) RETURN


    FTM => NCF%FTIME

    IF(FTM%NEXT_STKCNT .EQ. VNEXT%CURR_STKCNT) THEN

       IF(FTM%PREV_STKCNT .NE. VPREV%CURR_STKCNT) THEN
          CALL NC_READ_VAR(VPREV,FTM%PREV_STKCNT)

          IF(PRESENT(INTERP)) THEN
             CALL NC_POINT_VAR(VPREV,VVECP)
             CALL NC_POINT_VAR(VPREV,VARRP)              
             CALL INTERP_BILINEAR_P(VARRP,INTERP,VVECP)
          END IF

          !        status = -2
       END IF

    ELSE IF (FTM%PREV_STKCNT .EQ. VPREV%CURR_STKCNT) THEN
       ! NEXT MUST NEED UPDATE
       CALL NC_READ_VAR(VNEXT,FTM%NEXT_STKCNT)
       IF(PRESENT(INTERP)) THEN
          CALL NC_POINT_VAR(VNEXT,VVECP)
          CALL NC_POINT_VAR(VNEXT,VARRP)              
          CALL INTERP_BILINEAR_P(VARRP,INTERP,VVECP)
       END IF
       !     status = 2

    ELSE  ! IF THEY ARE BOTH WRONG: SEE IF WE CAN JUST INCRIMENT BY ONE CNT?

       IF(FTM%PREV_STKCNT .EQ. VNEXT%CURR_STKCNT) THEN
          VTMP => VPREV
          VPREV => VNEXT
          VNEXT => VTMP
          CALL NC_READ_VAR(VNEXT,FTM%NEXT_STKCNT)
          IF(PRESENT(INTERP)) THEN
             CALL NC_POINT_VAR(VNEXT,VVECP)
             CALL NC_POINT_VAR(VNEXT,VARRP)              
             CALL INTERP_BILINEAR_P(VARRP,INTERP,VVECP)
          END IF
          !        status = 2
       ELSE IF ( FTM%NEXT_STKCNT .EQ. VPREV%CURR_STKCNT) THEN
          VTMP => VNEXT
          VNEXT => VPREV
          VPREV => VTMP
          CALL NC_READ_VAR(VPREV,FTM%PREV_STKCNT)
          IF(PRESENT(INTERP)) THEN
             CALL NC_POINT_VAR(VPREV,VVECP)
             CALL NC_POINT_VAR(VPREV,VARRP) 
             CALL INTERP_BILINEAR_P(VARRP,INTERP,VVECP)
          END IF
       ELSE ! MUST LOAD BOTH FROM THE FILE

          CALL NC_READ_VAR(VNEXT,FTM%NEXT_STKCNT)

          IF(PRESENT(INTERP)) THEN
             CALL NC_POINT_VAR(VNEXT,VVECP)
             CALL NC_POINT_VAR(VNEXT,VARRP)              
             CALL INTERP_BILINEAR_P(VARRP,INTERP,VVECP)
          END IF

          CALL NC_READ_VAR(VPREV,FTM%PREV_STKCNT)

          IF(PRESENT(INTERP)) THEN
             CALL NC_POINT_VAR(VPREV,VVECP)
             CALL NC_POINT_VAR(VPREV,VARRP) 
             CALL INTERP_BILINEAR_P(VARRP,INTERP,VVECP)
          END IF

       END IF
    END IF

  END SUBROUTINE UPDATE_VAR_BRACKET
!====================================================================
!====================================================================
FUNCTION SET_FILE_TIME_TYPE(NCF) RESULT(RES)
! RETURN VALUES:
!  -1 UNKNOWN TIME VARIABLE TYPE
!  -2 NOT AN IDENTIFABLE TIME VARIABLE FILE
!  0 SUCESSES!
  IMPLICIT NONE
  TYPE(NCFILE), pointer :: NCF
  TYPE(NCFTIME), POINTER :: FTM
  INTEGER :: RES
  
  type(time) :: TEST
  Character(len=80) :: tzone
  TYPE(NCDIM), pointer :: DIM
  TYPE(NCATT), pointer :: ATT
  TYPE(NCVAR), pointer :: VAR1, VAR2
  LOGICAL :: FOUND, VALID
  
  IF(.NOT. ASSOCIATED(NCF)) CALL FATAL_ERROR &
       & ("SET_FILE_TIME_TYPE: FILE object argument is not assocaited!")
  
  RES = -1
  

  IF(ASSOCIATED(NCF%FTIME))THEN
     IF(NCF%FTIME%TMTYPE /= TMtype_UNKNOWN) THEN
        RES = 0
        RETURN
     END IF
  END IF

  ! LOOK FOR CHARACTER STRING DATE
  VAR1  => FIND_VAR(NCF,'Times',FOUND)
  IF (FOUND) THEN

     VALID = IS_VALID_DATETIME(VAR1,tzone)
     IF (VALID) THEN

        ! WE HAVE AN IDENTIFIED TIME VARIABLE FILE 
        IF(.NOT. ASSOCIATED(NCF%FTIME)) NCF%FTIME=>NEW_FTIME()


        ! SET THE LEGNTH OF THE TIME DIMENSION
        DIM => FIND_DIM(VAR1,'Time',FOUND)
        IF (.not. FOUND) DIM => FIND_DIM(VAR1,'time',FOUND)
        
        IF (.not. FOUND) THEN
           CALL PRINT_FILE(NCF)
           CALL PRINT_DIM_LIST(NCF)
           CALL FATAL_ERROR("SET_FILE_TIME_TYPE: CAN NOT FIND THE TIME DIMENSION")
        END IF
        
        NCF%FTIME%STK_LEN = DIM%DIM
        
        ! COPY THE TIME VARIABLE FOR SPACE TO LOAD IT
        NCF%FTIME%TM1 => REFERENCE_VAR(VAR1)
        
        ! SET THE TYPE
        NCF%FTIME%TMTYPE = TMtype_CHAR_DATE
        NCF%FTIME%TimeZone=TZone

        RES = 0
        
        IF (DIM%DIM == 0) THEN
           CALL WARNING("THE NETCDF FILE"//TRIM(NCF%FNAME),&
                & "Has a Time variable but the dimension is zero!")
           RETURN
        END IF

        if(dbg_set(dbg_io)) write(ipt,*) "Testing get time for file:"//trim(NCF%FNAME)
        TEST=GET_FILE_TIME(NCF,1)
        if(dbg_set(dbg_io)) write(ipt,*) "Test Passed!"
        RETURN
     ELSE
        CALL WARNING ("SET_FILE_TIME_TYPE: FOUND VARIABLE NAMED: 'Times'&
             & BUT IT IS NOT VALID ACCORDING TO IS_VALID_DATETIME")
     END IF
  END IF

  
  ! LOOK FOR 2 INTEGER MODIFIED JULIAN DAY SECOND
  VAR1  => FIND_VAR(NCF,'Itime',FOUND)
  IF (FOUND) THEN
     VAR2 => FIND_VAR(NCF,'Itime2',FOUND)
     
     IF (FOUND) THEN
        VALID = IS_VALID_ITIME( VAR1, VAR2,TZone)

        IF (VALID) THEN
           ! WE HAVE A TIME VARIABLE FILE 
           IF(.NOT. ASSOCIATED(NCF%FTIME)) NCF%FTIME=>NEW_FTIME()
    
           ! SET THE LEGNTH OF THE TIME DIMENSION
           DIM => FIND_DIM(VAR1,'time',FOUND)           
           IF (FOUND) THEN
              DIM => FIND_DIM(VAR2,'time',FOUND) ! CHECK BOTH VARIABLES
              IF (.not. FOUND) THEN
                 CALL PRINT_FILE(NCF)
                 CALL PRINT_DIM_LIST(NCF)
                 CALL FATAL_ERROR("SET_FILE_TIME_TYPE: CAN NOT FIND THE TIME DIMENSION")
              END IF
           END IF
           NCF%FTIME%STK_LEN = DIM%DIM

           ! COPY THE TIME VARIABLE FOR SPACE TO LOAD IT
           NCF%FTIME%TM1 => REFERENCE_VAR(VAR1)
           NCF%FTIME%TM2 => REFERENCE_VAR(VAR2)

           ! SET THE TYPE
           NCF%FTIME%TMTYPE = TMtype_INT2_MJD 

           NCF%FTIME%TimeZone=TZone

           RES = 0

           IF (DIM%DIM == 0) THEN
              CALL WARNING("THE NETCDF FILE"//TRIM(NCF%FNAME),&
                   & "Has a Time variable but the dimension is zero!")
              RETURN
           END IF

           if(dbg_set(dbg_io)) write(ipt,*) "Testing get time for file:"//trim(NCF%FNAME)
           TEST=GET_FILE_TIME(NCF,1)
           if(dbg_set(dbg_io)) write(ipt,*) "Test Passed!"
           RETURN
        ELSE
           CALL WARNING ("SET_FILE_TIME_TYPE: FOUND VARIABLE NAMES: Itime&
                & and Itime2 BUT THEY ARE NOT VALID ACCORDING TO IS_VALID_ITIME")
        END IF

     ELSE 
        CALL WARNING ("SET_FILE_TIME_TYPE: FOUND ONLY ONE OF TWO INTEGER TIME VARIABLES?")
     END IF
  END IF
       
  ! LOOK FOR FLOATING POINT MODIFIED JULIAN DAY
  VAR1 => FIND_VAR(NCF,'time',FOUND)
  IF (FOUND) THEN
     
     IF (IS_VALID_FLOAT_DAYS( VAR1, TZone) ) THEN
        ! WE HAVE A TIME VARIABLE FILE 
        IF(.NOT. ASSOCIATED(NCF%FTIME)) NCF%FTIME=>NEW_FTIME()

        ! SET THE LEGNTH OF THE TIME DIMENSION
        DIM => FIND_DIM(VAR1,'time',FOUND)           
        IF (.NOT.FOUND) THEN
           CALL PRINT_FILE(NCF)
           CALL PRINT_DIM_LIST(NCF)
           CALL FATAL_ERROR("SET_FILE_TIME_TYPE: CAN NOT FIND THE TIME DIMENSION")
        END IF
        NCF%FTIME%STK_LEN = DIM%DIM
        
        ! COPY THE TIME VARIABLE FOR SPACE TO LOAD IT
        NCF%FTIME%TM1 => REFERENCE_VAR(VAR1)
        
        NCF%FTIME%TMTYPE = TMtype_FLOAT_DAYS

        NCF%FTIME%TimeZone=TZone

        RES = 0

        IF (DIM%DIM == 0) THEN
           CALL WARNING("THE NETCDF FILE"//TRIM(NCF%FNAME),&
                & "Has a Time variable but the dimension is zero!")
           RETURN
        END IF

        if(dbg_set(dbg_io)) write(ipt,*) "Testing get time for file:"//trim(NCF%FNAME)
        TEST=GET_FILE_TIME(NCF,1)
        if(dbg_set(dbg_io)) write(ipt,*) "Test Passed!"
        RETURN

     ELSE IF(IS_VALID_FLOAT_SECONDS( VAR1, TZone) ) THEN
        
        ! WE HAVE A TIME VARIABLE FILE 
        IF(.NOT. ASSOCIATED(NCF%FTIME)) NCF%FTIME=>NEW_FTIME()
        
        ! SET THE LEGNTH OF THE TIME DIMENSION
        DIM => FIND_DIM(VAR1,'time',FOUND)           
        IF (.NOT.FOUND) THEN
           CALL PRINT_FILE(NCF)
           CALL PRINT_DIM_LIST(NCF)
           CALL FATAL_ERROR("SET_FILE_TIME_TYPE: CAN NOT FIND THE TIME DIMENSION")
        END IF
        NCF%FTIME%STK_LEN = DIM%DIM
        
        !COPY THE TIME VARIABLE FOR SPACE TO LOAD IT
        NCF%FTIME%TM1 => REFERENCE_VAR(VAR1)
        
        NCF%FTIME%TMTYPE = TMtype_FLOAT_SECONDS

        NCF%FTIME%TimeZone=TZone

        RES = 0
        IF (DIM%DIM == 0) THEN
           CALL WARNING("THE NETCDF FILE"//TRIM(NCF%FNAME),&
                & "Has a Time variable but the dimension is zero!")
           RETURN
        END IF


        if(dbg_set(dbg_io)) write(ipt,*) "Testing get time for file:"//trim(NCF%FNAME)
        TEST=GET_FILE_TIME(NCF,1)
        if(dbg_set(dbg_io)) write(ipt,*) "Test Passed!"
        RETURN
     ELSE
        CALL WARNING ("SET_FILE_TIME_TYPE: FOUND VARIABLE NAMEd: 'time'&
             & BUT IT IS NOT VALID ACCORDING TO IS_VALID_FLOAT_DAYS/SECONDS")
     END IF
  END IF


END FUNCTION SET_FILE_TIME_TYPE

FUNCTION GET_FILE_TIME_NCF(NCF,STKCNT) RESULT(THETIME)
  IMPLICIT NONE
  TYPE(TIME)   :: THETIME
  TYPE(NCFILE),POINTER :: NCF
  INTEGER, intent(in) :: STKCNT
  TYPE(NCFTIME), POINTER :: FTM

   IF(.NOT. ASSOCIATED(NCF)) CALL FATAL_ERROR &
       & ("GET_FILE_TIME_NCF: FILE object argument is not assocaited!")
  

   IF(.NOT. ASSOCIATED(NCF%FTIME)) THEN
      CALL PRINT_FILE(NCF)
      CALL FATAL_ERROR("GET_FILE_TIME_NCF: FILE object's FTIME is not assocaited!")
   END IF

   FTM =>NCF%FTIME

   THETIME = GET_FILE_TIME(FTM,STKCNT)

 END FUNCTION GET_FILE_TIME_NCF


FUNCTION GET_FILE_TIME_NCFTIME(FTM,STKCNT) RESULT(THETIME)
  IMPLICIT NONE
  TYPE(TIME)   :: THETIME
  TYPE(NCFTIME),POINTER :: FTM
  TYPE(NCVAR), POINTER :: VAR1, VAR2
  TYPE(NCATT), POINTER :: ATT
  INTEGER, intent(in) :: STKCNT
  REAL(SP),target :: float_time
  INTEGER, target :: mjd,msec
  Character(len=80), target :: dstring
  INTEGER :: STATUS
  Character(len=80) :: dformat
  LOGICAL FOUND

   IF(.NOT. ASSOCIATED(FTM)) THEN
      CALL FATAL_ERROR("GET_FILE_TIME_NCFTIME: THE ARGUMENT FTM is not assocaited!")
   END IF


   IF(STKCNT .LE. 0) CALL FATAL_ERROR&
        &('GET_FILE_TIME: YOU CAN NOT GET TIME FOR A STACK COUNT LESS THAN ONE!')

   SELECT CASE(FTM%TMTYPE)
   CASE(TMtype_FLOAT_DAYS)

      var1 => FTM%TM1
      
      CALL NC_CONNECT_AVAR(VAR1,float_time)
      CALL NC_READ_VAR(VAR1,stkcnt)

      THETIME = DAYS2TIME(float_time) - TIME_ZONE(FTM%TIMEZONE,status)

   CASE(TMtype_FLOAT_SECONDS)

      var1 => FTM%TM1
      
      CALL NC_CONNECT_AVAR(VAR1,float_time)
      CALL NC_READ_VAR(VAR1,stkcnt)

      THETIME = SECONDS2TIME(float_time) - TIME_ZONE(FTM%TIMEZONE,status)

   CASE(TMtype_INT2_MJD)

      var1 => FTM%TM1
      
      CALL NC_CONNECT_AVAR(VAR1,mjd)
      CALL NC_READ_VAR(VAR1,stkcnt)

      var2 => FTM%TM2

      CALL NC_CONNECT_AVAR(VAR2,msec)
      CALL NC_READ_VAR(VAR2,stkcnt)

      THETIME%MJD= mjd
      THETIME%MuSod= INT(msec,ITIME)* INT(1000,ITIME)
 
      THETIME = THETIME - TIME_ZONE(FTM%TIMEZONE,status)

   CASE(TMtype_CHAR_DATE)

      dformat = "ymd"

      var1 => FTM%TM1
      
      CALL NC_CONNECT_AVAR(VAR1,dstring)
      CALL NC_READ_VAR(VAR1,stkcnt)

      THETIME = READ_DATETIME(dstring,dFORMAT,FTM%TIMEZONE,status)
      if(STATUS == 0) THEN
         CALL PRINT_FTIME(FTM)
         CALL FATAL_ERROR("GET_FILE_TIME: COULD NOT TRANSLATE TIME STRING",&
              "FILE RETURNED: "//TRIM(dstring))
      END if


   CASE DEFAULT
      CALL PRINT_FTIME(FTM)
      CALL FATAL_ERROR("GET_FILE_TIME: FILE object's FTIME is an unrecognized type")
   END SELECT

  

 END FUNCTION GET_FILE_TIME_NCFTIME
 
 !========================================================================!
 !
 !========================================================================!
 SUBROUTINE Set_File_Stack(NCF_IN,NOW,INTERVAL)
   ! THIS IS A STRICT FILE TIME MATCHING PROGRAM. IT WILL MATCH A
   ! TIME IN A EXISTING FILE TO THE NCF WITH THE SAME FILE NAME.
   ! IF THE EXACT TIME PASSED DOES NOT EXIST IN THE FILE, IT RETURNS
   ! AN ERROR.
   IMPLICIT NONE
   TYPE(NCFILE), POINTER :: NCF_IN,NCF_TMP
   TYPE(TIME) :: Ttest,NOW,INTERVAL
   INTEGER :: IDX
   
   NCF_TMP => NEW_FILE()
   NCF_TMP%FNAME = NCF_IN%FNAME
   
   Call NC_OPEN(NCF_TMP)
   CALL NC_LOAD(NCF_TMP)
   
   IDX = NCF_TMP%FTIME%STK_LEN
   
   IF(IDX == 0) CALL FATAL_ERROR&
        &("THERE IS NO DATA IN THE FILE:"//TRIM(NCF_IN%FNAME),&
        & "(STK_LEN == 0) WHEN THE MODEL CRASHED, SO HOTSTART IS IMPOSSIBLE.",&
        & "BETTER LUCK NEXT TIME!")
   
   Ttest = GET_FILE_TIME(NCF_TMP,IDX)
   
   DO WHILE(Ttest > NOW)
      
      IDX = IDX -1
      IF(IDX ==0) THEN
         CALL PRINT_TIME(NOW,IPT,"NOW")
         CALL PRINT_TIME(TTEST,IPT,"EARLIEST FILE TIME")
         
         CALL FATAL_ERROR&
              &("SET_FILE_STACK: SEARCHING FOR TIME IN FILE NAME:",&
              & TRIM(NCF_IN%FNAME),&
              & "ALL TIMES IN FILE ARE GREATER THAN NOW?")
      END IF
      
      Ttest = GET_FILE_TIME(NCF_TMP,IDX)
      
   END DO
   
   IF(TTEST /= NOW) THEN
      CALL PRINT_TIME(NOW,IPT,"NOW")
      CALL PRINT_TIME(TTEST,IPT,"NEAREST TIME IN FILE")
      
      CALL FATAL_ERROR&
           &("SET_FILE_STACK: SEARCHING FOR TIME IN FILE NAME:",&
           & TRIM(NCF_IN%FNAME),&
           & "CAN'T FIND EXACT MATCHING TIME IN FILE")
   END IF
   
   CALL KILL_FILE(NCF_TMP)
   
   NCF_IN%FTIME%PREV_STKCNT=IDX
   NCF_IN%FTIME%NEXT_STKCNT=IDX+1
   NCF_IN%FTIME%PREV_IO=TTEST
   NCF_IN%FTIME%NEXT_IO=TTEST+INTERVAL
   
   
 END SUBROUTINE Set_File_Stack
 
 

END MODULE MOD_NCTOOLS

