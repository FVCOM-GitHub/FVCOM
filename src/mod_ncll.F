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

MODULE MOD_NCLL ! NETCDF FILE CROSS LINKED LIST
  USE NETCDF
  USE MOD_TIME
  USE MOD_UTILS
  USE MOD_INTERP
  IMPLICIT NONE


  INTEGER, PARAMETER :: TMtype_UNKNOWN = 0
  INTEGER, PARAMETER :: TMtype_CHAR_DATE = 1
  INTEGER, PARAMETER :: TMtype_INT2_MJD  = 2
  INTEGER, PARAMETER :: TMtype_FLOAT_DAYS  = 3
  INTEGER, PARAMETER :: TMtype_FLOAT_SECONDS  = 4



  TYPE NCDIM
     INTEGER :: DIMID
     CHARACTER(Len=NF90_MAX_NAME+1) :: DIMNAME
     INTEGER :: DIM
     LOGICAL :: UNLIMITED
!     LOGICAL :: POINTING    ! To test if two DIM POINTERS go to the same memory
  END TYPE NCDIM

  TYPE NCDIMP
     TYPE(NCDIM), POINTER  :: DIM
     TYPE(NCDIMP), POINTER :: NEXT
  END TYPE NCDIMP

  integer, parameter :: Char_max_attlen = 160
  TYPE NCATT
     CHARACTER(Len=NF90_MAX_NAME+1)  ATTNAME
     INTEGER    XTYPE
     INTEGER    LEN
     INTEGER    ATTID
     INTEGER,   ALLOCATABLE, DIMENSION(:) :: int
     REAL(SPA), ALLOCATABLE, DIMENSION(:) :: flt
     REAL(DP), ALLOCATABLE, DIMENSION(:) :: dbl
     CHARACTER(LEN=char_max_attlen), ALLOCATABLE, DIMENSION(:) :: CHR
  END TYPE NCATT
  
  TYPE NCATTP
     TYPE(NCATT), POINTER  :: ATT
     TYPE(NCATTP), POINTER :: NEXT
  END TYPE NCATTP
  

  ! NCVAR IS THE HEAD NODE FOR THE LIST: DIMS, ATTS
  TYPE NCVAR
     LOGICAL CONNECTED
     INTEGER, POINTER :: NCID ! THIS IS NEVER ALLOCATED, ONLY ASSOCIATED
     INTEGER VARID
     CHARACTER(Len=NF90_MAX_NAME+1) VARNAME
     INTEGER XTYPE
     ! THE DIMIDS ARE DETERMINED BY WHICH DIMS ARE ASSOCIATED AND IN
     ! WHAT ORDER - THERE IS NO EXPLICIT DIMIDS VARIABLE
!     INTEGER, allocatable :: dimids(:)
     TYPE(NCDIMP), POINTER :: DIMS
     TYPE(NCATTP), POINTER :: ATTS

     ! Data pointers
     INTEGER :: CURR_STKCNT
     
     ! THE DATA POINTERS ARE NEVER ALLOCATED, ONLY ASSOCIATED!

     INTEGER, POINTER                        :: SCL_INT
     INTEGER, POINTER,DIMENSION(:)           :: VEC_INT
     INTEGER, POINTER,DIMENSION(:,:)         :: ARR_INT
     INTEGER, POINTER,DIMENSION(:,:,:)       :: CUB_INT
     INTEGER, POINTER,DIMENSION(:,:,:,:)     :: FDA_INT

     REAL(SPA), POINTER                      :: SCL_FLT
     REAL(SPA), POINTER,DIMENSION(:)         :: VEC_FLT
     REAL(SPA), POINTER,DIMENSION(:,:)       :: ARR_FLT
     REAL(SPA), POINTER,DIMENSION(:,:,:)     :: CUB_FLT
     REAL(SPA), POINTER,DIMENSION(:,:,:,:)   :: FDA_FLT

     REAL(DP), POINTER                       :: SCL_DBL
     REAL(DP), POINTER,DIMENSION(:)          :: VEC_DBL
     REAL(DP), POINTER,DIMENSION(:,:)        :: ARR_DBL
     REAL(DP), POINTER,DIMENSION(:,:,:)      :: CUB_DBL
     REAL(DP), POINTER,DIMENSION(:,:,:,:)    :: FDA_DBL

     CHARACTER(LEN=80), POINTER              :: SCL_CHR
     CHARACTER(LEN=80), POINTER,DIMENSION(:) :: VEC_CHR
  END TYPE NCVAR
  
  TYPE NCVARP
     TYPE(NCVAR), POINTER :: VAR
     TYPE(NCVARP), POINTER :: NEXT
  END TYPE NCVARP
  
  TYPE NCFTIME
     INTEGER :: TMTYPE=0
     TYPE(NCVAR), POINTER :: TM1
     TYPE(NCVAR), POINTER :: TM2

     INTEGER :: STK_LEN=0

     INTEGER    ::  PREV_STKCNT=0
     INTEGER    ::  NEXT_STKCNT=0
     INTEGER    ::  MAX_STKCNT=0
     TYPE(TIME) :: PREV_IO
     TYPE(TIME) :: NEXT_IO  
     REAL(SP)   :: PREV_WGHT=0.0_SP
     REAL(SP)   :: NEXT_WGHT=0.0_SP
     TYPE(TIME) :: INTERVAL
     CHARACTER(len=80)  :: TIMEZONE ="none"

  END TYPE NCFTIME


  ! NCFILE IS THE FIRST HEAD NODE FOR THE LISTS: DIMS,ATTS,VARS
  TYPE NCFILE
     INTEGER, POINTER :: NCID ! EACH VARIABLE NEED TO POINT AT THIS!
     CHARACTER(LEN=160) FNAME
     LOGICAL WRITABLE
     LOGICAL OPEN
     LOGICAL CONNECTED
     LOGICAL INDEFMODE
     INTEGER :: unlimdimid
     
     TYPE(NCDIMP), POINTER :: DIMS
     TYPE(NCATTP), POINTER :: ATTS
     TYPE(NCVARP), POINTER :: VARS

     TYPE(NCFTIME), POINTER :: FTIME
     TYPE(INTERP_WEIGHTS),POINTER :: INTERP_N
     TYPE(INTERP_WEIGHTS),POINTER :: INTERP_C
  END TYPE NCFILE
    
  TYPE NCFILEP
     TYPE(NCFILE), pointer :: NCF
     TYPE(NCFILEP), pointer :: NEXT
  END TYPE NCFILEP

  TYPE NCFILELIST
     TYPE(NCFILEP), pointer :: FIRST
  END TYPE NCFILELIST
  

!DO NOT OVERLOAD THE OPERATORS - BE EXPLICIT ABOUT HOW DIMS,VARS,ATTS AND FILES ARE FOUND
!!$  INTERFACE OPERATOR(==)
!!$     MODULE PROCEDURE NAME_EQ_NCF
!!$     MODULE PROCEDURE NCF_EQ_NAME
!!$  END INTERFACE


!!!!!!
!!!!!!  list of module procedures
!!!!!!
! NEW_FILE
! NEW_FILE_LIST
  INTERFACE DELETE_FILE_LINK
     MODULE PROCEDURE DELETE_FILEP_BYNAME
     MODULE PROCEDURE DELETE_FILEP_BYNCID
  END INTERFACE
  ! DELETE_FILE_LIST
  ! KILL_VAR

  INTERFACE INSERT_FILE_LINK
!     MODULE PROCEDURE INSERT_FILEP_BYNAME
     MODULE PROCEDURE INSERT_FILEP_BYNCF
  END INTERFACE

  INTERFACE FIND_FILE
     MODULE PROCEDURE FIND_FILE_BYNAME
     MODULE PROCEDURE FIND_FILE_BYNCID
  END INTERFACE
! COUNT_FILE_LIST
! PRINT_FILE_LIST


! NEW_VAR
! NEW_VARP
! Reference_Var
! COPY_VAR
  INTERFACE DELETE_VAR_LINK
     MODULE PROCEDURE DELETE_VARP_BYNAME
     MODULE PROCEDURE DELETE_VARP_BYVARID
  END INTERFACE
  ! DELETE_VAR_LIST
  ! KILL_VAR

  INTERFACE INSERT_VAR_LINK
!     MODULE PROCEDURE INSERT_VARP_BYNAME
     MODULE PROCEDURE INSERT_VARP_BYVAR
  END INTERFACE

  INTERFACE FIND_VAR
     MODULE PROCEDURE FIND_VAR_BYNAME
     MODULE PROCEDURE FIND_VAR_BYVARID
  END INTERFACE

! COUNT_VAR_LIST
! PRINT_VAR_LIST
  

! NEW_ATT
! NEW_ATTP
! COPY_ATT
! COPY_ATT_LIST
  INTERFACE DELETE_ATT_LINK
     MODULE PROCEDURE DELETE_NCF_ATTP_BYNAME
     MODULE PROCEDURE DELETE_NCF_ATTP_BYATTID
     MODULE PROCEDURE DELETE_VAR_ATTP_BYNAME
     MODULE PROCEDURE DELETE_VAR_ATTP_BYATTID     
  END INTERFACE

  INTERFACE DELETE_ATT_LIST
     MODULE PROCEDURE DELETE_NCF_ATTP_LIST
     MODULE PROCEDURE DELETE_VAR_ATTP_LIST
  END INTERFACE
  ! KILL_VAR

  INTERFACE INSERT_ATT_LINK
!     MODULE PROCEDURE INSERT_NCF_ATTP_BYNAME
     MODULE PROCEDURE INSERT_NCF_ATTP_BYATT
!     MODULE PROCEDURE INSERT_VAR_ATTP_BYNAME
     MODULE PROCEDURE INSERT_VAR_ATTP_BYATT
  END INTERFACE

  INTERFACE FIND_ATT
     MODULE PROCEDURE FIND_NCF_ATT_BYNAME
     MODULE PROCEDURE FIND_NCF_ATT_BYATTID
     MODULE PROCEDURE FIND_VAR_ATT_BYNAME
     MODULE PROCEDURE FIND_VAR_ATT_BYATTID
  END INTERFACE

  INTERFACE COUNT_ATT_LIST
     MODULE PROCEDURE COUNT_NCF_ATT_LIST
     MODULE PROCEDURE COUNT_VAR_ATT_LIST
  END INTERFACE

  INTERFACE PRINT_ATT_LIST
     MODULE PROCEDURE PRINT_NCF_ATT_LIST
     MODULE PROCEDURE PRINT_VAR_ATT_LIST
  END INTERFACE


! NEW_DIM
! NEW_DIMP
! COPY_DIM
  INTERFACE DELETE_DIM_LINK
     MODULE PROCEDURE DELETE_NCF_DIMP_BYNAME
     MODULE PROCEDURE DELETE_NCF_DIMP_BYDIMID
     MODULE PROCEDURE DELETE_VAR_DIMP_BYNAME
     MODULE PROCEDURE DELETE_VAR_DIMP_BYDIMID     
  END INTERFACE

  INTERFACE DELETE_DIM_LIST
     MODULE PROCEDURE DELETE_NCF_DIMP_LIST
     MODULE PROCEDURE DELETE_VAR_DIMP_LIST
  END INTERFACE
! KILL_DIM

  INTERFACE INSERT_DIM_LINK
!     MODULE PROCEDURE INSERT_NCF_DIMP_BYNAME - NOT IMPLIMENTED- DANGEROUS
     MODULE PROCEDURE INSERT_NCF_DIMP_BYDIM
!     MODULE PROCEDURE INSERT_VAR_DIMP_BYNAME - NOT IMPLIMENTED- DANGEROUS
     MODULE PROCEDURE INSERT_VAR_DIMP_BYDIM
  END INTERFACE

  INTERFACE FIND_DIM
     MODULE PROCEDURE FIND_NCF_DIM_BYNAME
     MODULE PROCEDURE FIND_NCF_DIM_BYDIMID
     MODULE PROCEDURE FIND_VAR_DIM_BYNAME
     MODULE PROCEDURE FIND_VAR_DIM_BYDIMID
  END INTERFACE

  INTERFACE FIND_UNLIMITED
     MODULE PROCEDURE FIND_NCF_DIM_UNLIMITED
     MODULE PROCEDURE FIND_VAR_DIM_UNLIMITED
  END INTERFACE
  
  INTERFACE HAS_UNLIMITED
     MODULE PROCEDURE HAS_UNLIMITED_NCF
     MODULE PROCEDURE HAS_UNLIMITED_VAR
  END INTERFACE

  INTERFACE COUNT_DIM_LIST
     MODULE PROCEDURE COUNT_NCF_DIM_LIST
     MODULE PROCEDURE COUNT_VAR_DIM_LIST
  END INTERFACE

  INTERFACE COUNT_NONSINGLETON_DIM_LIST
     MODULE PROCEDURE COUNT_NCF_NS_DIM_LIST
     MODULE PROCEDURE COUNT_VAR_NS_DIM_LIST
  END INTERFACE

  INTERFACE PRINT_DIM_LIST
     MODULE PROCEDURE PRINT_NCF_DIM_LIST
     MODULE PROCEDURE PRINT_VAR_DIM_LIST
  END INTERFACE


!!!
!!! Note: tried to add interface operator to addition(+) and
!!! subtraction(-) of nc objects but intent(in) is required 
!!! for both arguments of the function. That is not practical. 
!!!
  INTERFACE ADD
     MODULE PROCEDURE VAR_PLUS_DIM
     MODULE PROCEDURE VAR_PLUS_ATT
     MODULE PROCEDURE NCF_PLUS_VAR
     MODULE PROCEDURE NCF_PLUS_ATT
     MODULE PROCEDURE NCF_PLUS_DIM
     MODULE PROCEDURE NCF_PLUS_NCF
     MODULE PROCEDURE NCFLIST_PLUS_NCF
  END INTERFACE
  
  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE COPY_DIM_LIST
     MODULE PROCEDURE COPY_ATT_LIST
     MODULE PROCEDURE COPY_FTIME
  END INTERFACE

CONTAINS
!====================================================================
!====================================================================
!============================ FILES =================================
!====================================================================
!====================================================================
  FUNCTION NEW_FILE(fname) RESULT(NCF)
    IMPLICIT NONE
!    CHARACTER(LEN=160), INTENT(IN)  :: NAME
    TYPE(NCFILE), POINTER :: NCF
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: FNAME
    integer status
    
    nullify(NCF)
    ALLOCATE(NCF,stat=status)
    if(status/=0) CALL FATAL_ERROR("NEW_FILE: COULD NOT ALLOCATE!")
    
    Nullify(NCF%NCID)
    ALLOCATE(NCF%NCID)
    NCF%NCID      = -1
    NCF%FNAME     = ' '
    IF (PRESENT(FNAME)) NCF%FNAME=FNAME
    NCF%WRITABLE = .false.
    NCF%OPEN      = .false.
    NCF%CONNECTED = .false.
    NCF%UNLIMDIMID = -1

    NULLIFY(NCF%FTIME)
    NULLIFY(NCF%INTERP_N)
    NULLIFY(NCF%INTERP_C)

    ! TIME VARIABLES ARE ALL INTEGER => AUTOMATICALLY INITIALIZED TO ZERO

    ! MAKE NEW LISTS FOR THE NCFILE
    NCF%DIMS => NEW_DIMP()
    NCF%ATTS => NEW_ATTP()
    NCF%VARS => NEW_VARP()
  END FUNCTION NEW_FILE
!====================================================================
!====================================================================
  FUNCTION COPY_FILE(NCFIN) RESULT(NCFOUT)
    IMPLICIT NONE
!    CHARACTER(LEN=160), INTENT(IN)  :: NAME
    TYPE(NCFILE), POINTER,INTENT(IN) :: NCFIN
    TYPE(NCFILE), POINTER :: NCFOUT
    TYPE(NCVAR), POINTER :: NCFOUT1
    TYPE(NCVARP), POINTER :: CURRENT
    integer status

    IF(.not. Associated(NCFIN)) Call Fatal_Error &
         ("COPY_NCF: INPUT FILE IS NOT ASSOCIATED!")
    
    NCFOUT => NEW_FILE()
    
    NCFOUT%NCID      = -1
    NCFOUT%FNAME     = NCFIN%FNAME  
    NCFOUT%WRITABLE = .FALSE.
    NCFOUT%OPEN      = .false.
    NCFOUT%CONNECTED = .false.
    NCFOUT%UNLIMDIMID = NCFIN%UNLIMDIMID

!!$
!!$    IF(Associated(NCFIN%INTERP_N)) NCFOUT%INTERP_N => COPY_INTERP(NCFIN%INTERP_N)
!!$
!!$    IF(Associated(NCFIN%INTERP_C)) NCFOUT%INTERP_C => COPY_INTERP(NCFIN%INTERP_C)
!!$

    ! TIME VARIABLES ARE ALL INTEGER => AUTOMATICALLY INITIALIZED TO ZERO

    ! MAKE NEW LISTS FOR THE NCFILE
    IF(.not. Associated(NCFIN%DIMS)) Call Fatal_Error &
         ("COPY_NCF: INPUT FILE DIMS LIST IS NOT ASSOCIATED!")

    IF(.not. Associated(NCFIN%ATTS)) Call Fatal_Error &
         ("COPY_NCF: INPUT FILE ATTS LIST IS NOT ASSOCIATED!")

    IF(.not. Associated(NCFIN%VARS)) Call Fatal_Error &
         ("COPY_NCF: INPUT FILE VARS LIST IS NOT ASSOCIATED!")

    ! USE OVERLOADED ASSIGNMENT OPERATOR
    NCFOUT%DIMS = NCFIN%DIMS
    NCFOUT%ATTS = NCFIN%ATTS

    CURRENT => NCFIN%VARS%NEXT
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) THEN
          EXIT ! END OF VAR LIST
       END IF

       IF(.NOT. ASSOCIATED(CURRENT%VAR)) THEN
          CALL FATAL_ERROR("COPY_FILE: FOUND NULL VAR POINTER IN THE LIST")
       END IF

       NCFOUT1 => COPY_VAR(CURRENT%VAR)
       NCFOUT => ADD(NCFOUT,NCFOUT1)
!!$       NCFOUT => ADD(NCFOUT,COPY_VAR(CURRENT%VAR))

       CURRENT => CURRENT%NEXT

    END DO


    ! NOTE - THIS ONLY SETS THE INTEGER VALUES THE TIME VARIABLE
    ! POINTERS ARE NOT SET HERE!
    IF(Associated(NCFIN%FTIME)) THEN
       NCFOUT%FTIME => NEW_FTIME()
       NCFOUT%FTIME = NCFIN%FTIME
    END IF


  END FUNCTION COPY_FILE
!====================================================================
!====================================================================
  FUNCTION NEW_FTIME() RESULT(FTM)
    IMPLICIT NONE
    TYPE(NCFTIME), POINTER :: FTM
    integer status
    
    ALLOCATE(FTM,stat=status)
    if(status/=0) CALL FATAL_ERROR("NEW_FTM: COULD NOT ALLOCATE FTM!")

    nullify(FTM%tm1)
    nullify(FTM%tm2)

    FTM%TIMEZONE="none"

  END FUNCTION NEW_FTIME
!====================================================================
!====================================================================
  SUBROUTINE COPY_FTIME(FTIME_OUT,FTIME_IN)
    IMPLICIT NONE
    TYPE(NCFTIME),INTENT(OUT):: FTIME_OUT
    TYPE(NCFTIME),INTENT(IN) :: FTIME_IN
    
    FTIME_OUT%TMTYPE = FTIME_IN%TMTYPE
    FTIME_OUT%STK_LEN = FTIME_IN%STK_LEN
    FTIME_OUT%PREV_STKCNT = FTIME_IN%PREV_STKCNT
    FTIME_OUT%NEXT_STKCNT = FTIME_IN%NEXT_STKCNT
    FTIME_OUT%MAX_STKCNT = FTIME_IN%MAX_STKCNT
    FTIME_OUT%PREV_IO = FTIME_IN%PREV_IO
    FTIME_OUT%NEXT_IO = FTIME_IN%NEXT_IO
    FTIME_OUT%PREV_WGHT = FTIME_IN%PREV_WGHT
    FTIME_OUT%NEXT_WGHT = FTIME_IN%NEXT_WGHT
    FTIME_OUT%INTERVAL = FTIME_IN%INTERVAL
    FTIME_OUT%TIMEZONE = FTIME_IN%TIMEZONE

  END SUBROUTINE COPY_FTIME
!====================================================================
!====================================================================
  FUNCTION NEW_FILEP() RESULT(NCFP)
    IMPLICIT NONE
    TYPE(NCFILEP), POINTER :: NCFP
    integer status
    
    ALLOCATE(NCFP,stat=status)
    if(status/=0) CALL FATAL_ERROR("ALLOC_FILEP: COULD NOT ALLOCATE!")
    NULLIFY(NCFP%NEXT)
    NULLIFY(NCFP%NCF)
    

  END FUNCTION NEW_FILEP
!====================================================================
!====================================================================
  FUNCTION NEW_FILEHEAD() RESULT(FILEHEAD)
    IMPLICIT NONE
    TYPE(NCFILELIST), POINTER :: FILEHEAD
    integer status
    
    ALLOCATE(FILEHEAD,stat=status)
    if(status/=0) CALL FATAL_ERROR("ALLOC_FILEHEAD: COULD NOT ALLOCATE!")
    FILEHEAD%FIRST => NEW_FILEP()

  END FUNCTION NEW_FILEHEAD
!====================================================================
!====================================================================
  SUBROUTINE KILL_FILEHEAD(FILEHEAD)
    IMPLICIT NONE
    TYPE(NCFILELIST), POINTER :: FILEHEAD

    call delete_FILE_list(FILEHEAD)
    DEALLOCATE(FILEHEAD)
    NULLIFY(FILEHEAD)
    
  END SUBROUTINE KILL_FILEHEAD
!====================================================================
!====================================================================
  SUBROUTINE DELETE_FILEP_BYNAME(LIST,NAME,FOUND)
    IMPLICIT NONE
    TYPE(NCFILELIST),   INTENT(INOUT) :: LIST
    CHARACTER(LEN=*), INTENT(IN)    :: NAME
    LOGICAL,            INTENT(OUT)   :: FOUND
    TYPE(NCFILEP), pointer            :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%FIRST
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( NAME == CURRENT%NCF%FNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_FILE(CURRENT%NCF)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_FILEP_BYNAME
!====================================================================
!====================================================================
  SUBROUTINE DELETE_FILEP_BYNCID(LIST,NCID,FOUND)
    IMPLICIT NONE
    TYPE(NCFILELIST),   INTENT(INOUT) :: LIST
    INTEGER, INTENT(IN)               :: NCID
    LOGICAL,            INTENT(OUT)   :: FOUND
    TYPE(NCFILEP),pointer             :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%FIRST
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( NCID == CURRENT%NCF%NCID ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_FILE(CURRENT%NCF)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_FILEP_BYNCID
!====================================================================
!====================================================================
  SUBROUTINE DELETE_FILE_LIST(LIST)
    IMPLICIT NONE
    TYPE(NCFILELIST),   INTENT(INOUT) :: LIST
    TYPE(NCFILEP),pointer             :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%FIRST
    CURRENT  => PREVIOUS%NEXT
 
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       PREVIOUS%NEXT => CURRENT%NEXT
       CALL KILL_FILE(CURRENT%NCF)
       DEALLOCATE(CURRENT)
       
       CURRENT => PREVIOUS%NEXT

    END DO
   
  END SUBROUTINE DELETE_FILE_LIST
!====================================================================
!====================================================================
 SUBROUTINE KILL_FILE(NCF)
   TYPE(NCFILE),POINTER :: NCF
   INTEGER STATUS

   IF(.not. ASSOCIATED(NCF))THEN
      
      IF(DBG_SET(DBG_LOG)) CALL WARNING &
           & ("CALLED KILL FILL BUT FILE OBJECT IS NOT ASSOCIATED?")
      RETURN
   END IF

   IF (NCF%OPEN) THEN
      status = nf90_close(NCF%ncid)
   END IF

   CALL DELETE_VAR_LIST(NCF)
   DEALLOCATE(NCF%VARS)

   CALL DELETE_DIM_LIST(NCF)
   DEALLOCATE(NCF%DIMS)

   CALL DELETE_ATT_LIST(NCF)
   DEALLOCATE(NCF%ATTS)
   
   IF(ASSOCIATED(NCF%FTIME)) THEN
      NULLIFY(NCF%FTIME%TM1)
      NULLIFY(NCF%FTIME%TM2)
      DEALLOCATE(NCF%FTIME)
      NULLIFY (NCF%FTIME)
   END IF

   ! INTERP MEMORY IS NEVER BELLONGS TO THE FILE POINTER
   IF(ASSOCIATED(NCF%INTERP_N)) Nullify(NCF%INTERP_N)
   IF(ASSOCIATED(NCF%INTERP_C)) Nullify(NCF%INTERP_C)

   DEALLOCATE(NCF%NCID)
   NULLIFY(NCF%NCID)
   
   DEALLOCATE(NCF,STAT=STATUS)
   IF(STATUS /= 0) CALL FATAL_ERROR("KILL_FILE: COULD NOT DEALLOCATE")
   NULLIFY(NCF)

 END SUBROUTINE KILL_FILE
!====================================================================
!====================================================================
!!$  SUBROUTINE INSERT_FILEP_BYNAME(LIST,NAME,FOUND)
!!$    ! ONLY INSERT NEW FILE IF NOT FOUND
!!$    ! ALWAYS INSERT NCFILEP AT THE END OF THE LIST
!!$    IMPLICIT NONE
!!$    CHARACTER(LEN=160), INTENT(IN)    :: NAME
!!$    LOGICAL,            INTENT(OUT)   :: FOUND
!!$    TYPE(NCFILELIST),   INTENT(INOUT) :: LIST
!!$    TYPE(NCFILEP), pointer            :: CURRENT, PREVIOUS
!!$    
!!$    PREVIOUS => LIST%FIRST
!!$    CURRENT  => PREVIOUS%NEXT
!!$    FOUND = .FALSE.
!!$    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
!!$    DO
!!$       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
!!$       IF( NAME == CURRENT%NCF%FNAME ) THEN ! FOUND THE FILE LINK
!!$          FOUND = .TRUE.
!!$          RETURN
!!$       ELSE
!!$          PREVIOUS => PREVIOUS%NEXT
!!$          CURRENT  => CURRENT%NEXT
!!$       END IF
!!$    END DO
!!$    ! NOT FOUND - ADD NEW FILE TO END OF LIST
!!$    CURRENT => NEW_FILEP()
!!$    CURRENT%NCF => NEW_FILE(NAME)
!!$    
!!$  END SUBROUTINE INSERT_FILEP_BYNAME
!====================================================================
!====================================================================
  SUBROUTINE INSERT_FILEP_BYNCF(LIST,NCF,FOUND)
    ! ONLY INSERT NEW FILE IF NOT FOUND
    ! ALWAYS INSERT NCFILEP AT THE END OF THE LIST
    IMPLICIT NONE
    TYPE(NCFILE),POINTER              :: NCF
    LOGICAL,            INTENT(OUT)   :: FOUND
    TYPE(NCFILELIST),   INTENT(INOUT) :: LIST
    TYPE(NCFILEP), pointer            :: CURRENT, PREVIOUS
    
    IF(.NOT.ASSOCIATED(NCF))&
         & CALL FATAL_ERROR("INSERT_FILEP_BYNCF: NCF NOT ASSOCIATED!")

    PREVIOUS => LIST%FIRST
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.
    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       IF( NCF%FNAME == CURRENT%NCF%FNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          RETURN
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO
    ! NOT FOUND - ADD NEW FILE TO END OF LIST

    PREVIOUS%NEXT => NEW_FILEP()
    PREVIOUS%NEXT%NCF => NCF
    
  END SUBROUTINE INSERT_FILEP_BYNCF
!====================================================================
!====================================================================
  FUNCTION FIND_FILE_BYNAME(LIST,NAME,FOUND) RESULT(NCF)
    IMPLICIT NONE
    TYPE(NCFILE),POINTER               :: NCF
    TYPE(NCFILELIST),    INTENT(IN)    :: LIST
    CHARACTER(LEN=*),  INTENT(IN)    :: NAME
    LOGICAL,             INTENT(OUT)   :: FOUND
    TYPE(NCFILEP), pointer             :: CURRENT, PREVIOUS
    integer                            :: idx

    PREVIOUS => LIST%FIRST
    CURRENT  => PREVIOUS%NEXT
    NULLIFY(NCF)
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       if(.not. associated(current%ncf)) Call Fatal_Error&
            & ("FIND_FILE: Link in file list has file pointer that is &
            &not associated!")

       idx = index(CURRENT%NCF%FNAME,NAME)

       IF( idx /= 0 ) THEN ! FOUND THE FILE LINK
          NCF => CURRENT%NCF
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO
    
  END FUNCTION FIND_FILE_BYNAME
!====================================================================
!====================================================================
  FUNCTION FIND_FILE_BYNCID(LIST,NCID,FOUND) RESULT(NCF)
    IMPLICIT NONE
    TYPE(NCFILE),POINTER               :: NCF
    TYPE(NCFILELIST),    INTENT(IN) :: LIST
    INTEGER,             INTENT(IN)    :: NCID
    LOGICAL,             INTENT(OUT)   :: FOUND
    TYPE(NCFILEP), pointer             :: CURRENT, PREVIOUS
    
    NULLIFY(NCF)
    PREVIOUS => LIST%FIRST
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( NCID == CURRENT%NCF%NCID ) THEN ! FOUND THE FILE LINK
          NCF => CURRENT%NCF
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO
    
  END FUNCTION FIND_FILE_BYNCID
!====================================================================
!====================================================================
  FUNCTION COUNT_FILE_LIST(LIST) RESULT(CNT)
    IMPLICIT NONE
    INTEGER                            :: CNT
    TYPE(NCFILELIST),    INTENT(IN)    :: LIST
    TYPE(NCFILEP), pointer             :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%FIRST
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       CNT = CNT + 1
    END DO
    
  END FUNCTION COUNT_FILE_LIST
!====================================================================
!====================================================================
  SUBROUTINE PRINT_FILE_LIST(LIST)
    IMPLICIT NONE
    type(NCFILELIST), intent(IN)       :: LIST
    TYPE(NCFILEP), pointer             :: CURRENT, PREVIOUS
    INTEGER                            :: CNT
    Character(len=4)                   :: chr

    PREVIOUS => LIST%FIRST
    CURRENT  => PREVIOUS%NEXT

    IF(.NOT. ASSOCIATED(CURRENT)) THEN ! EMPTY LIST
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%%%%%% FILE LIST IS EMPTY %%%%%%%%%%%%%"
       RETURN
    ELSE
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%%%%%% PRINTING FILE LIST %%%%%%%%%%%%%"
    END IF
    
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       CNT = CNT + 1
       write(chr,'(I4.4)')CNT
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"! PRINTING FILE LIST ENTRY #"//CHR
       CALL PRINT_FILE(CURRENT%NCF)

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
    END DO
    if(DBG_SET(DBG_LOG)) &
         & write(IPT,*)"%%%%%%%%%%% END OF FILE LIST %%%%%%%%%%%%%"
  END SUBROUTINE PRINT_FILE_LIST
!====================================================================
!====================================================================
  SUBROUTINE PRINT_FILE(NCF)
    implicit none
    type(NCFILE), pointer :: NCF
    type(NCFTIME), pointer :: TM
    
    if(DBG_SET(dbg_log)) then
       WRITE(IPT,*) "======== PRINT NCFILE TYPE ======="
       if(.not. associated(NCF)) then
          WRITE(IPT,*) "THIS NCFILE HAS NOT BEEN ASSOCIATED"
          WRITE(IPT,*) "======= PRINTED NCFILE TYPE ======"
          return
       end if
       WRITE(IPT,*) "=FILENAME    ::"//trim(NCF%FNAME)
       WRITE(IPT,*) "=NCID        ::",NCF%NCID
       WRITE(IPT,*) "=WRITABLE    ::",NCF%WRITABLE
       WRITE(IPT,*) "=OPEN        ::",NCF%OPEN
       WRITE(IPT,*) "=CONNECTED   ::",NCF%CONNECTED
       WRITE(IPT,*) "=INDEFMODE   ::",NCF%INDEFMODE
       WRITE(IPT,*) "=UNLIMDIMID  ::",NCF%UNLIMDIMID
       WRITE(IPT,*) "= "
       TM => NCF%FTIME
       CALL PRINT_FTIME(TM)
       WRITE(IPT,*) "= "
       WRITE(IPT,*) "=     FILE OBJECT COUNTS"
       WRITE(IPT,*) "=nDIMS       ::",count_dim_list(NCF)
       WRITE(IPT,*) "=nATTS       ::",count_att_list(NCF)
       WRITE(IPT,*) "=nVARS       ::",count_var_list(NCF)
       IF(ASSOCIATED(NCF%INTERP_N))WRITE(IPT,*) "= HAS INTERP COEF'S TO NODES"
       IF(ASSOCIATED(NCF%INTERP_C))WRITE(IPT,*) "= HAS INTERP COEF'S TO CELLS"
       WRITE(IPT,*) "====== PRINTED NCFILE TYPE ======"
       
    end if
END SUBROUTINE PRINT_FILE
!===================================================================
!===================================================================
SUBROUTINE PRINT_FTIME(FTIME)
  USE CONTROL, only : USE_REAL_WORLD_TIME
  IMPLICIT NONE
  TYPE(NCFTIME), POINTER :: FTIME
  

  WRITE(IPT,*) "===== FILE IO TIME INFO ===="

  IF (ASSOCIATED(FTIME)) THEN
     select case(FTIME%TMtype)
     CASE(TMtype_UNKNOWN)
        WRITE(IPT,*) "=    TMTYPE  :: UNKNOWN"
     CASE(TMtype_CHAR_DATE)
        WRITE(IPT,*) "=    TMTYPE  :: CHARACTER STRING DATE"
     CASE(TMtype_INT2_MJD)
        WRITE(IPT,*) "=    TMTYPE  :: TWO INTEGER MJD"
     CASE(TMtype_FLOAT_SECONDS)
        WRITE(IPT,*) "=    TMTYPE  :: FLOATING POINT SECONDS"
     CASE(TMtype_FLOAT_DAYS)
        WRITE(IPT,*) "=    TMTYPE  :: FLOATING POINT DAYS"
     END select
     CALL PRINT_VAR(FTIME%TM1)
     CALL PRINT_VAR(FTIME%TM2)
          
     IF(USE_REAL_WORLD_TIME) THEN
        CALL PRINT_REAL_TIME(FTIME%PREV_IO,IPT,"PREV IO")
        CALL PRINT_REAL_TIME(FTIME%NEXT_IO,IPT,"NEXT IO")
     ELSE
        CALL PRINT_TIME(FTIME%PREV_IO,IPT,"PREV IO")
        CALL PRINT_TIME(FTIME%NEXT_IO,IPT,"NEXT IO")
     END IF
     CALL PRINT_TIME(FTIME%INTERVAL,IPT,"IO INTERVAL")
     
     WRITE(IPT,*) "=STK_LEN     ::",FTIME%STK_LEN

     WRITE(IPT,*) "=PREV_STKCNT ::",FTIME%PREV_STKCNT
     WRITE(IPT,*) "=NEXT_STKCNT ::",FTIME%NEXT_STKCNT
     WRITE(IPT,*) "=MAX_STKCNT  ::",FTIME%MAX_STKCNT
     WRITE(IPT,*) "=PREV_WGHT   ::",FTIME%PREV_WGHT
     WRITE(IPT,*) "=NEXT_WGHT   ::",FTIME%NEXT_WGHT
     WRITE(IPT,*) "=TIMEZONE    ::"//TRIM(FTIME%TIMEZONE)
     
  ELSE
     WRITE(IPT,*) "=  FTIME NOT ALLOCATED: THE FILE HAS NO RECOGNIZED TIME VARIABLE"
  END IF
  WRITE(IPT,*) "=END  FILE IO TIME TYPE ===="
  
END SUBROUTINE PRINT_FTIME
!====================================================================
!====================================================================
!=========================== VARIABLES ==============================
!====================================================================
!====================================================================
  FUNCTION NEW_VAR() RESULT(VAR)
    IMPLICIT NONE
!    CHARACTER(Len=*), INTENT(IN) :: NAME
    TYPE(NCVAR), POINTER :: VAR
    integer status
    
    ALLOCATE(VAR,stat=status)
    if(status/=0) CALL FATAL_ERROR("NEW_VAR: COULD NOT ALLOCATE!")

    NULLIFY(VAR%NCID)
    VAR%CONNECTED = .FALSE.
    VAR%VARID     = -1
!    VAR%VARNAME   = NAME
    VAR%VARNAME   = ' '
    VAR%XTYPE     = -1
    nullify(VAR%DIMS)
    nullify(VAR%ATTS)

    VAR%DIMS => NEW_DIMP()
    VAR%ATTS => NEW_ATTP()

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

  END FUNCTION NEW_VAR
!====================================================================
!====================================================================
  FUNCTION NEW_VARP() RESULT(VARP)
    IMPLICIT NONE
    TYPE(NCVARP), POINTER :: VARP
    integer status
    
    ALLOCATE(VARP,stat=status)
    if(status/=0) CALL FATAL_ERROR("NEW_VARP: COULD NOT ALLOCATE!")
    NULLIFY(VARP%NEXT)
    NULLIFY(VARP%VAR)

  END FUNCTION NEW_VARP
!====================================================================
!====================================================================
  FUNCTION COPY_VAR(VARIN) RESULT(VAROUT)
    IMPLICIT NONE
    TYPE(NCVAR), POINTER :: VARIN, VAROUT
    integer status

    VAROUT => NEW_VAR()

    ! DO NOT TAKE VALUES THAT BELONG TO THIS FILE
!    VAROUT%NCID      => VARIN%NCID
!    VAROUT%CONNECTED = VARIN%CONNECTED
!    VAROUT%VARID     = VARIN%VARID
    VAROUT%VARNAME   = VARIN%VARNAME
    VAROUT%XTYPE     = VARIN%XTYPE
    
    IF(Associated(VARIN%SCL_INT)) VAROUT%SCL_INT=>VARIN%SCL_INT
    IF(Associated(VARIN%VEC_INT)) VAROUT%VEC_INT=>VARIN%VEC_INT
    IF(Associated(VARIN%ARR_INT)) VAROUT%ARR_INT=>VARIN%ARR_INT
    IF(Associated(VARIN%CUB_INT)) VAROUT%CUB_INT=>VARIN%CUB_INT
    IF(Associated(VARIN%FDA_INT)) VAROUT%FDA_INT=>VARIN%FDA_INT

    IF(Associated(VARIN%SCL_FLT)) VAROUT%SCL_FLT=>VARIN%SCL_FLT
    IF(Associated(VARIN%VEC_FLT)) VAROUT%VEC_FLT=>VARIN%VEC_FLT
    IF(Associated(VARIN%ARR_FLT)) VAROUT%ARR_FLT=>VARIN%ARR_FLT
    IF(Associated(VARIN%CUB_FLT)) VAROUT%CUB_FLT=>VARIN%CUB_FLT
    IF(Associated(VARIN%FDA_FLT)) VAROUT%FDA_FLT=>VARIN%FDA_FLT

    IF(Associated(VARIN%SCL_DBL)) VAROUT%SCL_DBL=>VARIN%SCL_DBL
    IF(Associated(VARIN%VEC_DBL)) VAROUT%VEC_DBL=>VARIN%VEC_DBL
    IF(Associated(VARIN%ARR_DBL)) VAROUT%ARR_DBL=>VARIN%ARR_DBL
    IF(Associated(VARIN%CUB_DBL)) VAROUT%CUB_DBL=>VARIN%CUB_DBL
    IF(Associated(VARIN%FDA_DBL)) VAROUT%FDA_DBL=>VARIN%FDA_DBL

    IF(Associated(VARIN%SCL_CHR)) VAROUT%SCL_CHR=>VARIN%SCL_CHR
    IF(Associated(VARIN%VEC_CHR)) VAROUT%VEC_CHR=>VARIN%VEC_CHR



    IF(.not. Associated(VARIN%DIMS)) Call Fatal_Error &
         ("COPY_VAR: INPUT VARIABLE DIMS LIST IS NOT ASSOCIATED!")

    IF(.not. Associated(VARIN%ATTS)) Call Fatal_Error &
         ("COPY_VAR: INPUT VARIABLE ATTS LIST IS NOT ASSOCIATED!")

    ! USING OVERLOADED ASSIGNMENT OPERATOR
!    call copy_dim_list(VAROUT%DIMS, VARIN%DIMS)
!    call copy_att_list(VAROUT%ATTS, VARIN%ATTS)

    VAROUT%DIMS = VARIN%DIMS
    VAROUT%ATTS = VARIN%ATTS

  END FUNCTION COPY_VAR
!====================================================================
!====================================================================
  FUNCTION REFERENCE_VAR(VARIN) RESULT(VAROUT)
    IMPLICIT NONE
    TYPE(NCVAR), POINTER :: VARIN, VAROUT
    integer status

    VAROUT => NEW_VAR()

    VAROUT%NCID      => VARIN%NCID
    VAROUT%CONNECTED = VARIN%CONNECTED
    VAROUT%VARID     = VARIN%VARID
    VAROUT%VARNAME   = VARIN%VARNAME
    VAROUT%XTYPE     = VARIN%XTYPE
    VAROUT%DIMS%NEXT => VARIN%DIMS%NEXT
    VAROUT%ATTS%NEXT => VARIN%ATTS%NEXT

  END FUNCTION REFERENCE_VAR
!====================================================================
!====================================================================
  SUBROUTINE DELETE_VARP_BYNAME(LIST,NAME,FOUND)
    IMPLICIT NONE
    ! NCFILE IS ALWAYS THE HEAD FOR A VAR LINK LIST
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    CHARACTER(LEN=*),INTENT(IN)                 :: NAME
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCVARP), pointer                       :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( NAME == CURRENT%VAR%VARNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_VAR(CURRENT%VAR)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_VARP_BYNAME
!====================================================================
!====================================================================
  SUBROUTINE DELETE_VARP_BYVARID(LIST,VARID,FOUND)
    IMPLICIT NONE
    ! NCFILE IS ALWAYS THE HEAD FOR A VAR LINK LIST
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    INTEGER,                       INTENT(IN)   :: VARID
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCVARP), pointer                       :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( VARID == CURRENT%VAR%VARID ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_VAR(CURRENT%VAR)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_VARP_BYVARID
!====================================================================
!====================================================================
  SUBROUTINE DELETE_VAR_LIST(LIST)
    IMPLICIT NONE
    ! NCFILE IS ALWAYS THE HEAD FOR A VAR LINK LIST
    TYPE(NCFILE), POINTER                       :: LIST
    TYPE(NCVARP), pointer                       :: CURRENT, PREVIOUS

    IF(.not.ASSOCIATED(LIST)) RETURN
    
    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST

       PREVIOUS%NEXT => CURRENT%NEXT
       CALL KILL_VAR(CURRENT%VAR)
       DEALLOCATE(CURRENT)
       CURRENT => PREVIOUS%NEXT

    END DO

   
  END SUBROUTINE DELETE_VAR_LIST
!====================================================================
!====================================================================
 SUBROUTINE KILL_VAR(VAR)
   TYPE(NCVAR),POINTER :: VAR
   INTEGER STATUS

   IF(.not.ASSOCIATED(VAR)) RETURN
   
   CALL DELETE_DIM_LIST(VAR)
   DEALLOCATE(VAR%DIMS)
   
   CALL DELETE_ATT_LIST(VAR)
   DEALLOCATE(VAR%ATTS)

   DEALLOCATE(VAR,STAT=STATUS)
   IF(STATUS /= 0) CALL FATAL_ERROR("KILL_VAR: COULD NOT DEALLOCATE")
   NULLIFY(VAR)
   
 END SUBROUTINE KILL_VAR
!====================================================================
!====================================================================
!!$ SUBROUTINE INSERT_VARP_BYNAME(LIST,NAME,FOUND)
!!$    ! IF FOUND DO NOT INSERT DUPLICATE, RETURN FOUND
!!$    ! ALWAYS INSERT NCVARP AT THE END OF THE LIST
!!$    IMPLICIT NONE
!!$    CHARACTER(LEN=NF90_MAX_NAME+1),INTENT(IN)   :: NAME
!!$    LOGICAL,                       INTENT(OUT)  :: FOUND
!!$    TYPE(NCFILE),                  INTENT(INOUT):: LIST
!!$    TYPE(NCVARP), pointer                       :: CURRENT, PREVIOUS
!!$    
!!$    PREVIOUS => LIST%VARS
!!$    CURRENT  => PREVIOUS%NEXT
!!$    FOUND = .FALSE.
!!$    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
!!$    DO
!!$       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
!!$       IF( NAME == CURRENT%VAR%VARNAME ) THEN ! FOUND THE FILE LINK
!!$          FOUND = .TRUE.
!!$          RETURN
!!$       ELSE
!!$          PREVIOUS => PREVIOUS%NEXT
!!$          CURRENT  => CURRENT%NEXT
!!$       END IF
!!$    END DO
!!$    ! NOT FOUND - ADD NEW FILE TO END OF LIST
!!$    CURRENT => NEW_VARP()
!!$    CURRENT%VAR => NEW_VAR(NAME)
!!$    
!!$  END SUBROUTINE INSERT_VARP_BYNAME
!====================================================================
!====================================================================
 SUBROUTINE INSERT_VARP_BYVAR(LIST,VAR,FOUND)
    ! IF FOUND DO NOT INSERT DUPLICATE, RETURN FOUND
    ! ALWAYS INSERT NCVARP AT THE END OF THE LIST
    IMPLICIT NONE
    TYPE(NCVAR),POINTER                         :: VAR
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    TYPE(NCVARP),POINTER                        :: CURRENT, PREVIOUS
    INTEGER CNT

    IF(.NOT.ASSOCIATED(VAR))&
         & CALL FATAL_ERROR("INSERT_VARP_BYVAR: VAR NOT ASSOCIATED!")

    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.
    CNT = 1
    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       IF( VAR%VARNAME == CURRENT%VAR%VARNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          RETURN
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
          CNT = CNT + 1
       END IF
    END DO
    ! NOT FOUND - ADD NEW FILE TO END OF LIST
    PREVIOUS%NEXT => NEW_VARP()
    PREVIOUS%NEXT%VAR => VAR
    VAR%VARID = CNT

  END SUBROUTINE INSERT_VARP_BYVAR
!====================================================================
!====================================================================
  FUNCTION FIND_VAR_BYNAME(LIST,NAME,FOUND) RESULT(VAR)
    IMPLICIT NONE
    TYPE(NCVAR), POINTER                :: VAR
    TYPE(NCFILE),         INTENT(IN) :: LIST
    CHARACTER(LEN=*),INTENT(IN)    :: NAME
    LOGICAL,              INTENT(OUT)   :: FOUND
    TYPE(NCVARP)                ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(VAR)
    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( NAME == CURRENT%VAR%VARNAME ) THEN ! FOUND THE FILE LINK
          VAR => CURRENT%VAR
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_VAR_BYNAME
!====================================================================
!====================================================================
  FUNCTION FIND_VAR_BYVARID(LIST,VARID,FOUND) RESULT(VAR)
    IMPLICIT NONE
    TYPE(NCVAR), POINTER                :: VAR
    TYPE(NCFILE),         INTENT(IN) :: LIST
    INTEGER,              INTENT(IN)    :: VARID
    LOGICAL,              INTENT(OUT)   :: FOUND
    TYPE(NCVARP)                ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(VAR)
    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( VARID == CURRENT%VAR%VARID ) THEN ! FOUND THE FILE LINK
          VAR => CURRENT%VAR
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO
    
  END FUNCTION FIND_VAR_BYVARID    
!====================================================================
!====================================================================
  FUNCTION COUNT_VAR_LIST(LIST) RESULT(CNT)
    IMPLICIT NONE
    INTEGER                            :: CNT
    TYPE(NCFILE),        INTENT(IN)    :: LIST
    TYPE(NCVARP)              ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       CNT = CNT + 1
    END DO
    
  END FUNCTION COUNT_VAR_LIST
!====================================================================
!====================================================================
  FUNCTION COUNT_UNLIMITED_VARS(LIST) RESULT(CNT)
    IMPLICIT NONE
    INTEGER                            :: CNT
    TYPE(NCFILE),        INTENT(IN)    :: LIST
    TYPE(NCVARP)              ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF(HAS_UNLIMITED(CURRENT%VAR)) CNT = CNT + 1

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
    END DO
    
  END FUNCTION COUNT_UNLIMITED_VARS
!====================================================================
!====================================================================
  SUBROUTINE PRINT_VAR_LIST(LIST)
    IMPLICIT NONE
    type(NCFILE), intent(IN)       :: LIST
    TYPE(NCVARP)              ,POINTER :: CURRENT, PREVIOUS
    INTEGER                            :: CNT
    Character(len=4)                   :: chr

    PREVIOUS => LIST%VARS
    CURRENT  => PREVIOUS%NEXT

    IF(.NOT. ASSOCIATED(CURRENT)) THEN ! EMPTY LIST
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%%%%%% VARIABLE LIST IS EMPTY %%%%%%%%%%%%%"
       RETURN
    ELSE
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%%%%%% PRINTING VARIABLE LIST %%%%%%%%%%%%%"
    END IF
    
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       CNT = CNT + 1
       write(chr,'(I4.4)')CNT
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"! PRINTING VARIABLE LIST ENTRY #"//CHR
       CALL PRINT_VAR(CURRENT%VAR)

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
    END DO
    if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%%%%%% END OF VARIABLE LIST %%%%%%%%%%%%%"
  END SUBROUTINE PRINT_VAR_LIST
!====================================================================
!====================================================================
SUBROUTINE PRINT_VAR(VAR)
  implicit none
  type(NCVAR), POINTER, intent(IN) :: VAR
  
  if(DBG_SET(dbg_log)) then
     WRITE(IPT,*) "======== PRINT NCVAR TYPE ======="

     if(.not. associated(VAR)) then
        WRITE(IPT,*) "THIS NCVAR HAS NOT BEEN ASSOCIATED"
        WRITE(IPT,*) "======= PRINTED NCVAR TYPE ======"
        return
     end if

     WRITE(IPT,*) "VARNAME     ::"//trim(VAR%VARNAME)
     WRITE(IPT,*) "VARID       ::",VAR%VARID
      if(.not. associated(VAR%NCID)) then
         WRITE(IPT,*) "NCID        :: NOT ASSOCIATED"
      else
         WRITE(IPT,*) "NCID        ::",VAR%NCID
     end if
     WRITE(IPT,*) "CONNECTED   ::",VAR%CONNECTED
     WRITE(IPT,*) "CURR_STKCNT ::",VAR%CURR_STKCNT
     
     select case(VAR%XTYPE)
     case(NF90_CHAR)
        WRITE(IPT,*) "XYTPE       :: CHAR"
     case(NF90_BYTE)
        WRITE(IPT,*) "XYTPE       :: BYTE"
     case(NF90_SHORT)
        WRITE(IPT,*) "XYTPE       :: SHORT"
     case(NF90_INT)
        WRITE(IPT,*) "XYTPE       :: INT"
     case(NF90_FLOAT)
        WRITE(IPT,*) "XYTPE       :: FLOAT"
     case(NF90_DOUBLE)
        WRITE(IPT,*) "XYTPE       :: DOUBLE"
     END select

     IF(ASSOCIATED(VAR%SCL_INT)) WRITE(IPT,*) "ASSOCIATED  :: SCL_INT"
     IF(ASSOCIATED(VAR%VEC_INT)) WRITE(IPT,*) "ASSOCIATED  :: VEC_INT"
     IF(ASSOCIATED(VAR%ARR_INT)) WRITE(IPT,*) "ASSOCIATED  :: ARR_INT"
     IF(ASSOCIATED(VAR%CUB_INT)) WRITE(IPT,*) "ASSOCIATED  :: CUB_INT"
     IF(ASSOCIATED(VAR%FDA_INT)) WRITE(IPT,*) "ASSOCIATED  :: FDA_INT"

     IF(ASSOCIATED(VAR%SCL_FLT)) WRITE(IPT,*) "ASSOCIATED  :: SCL_FLT"
     IF(ASSOCIATED(VAR%VEC_FLT)) WRITE(IPT,*) "ASSOCIATED  :: VEC_FLT"
     IF(ASSOCIATED(VAR%ARR_FLT)) WRITE(IPT,*) "ASSOCIATED  :: ARR_FLT"
     IF(ASSOCIATED(VAR%CUB_FLT)) WRITE(IPT,*) "ASSOCIATED  :: CUB_FLT"
     IF(ASSOCIATED(VAR%FDA_FLT)) WRITE(IPT,*) "ASSOCIATED  :: FDA_FLT"

     IF(ASSOCIATED(VAR%SCL_DBL)) WRITE(IPT,*) "ASSOCIATED  :: SCL_DBL"
     IF(ASSOCIATED(VAR%VEC_DBL)) WRITE(IPT,*) "ASSOCIATED  :: VEC_DBL"
     IF(ASSOCIATED(VAR%ARR_DBL)) WRITE(IPT,*) "ASSOCIATED  :: ARR_DBL"
     IF(ASSOCIATED(VAR%CUB_DBL)) WRITE(IPT,*) "ASSOCIATED  :: CUB_DBL"
     IF(ASSOCIATED(VAR%FDA_DBL)) WRITE(IPT,*) "ASSOCIATED  :: FDA_DBL"

     IF(ASSOCIATED(VAR%SCL_CHR)) WRITE(IPT,*) "ASSOCIATED  :: SCL_CHR"
     IF(ASSOCIATED(VAR%VEC_CHR)) WRITE(IPT,*) "ASSOCIATED  :: VEC_CHR"



     WRITE(IPT,*) "DIMS        ::",count_dim_list(VAR)
     WRITE(IPT,*) "ATTS        ::",count_att_list(VAR)
     WRITE(IPT,*) "======= PRINTED NCVAR TYPE ======"

  end if
END SUBROUTINE PRINT_VAR
!====================================================================
!====================================================================
!====================== ATTRIBUTES ==================================
!====================================================================
!====================================================================
  FUNCTION NEW_ATT() RESULT(ATT)
    IMPLICIT NONE
!    CHARACTER(LEN=*),INTENT(IN)   :: NAME
    TYPE(NCATT), POINTER :: ATT
    integer status
    
    ALLOCATE(ATT,stat=status)
    if(status/=0) CALL FATAL_ERROR("NEW_ATT: COULD NOT ALLOCATE!")

!    ATT%ATTNAME   = NAME
    ATT%ATTNAME   = ' '
    ATT%LEN       = -1
    ATT%ATTID     = 0
    ATT%XTYPE     = -1
  END FUNCTION NEW_ATT
!====================================================================
!====================================================================
  FUNCTION NEW_ATTP() RESULT(ATTP)
    IMPLICIT NONE
    TYPE(NCATTP), POINTER :: ATTP
    integer status
    
    ALLOCATE(ATTP,stat=status)
    if(status/=0) CALL FATAL_ERROR("NEW_ATTP COULD NOT ALLOCATE!")
    NULLIFY(ATTP%NEXT)
    NULLIFY(ATTP%ATT)

  END FUNCTION NEW_ATTP
!====================================================================
!====================================================================
  FUNCTION COPY_ATT(ATTIN) RESULT(ATTOUT)
    IMPLICIT NONE
    TYPE(NCATT), POINTER, INTENT(IN) :: ATTIN
    TYPE(NCATT), POINTER :: ATTOUT
    integer status
    
    IF(.not. Associated(ATTIN)) CALL FATAL_ERROR("THE ARGUMENT MUST BE&
         & ASSOCIAED FOR COPY_ATT")


    ATTOUT=>NEW_ATT()

    ATTOUT%ATTNAME   = ATTIN%ATTNAME
    ATTOUT%LEN       = ATTIN%LEN 
    ATTOUT%ATTID     = ATTIN%ATTID
    ATTOUT%XTYPE     = ATTIN%XTYPE

    IF (Allocated(ATTIN%int)) THEN
       ALLOCATE(ATTOUT%int(ATTIN%LEN),stat=status)
       if(status/=0) CALL FATAL_ERROR("COPY_ATT COULD NOT ALLOCATE INT!")
       ATTOUT%int = ATTIN%int
    END IF

    IF (Allocated(ATTIN%flt)) THEN
       ALLOCATE(ATTOUT%flt(ATTIN%LEN),stat=status)
       if(status/=0) CALL FATAL_ERROR("COPY_ATT COULD NOT ALLOCATE FLT!")
       ATTOUT%flt = ATTIN%flt
    END IF

    IF (Allocated(ATTIN%dbl)) THEN
       ALLOCATE(ATTOUT%dbl(ATTIN%LEN),stat=status)
       if(status/=0) CALL FATAL_ERROR("COPY_ATT COULD NOT ALLOCATE DBL!")
       ATTOUT%dbl = ATTIN%dbl
    END IF

    IF (Allocated(ATTIN%chr)) THEN
       ALLOCATE(ATTOUT%chr(size(ATTIN%chr)),stat=status)
       if(status/=0) CALL FATAL_ERROR("COPY_ATT COULD NOT ALLOCATE CHR!")
       ATTOUT%chr = ATTIN%chr
    END IF
    
  END FUNCTION COPY_ATT
!====================================================================
!====================================================================
  SUBROUTINE COPY_ATT_LIST(ATTPOUT,ATTPIN)
    IMPLICIT NONE
!    TYPE(NCATTP), POINTER, INTENT(IN) :: ATTPIN
!    TYPE(NCATTP), POINTER, INTENT(INOUT) :: ATTPOUT
    TYPE(NCATTP), TARGET, INTENT(IN) :: ATTPIN
    TYPE(NCATTP), TARGET, INTENT(OUT) :: ATTPOUT
    TYPE(NCATT), POINTER              :: ATT
    TYPE(NCATTP),POINTER              :: CURRENT_IN, PREVIOUS_IN
    TYPE(NCATTP),POINTER              :: CURRENT_OUT, PREVIOUS_OUT

    integer status
    
!!$    IF(.not. Associated(ATTPIN)) CALL FATAL_ERROR("THE INPUT ARGUMENT MUST BE&
!!$         & ASSOCIAED FOR COPY_ATT_LIST?")
!!$
!!$    IF(.not. Associated(ATTPOUT)) CALL FATAL_ERROR("THE OUTPUT ARGUMENT MUST BE&
!!$         & ASSOCIAED FOR COPY_ATT_LIST?")

    PREVIOUS_IN => ATTPIN
    CURRENT_IN => PREVIOUS_IN%NEXT

    PREVIOUS_OUT => ATTPOUT
    CURRENT_OUT => PREVIOUS_OUT%NEXT

    ! I DON'T THINK THAT CURRENT_OUT SERVES ANY PURPOSE HERE?
    ! I DON'T THINK THAT PREVIOUS_IN SERVES ANY PURPOSE HERE?
    
    
    DO
       IF(.NOT. ASSOCIATED(CURRENT_IN)) THEN
          RETURN ! END OF ATT IN LIST
       END IF

       IF(.NOT. ASSOCIATED(CURRENT_IN%ATT)) THEN
          CALL FATAL_ERROR("COPY_ATT_LIST: FOUND NULL DIM POINTER IN THE LIST")
       END IF
       
       PREVIOUS_OUT%NEXT => NEW_ATTP()
       PREVIOUS_OUT%NEXT%ATT => COPY_ATT(CURRENT_IN%ATT)
       PREVIOUS_OUT%NEXT%NEXT => CURRENT_OUT
       
       ! INCRIMENT THROUGH THE LIST
       PREVIOUS_OUT => PREVIOUS_OUT%NEXT
       ! DO NOT INCRIMENT  CURRENT_OUT
       
       PREVIOUS_IN => PREVIOUS_IN%NEXT
       CURRENT_IN  => CURRENT_IN%NEXT
          
    END DO
    


  END SUBROUTINE COPY_ATT_LIST
!====================================================================
!====================================================================
  SUBROUTINE DELETE_NCF_ATTP_BYNAME(LIST,NAME,FOUND)
    IMPLICIT NONE
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    CHARACTER(LEN=*),INTENT(IN)   :: NAME
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCATTP)                        ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( NAME == CURRENT%ATT%ATTNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_ATT(CURRENT%ATT)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_NCF_ATTP_BYNAME
!====================================================================
!====================================================================
  SUBROUTINE DELETE_NCF_ATTP_BYATTID(LIST,ATTID,FOUND)
    IMPLICIT NONE
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    INTEGER,                       INTENT(IN)   :: ATTID
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCATTP)                        ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( ATTID == CURRENT%ATT%ATTID ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_ATT(CURRENT%ATT)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_NCF_ATTP_BYATTID
!====================================================================
!====================================================================
  SUBROUTINE DELETE_NCF_ATTP_LIST(LIST)
    IMPLICIT NONE
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    TYPE(NCATTP)                        ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST

       ! IF FOUND DELETE IT
       PREVIOUS%NEXT => CURRENT%NEXT
       CALL KILL_ATT(CURRENT%ATT)
       DEALLOCATE(CURRENT)
       CURRENT => PREVIOUS%NEXT

    END DO
   
  END SUBROUTINE DELETE_NCF_ATTP_LIST
!====================================================================
!====================================================================
  SUBROUTINE DELETE_VAR_ATTP_BYNAME(LIST,NAME,FOUND)
    IMPLICIT NONE
    TYPE(NCVAR),                   INTENT(INOUT):: LIST
!    CHARACTER(LEN=NF90_MAX_NAME+1),INTENT(IN)   :: NAME
    CHARACTER(LEN=*),INTENT(IN)   :: NAME
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCATTP)                        ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( NAME == CURRENT%ATT%ATTNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_ATT(CURRENT%ATT)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_VAR_ATTP_BYNAME
!====================================================================
!====================================================================
  SUBROUTINE DELETE_VAR_ATTP_BYATTID(LIST,ATTID,FOUND)
    IMPLICIT NONE
    TYPE(NCVAR),                   INTENT(INOUT):: LIST
    INTEGER,                       INTENT(IN)   :: ATTID
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCATTP), POINTER                       :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( ATTID == CURRENT%ATT%ATTID ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_ATT(CURRENT%ATT)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_VAR_ATTP_BYATTID
!====================================================================
!====================================================================
  SUBROUTINE DELETE_VAR_ATTP_LIST(LIST)
    IMPLICIT NONE
    TYPE(NCVAR), INTENT(INOUT):: LIST
    TYPE(NCATTP),POINTER      :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       PREVIOUS%NEXT => CURRENT%NEXT
       CALL KILL_ATT(CURRENT%ATT)
       DEALLOCATE(CURRENT)
       CURRENT => PREVIOUS%NEXT

    END DO

  END SUBROUTINE DELETE_VAR_ATTP_LIST
!====================================================================
!====================================================================
 SUBROUTINE KILL_ATT(ATT)
   TYPE(NCATT),POINTER :: ATT
   INTEGER STATUS
   IF(ASSOCIATED(ATT)) THEN
      IF(ALLOCATED(ATT%INT)) DEALLOCATE(ATT%INT)
      IF(ALLOCATED(ATT%FLT)) DEALLOCATE(ATT%FLT)
      IF(ALLOCATED(ATT%DBL)) DEALLOCATE(ATT%DBL)
      IF(ALLOCATED(ATT%CHR)) DEALLOCATE(ATT%CHR)

      DEALLOCATE(ATT,STAT=STATUS)
      IF(STATUS /= 0) CALL FATAL_ERROR("KILL_ATT: COULD NOT DEALLOCATE")
      NULLIFY(ATT)
   END IF

 END SUBROUTINE KILL_ATT
!====================================================================
!====================================================================
!!$  SUBROUTINE INSERT_NCF_ATTP_BYNAME(LIST,NAME,FOUND)
!!$    ! IF FOUND DO NOT INSERT DUPLICATE, RETURN FOUND
!!$    ! ALWAYS INSERT NCVARP AT THE END OF THE LIST
!!$    IMPLICIT NONE
!!$    CHARACTER(LEN=NF90_MAX_NAME+1),INTENT(IN)   :: NAME
!!$    LOGICAL,                       INTENT(OUT)  :: FOUND
!!$    TYPE(NCFILE),                  INTENT(INOUT):: LIST
!!$    TYPE(NCATTP)          ,POINTER              :: CURRENT, PREVIOUS
!!$    INTEGER CNT
!!$    PREVIOUS => LIST%ATTS
!!$    CURRENT  => PREVIOUS%NEXT
!!$    FOUND = .FALSE.
!!$    CNT = 1
!!$
!!$    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
!!$    DO
!!$       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
!!$       IF( NAME == CURRENT%ATT%ATTNAME ) THEN ! FOUND THE FILE LINK
!!$          FOUND = .TRUE.
!!$          RETURN
!!$       ELSE
!!$          PREVIOUS => PREVIOUS%NEXT
!!$          CURRENT  => CURRENT%NEXT
!!$          CNT = CNT + 1
!!$       END IF
!!$    END DO
!!$    ! NOT FOUND - ADD NEW FILE TO END OF LIST
!!$    CURRENT => NEW_ATTP()
!!$    CURRENT%ATT => NEW_ATT(NAME)
!!$    CURRENT%ATT%ATTID = CNT ! SET THE ATTID 
!!$
!!$  END SUBROUTINE INSERT_NCF_ATTP_BYNAME
!====================================================================
!====================================================================
  SUBROUTINE INSERT_NCF_ATTP_BYATT(LIST,ATT,FOUND)
    ! IF FOUND DO NOT INSERT DUPLICATE, RETURN FOUND
    ! ALWAYS INSERT NCVARP AT THE END OF THE LIST
    IMPLICIT NONE
    TYPE(NCATT), POINTER                        :: ATT
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    TYPE(NCATTP),POINTER                        :: CURRENT, PREVIOUS
    INTEGER CNT
    IF(.NOT.ASSOCIATED(ATT))&
         & CALL FATAL_ERROR("INSERT_NCF_ATTP_BYATT: ATT NOT ASSOCIATED!")

    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.
    CNT= 1

    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       IF( ATT%ATTNAME == CURRENT%ATT%ATTNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          RETURN
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
          CNT = CNT + 1
       END IF
    END DO
    ! NOT FOUND - ADD NEW FILE TO END OF LIST
     
    PREVIOUS%NEXT => NEW_ATTP()
    PREVIOUS%NEXT%ATT => ATT
    ATT%ATTID = CNT ! SET THE ATTID 

    PREVIOUS%NEXT%NEXT => CURRENT
    

  END SUBROUTINE INSERT_NCF_ATTP_BYATT
!====================================================================
!====================================================================
!!$  SUBROUTINE INSERT_VAR_ATTP_BYNAME(LIST,NAME,FOUND)
!!$    ! IF FOUND DO NOT INSERT DUPLICATE, RETURN FOUND
!!$    ! ALWAYS INSERT NCVARP AT THE END OF THE LIST
!!$    IMPLICIT NONE
!!$    CHARACTER(LEN=NF90_MAX_NAME+1),INTENT(IN)   :: NAME
!!$    LOGICAL,                       INTENT(OUT)  :: FOUND
!!$    TYPE(NCVAR),                   INTENT(INOUT):: LIST
!!$    TYPE(NCATTP)                        ,POINTER :: CURRENT, PREVIOUS
!!$    INTEGER CNT
!!$    PREVIOUS => LIST%ATTS
!!$    CURRENT  => PREVIOUS%NEXT
!!$    FOUND = .FALSE.
!!$    CNT = 1
!!$    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
!!$    DO
!!$       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
!!$       IF( NAME == CURRENT%ATT%ATTNAME ) THEN ! FOUND THE FILE LINK
!!$          FOUND = .TRUE.
!!$          RETURN
!!$       ELSE
!!$          PREVIOUS => PREVIOUS%NEXT
!!$          CURRENT  => CURRENT%NEXT
!!$          CNT = CNT + 1
!!$       END IF
!!$    END DO
!!$    ! NOT FOUND - ADD NEW FILE TO END OF LIST
!!$    CURRENT => NEW_ATTP()
!!$    CURRENT%ATT => NEW_ATT(NAME)
!!$    CURRENT%ATT%ATTID = CNT ! SET THE ATTID 
!!$  END SUBROUTINE INSERT_VAR_ATTP_BYNAME
!====================================================================
!====================================================================
  SUBROUTINE INSERT_VAR_ATTP_BYATT(LIST,ATT,FOUND)
    ! IF FOUND DO NOT INSERT DUPLICATE, RETURN FOUND
    ! ALWAYS INSERT NCVARP AT THE END OF THE LIST
    IMPLICIT NONE
    TYPE(NCATT),POINTER                         :: ATT
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCVAR),                   INTENT(INOUT):: LIST
    TYPE(NCATTP),POINTER                        :: CURRENT, PREVIOUS
    INTEGER CNT
    IF(.NOT.ASSOCIATED(ATT))&
         & CALL FATAL_ERROR("INSERT_VAR_ATTP_BYATT: ATT NOT ASSOCIATED!")

    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.
    CNT = 1
    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       IF( ATT%ATTNAME == CURRENT%ATT%ATTNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          RETURN
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
          CNT = CNT + 1
       END IF
    END DO
    ! NOT FOUND - ADD NEW FILE TO END OF LIST
    PREVIOUS%NEXT => NEW_ATTP()
    PREVIOUS%NEXT%ATT => ATT

    ATT%ATTID = CNT ! SET THE ATTID 

  END SUBROUTINE INSERT_VAR_ATTP_BYATT
!====================================================================
!====================================================================
  FUNCTION FIND_NCF_ATT_BYNAME(LIST,NAME,FOUND) RESULT(ATT)
    IMPLICIT NONE
    TYPE(NCATT), POINTER               :: ATT
    TYPE(NCFILE),         INTENT(IN):: LIST
    CHARACTER(LEN=*),INTENT(IN)    :: NAME
    LOGICAL,              INTENT(OUT)  :: FOUND
    TYPE(NCATTP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(ATT)
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( TRIM(NAME) == TRIM(CURRENT%ATT%ATTNAME) ) THEN ! FOUND THE FILE LINK
          ATT => CURRENT%ATT
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_NCF_ATT_BYNAME
!====================================================================
!====================================================================
  FUNCTION FIND_NCF_ATT_BYATTID(LIST,ATTID,FOUND) RESULT(ATT)
    IMPLICIT NONE
    TYPE(NCATT), POINTER               :: ATT
    INTEGER,              INTENT(IN)   :: ATTID
    TYPE(NCFILE),         INTENT(IN):: LIST
    LOGICAL,              INTENT(OUT)  :: FOUND
    TYPE(NCATTP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(ATT)
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( ATTID == CURRENT%ATT%ATTID ) THEN ! FOUND THE FILE LINK
          ATT => CURRENT%ATT
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_NCF_ATT_BYATTID
!====================================================================
!====================================================================
  FUNCTION FIND_VAR_ATT_BYNAME(LIST,NAME,FOUND) RESULT(ATT)
    IMPLICIT NONE
    TYPE(NCATT), POINTER               :: ATT
    TYPE(NCVAR),         INTENT(IN) :: LIST
    CHARACTER(LEN=*),INTENT(IN)    :: NAME
    LOGICAL,             INTENT(OUT)   :: FOUND
    TYPE(NCATTP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(ATT)
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( TRIM(NAME) == TRIM(CURRENT%ATT%ATTNAME) ) THEN ! FOUND THE FILE LINK
          ATT => CURRENT%ATT
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_VAR_ATT_BYNAME
!====================================================================
!====================================================================
  FUNCTION FIND_VAR_ATT_BYATTID(LIST,ATTID,FOUND) RESULT(ATT)
    IMPLICIT NONE
    TYPE(NCATT), POINTER               :: ATT
    INTEGER,              INTENT(IN)   :: ATTID
    TYPE(NCVAR),          INTENT(IN):: LIST
    LOGICAL,              INTENT(OUT)  :: FOUND
    TYPE(NCATTP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(ATT)
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( ATTID == CURRENT%ATT%ATTID ) THEN ! FOUND THE FILE LINK
          ATT => CURRENT%ATT
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_VAR_ATT_BYATTID
!====================================================================
!====================================================================
  FUNCTION COUNT_NCF_ATT_LIST(LIST) RESULT(CNT)
    IMPLICIT NONE
    INTEGER                            :: CNT
    TYPE(NCFILE),    INTENT(IN)        :: LIST
    TYPE(NCATTP)               ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       CNT = CNT + 1
    END DO
    
  END FUNCTION COUNT_NCF_ATT_LIST
!====================================================================
!====================================================================
  FUNCTION COUNT_VAR_ATT_LIST(LIST) RESULT(CNT)
    IMPLICIT NONE
    INTEGER                            :: CNT
    TYPE(NCVAR),    INTENT(IN)         :: LIST
    TYPE(NCATTP)               ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       CNT = CNT + 1
    END DO
    
  END FUNCTION COUNT_VAR_ATT_LIST
!====================================================================
!====================================================================
  SUBROUTINE PRINT_NCF_ATT_LIST(LIST)
    IMPLICIT NONE
    type(NCFILE), intent(IN)           :: LIST
    TYPE(NCATTP)               ,POINTER :: CURRENT, PREVIOUS
    INTEGER                            :: CNT
    Character(len=4)                   :: chr

    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT

    IF(.NOT. ASSOCIATED(CURRENT)) THEN ! EMPTY LIST
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%%%%%% FILE ATTRIBUTE LIST IS EMPTY %%%%%%%%%%%%%"
       RETURN
    ELSE
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%% PRINTING GLOBAL ATTRIBUTE LIST %%%%%%%%%"
    END IF
    
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       CNT = CNT + 1
       write(chr,'(I4.4)')CNT
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"! PRINTING ATTRIBUTE LIST ENTRY #"//CHR
       CALL PRINT_ATT(CURRENT%ATT)

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
    END DO
    if(DBG_SET(DBG_LOG)) &
         & write(IPT,*)"%%%%%%%%%%% END OF ATTRIBUTE LIST %%%%%%%%%%%%%"
  END SUBROUTINE PRINT_NCF_ATT_LIST
!====================================================================
!====================================================================
  SUBROUTINE PRINT_VAR_ATT_LIST(LIST)
    IMPLICIT NONE
    type(NCVAR), intent(IN)            :: LIST
    TYPE(NCATTP)               ,POINTER :: CURRENT, PREVIOUS
    INTEGER                            :: CNT
    Character(len=4)                   :: chr

    PREVIOUS => LIST%ATTS
    CURRENT  => PREVIOUS%NEXT

    IF(.NOT. ASSOCIATED(CURRENT)) THEN ! EMPTY LIST
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%%% VAIABLE ATTRIBUTE LIST IS EMPTY %%%%%%%%%%"
       RETURN
    ELSE
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%% PRINTING VARIALBE: "//TRIM(LIST%VARNAME)//"&
            &; ATTRIBUTE LIST %%%%%%%%"
    END IF
    
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       CNT = CNT + 1
       write(chr,'(I4.4)')CNT
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"! PRINTING ATTRIBUTE LIST ENTRY #"//CHR
       CALL PRINT_ATT(CURRENT%ATT)

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
    END DO
    if(DBG_SET(DBG_LOG)) &
         & write(IPT,*)"%%%%%%%%%%% END OF ATTRIBUTE LIST %%%%%%%%%%%%%"
  END SUBROUTINE PRINT_VAR_ATT_LIST
!====================================================================
!====================================================================
SUBROUTINE PRINT_ATT(ATT)
  implicit none
  type(NCATT), pointer, intent(IN) :: ATT
  integer I

  if(DBG_SET(dbg_log)) then
     WRITE(IPT,*) "======== PRINT NCATT TYPE ======="
     if(.not. associated(ATT)) then
        WRITE(IPT,*) "THIS NCATT HAS NOT BEEN ASSOCIATED"
        WRITE(IPT,*) "======= PRINTED NCATT TYPE ======"
        return
     end if
     WRITE(IPT,*) "ATTNAME::"//TRIM(ATT%ATTNAME)
     WRITE(IPT,*) "LEN    ::",ATT%LEN
     WRITE(IPT,*) "ATTID  ::",ATT%ATTID
     select case(ATT%XTYPE)
     case(NF90_CHAR)
        WRITE(IPT,*) "XYTPE  ::CHAR"
        IF (.not. Allocated(ATT%chr)) then
           WRITE(IPT,*) "CHAR   :: Not allocated!"
        ELSE
           DO I = 1,size(ATT%chr)
              WRITE(IPT,*) "CHAR   ::"//trim(ATT%chr(i))
           END DO
        END IF
     case(NF90_BYTE)
        WRITE(IPT,*) "XYTPE  ::BYTE - TYPE NOT DEFINED"
     case(NF90_SHORT)
        WRITE(IPT,*) "XYTPE  ::SHORT - TYPE NOT DEFINED"
     case(NF90_INT)
        WRITE(IPT,*) "XYTPE  ::INT"
        IF (.not. Allocated(ATT%int)) then
           WRITE(IPT,*) "INT    :: Not allocated!"
        ELSE
           write(IPT,'(I8)') ATT%int
        END IF
     case(NF90_FLOAT)
        WRITE(IPT,*) "XYTPE  ::FLOAT"
        IF (.not. Allocated(ATT%flt)) then
           WRITE(IPT,*) "FLOAT  :: Not allocated!"
        ELSE
           write(IPT,'(ES14.3)') ATT%flt
        END IF
     case(NF90_DOUBLE)
        WRITE(IPT,*) "XYTPE  ::DOUBLE"
        IF (.not. Allocated(ATT%DBL)) then
           WRITE(IPT,*) "DOUBLE :: Not allocated!"
        ELSE
           write(IPT,'(ES14.3)') ATT%dbl
        END IF
     END select
     WRITE(IPT,*) "======= PRINTED NCATT TYPE ======"
  end if
END SUBROUTINE PRINT_ATT
!====================================================================
!====================================================================
!====================== DIMENSIONS ==================================
!====================================================================
!====================================================================
  FUNCTION NEW_DIM() RESULT(DIM)
    IMPLICIT NONE
    TYPE(NCDIM), POINTER :: DIM
    integer status
    
    ALLOCATE(DIM,stat=status)
    if(status/=0) CALL FATAL_ERROR("ALLOC_DIM: COULD NOT ALLOCATE!")

    DIM%DIMID     = -1
    DIM%DIMNAME   =" "
    DIM%DIM       = -1
    DIM%UNLIMITED = .false.
  END FUNCTION NEW_DIM
!====================================================================
!====================================================================
  SUBROUTINE COPY_DIM_LIST(DIMPOUT,DIMPIN)
    IMPLICIT NONE
!    TYPE(NCDIMP),POINTER,INTENT(INOUT):: DIMPOUT
!    TYPE(NCDIMP), POINTER, INTENT(IN) :: DIMPIN
    TYPE(NCDIMP),TARGET, INTENT(OUT):: DIMPOUT
    TYPE(NCDIMP), TARGET, INTENT(IN) :: DIMPIN
    TYPE(NCDIM), POINTER              :: DIM
    TYPE(NCDIMP),POINTER              :: CURRENT_IN, PREVIOUS_IN
    TYPE(NCDIMP),POINTER              :: CURRENT_OUT, PREVIOUS_OUT

    integer status
    
!!$    IF(.not. Associated(DIMPIN)) CALL FATAL_ERROR("THE INPUT ARGUMENT MUST BE&
!!$         & ASSOCIAED FOR COPY_DIM_LIST?")
!!$
!!$    IF(.not. Associated(DIMPOUT)) CALL FATAL_ERROR("THE OUTPUT ARGUMENT MUST BE&
!!$         & ASSOCIAED FOR COPY_DIM_LIST?")

    PREVIOUS_IN => DIMPIN
    CURRENT_IN => PREVIOUS_IN%NEXT

!    DIMPOUT => NEW_DIMP()
    PREVIOUS_OUT => DIMPOUT
    CURRENT_OUT => PREVIOUS_OUT%NEXT

    ! I DON'T THINK THAT CURRENT_OUT SERVES ANY PURPOSE HERE?
    ! I DON'T THINK THAT PREVIOUS_IN SERVES ANY PURPOSE HERE?
    
    
    DO
       IF(.NOT. ASSOCIATED(CURRENT_IN)) THEN
          RETURN ! END OF DIM IN LIST
       END IF
       
       IF(.NOT. ASSOCIATED(CURRENT_IN%DIM)) THEN
          CALL FATAL_ERROR("COPY_DIM_LIST: FOUND NULL DIM POINTER IN THE LIST")
       END IF
       
       PREVIOUS_OUT%NEXT => NEW_DIMP()
       PREVIOUS_OUT%NEXT%DIM => COPY_DIM(CURRENT_IN%DIM)
       PREVIOUS_OUT%NEXT%NEXT => CURRENT_OUT
       
       ! INCRIMENT THROUGH THE LIST
       PREVIOUS_OUT => PREVIOUS_OUT%NEXT
       ! DO NOT INCRIMENT  CURRENT_OUT
       
       PREVIOUS_IN => PREVIOUS_IN%NEXT
       CURRENT_IN  => CURRENT_IN%NEXT
       
    END DO
    


  END SUBROUTINE COPY_DIM_LIST
!====================================================================
!====================================================================
  FUNCTION COPY_DIM(DIMIN) RESULT(DIMOUT)
    IMPLICIT NONE
    TYPE(NCDIM), POINTER, INTENT(IN) :: DIMIN
    TYPE(NCDIM), POINTER :: DIMOUT
    integer status
    
    IF(.not. Associated(DIMIN)) CALL FATAL_ERROR("THE ARGUMENT MUST BE&
         & ASSOCIAED FOR COPY_DIM?")

    DIMOUT => NEW_DIM()

    DIMOUT%DIMID     = DIMIN%DIMID
    DIMOUT%DIMNAME   = DIMIN%DIMNAME
    DIMOUT%DIM       = DIMIN%DIM
    DIMOUT%UNLIMITED = DIMIN%UNLIMITED
  END FUNCTION COPY_DIM
!====================================================================
!====================================================================
  FUNCTION NEW_DIMP() RESULT(DIMP)
    IMPLICIT NONE
    TYPE(NCDIMP), POINTER :: DIMP
    integer status
    
    NULLIFY(DIMP)
    ALLOCATE(DIMP,stat=status)
    if(status/=0) CALL FATAL_ERROR("ALLOC_NCDIMP COULD NOT ALLOCATE!")
    NULLIFY(DIMP%NEXT)
    NULLIFY(DIMP%DIM)

  END FUNCTION NEW_DIMP
!====================================================================
!====================================================================
  SUBROUTINE DELETE_NCF_DIMP_BYNAME(LIST,NAME,FOUND)
    IMPLICIT NONE
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    CHARACTER(LEN=*),INTENT(IN)   :: NAME
    LOGICAL,                       INTENT(OUT)  :: FOUND
    LOGICAL                                     :: F1
    TYPE(NCDIMP)                        ,POINTER :: CURRENT, PREVIOUS
    TYPE(NCVARP)                        ,POINTER :: CURRENT_VAR
    
    CURRENT_VAR => LIST%VARS%NEXT    
    DO 
       IF(.NOT. ASSOCIATED(CURRENT_VAR)) RETURN !END OF LIST

       IF(.NOT. ASSOCIATED(CURRENT_VAR%VAR)) THEN
          CALL FATAL_ERROR("DELETE_NCF_DIMP_BYNAME: NULL VAR POINTER IN FILE LIST?")
       END IF

      CALL DELETE_DIM_LINK(CURRENT_VAR%VAR,NAME,F1)

      CURRENT_VAR => CURRENT_VAR%NEXT
   END DO


    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( NAME == CURRENT%DIM%DIMNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_DIM(CURRENT%DIM)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_NCF_DIMP_BYNAME
!====================================================================
!====================================================================
  SUBROUTINE DELETE_NCF_DIMP_BYDIMID(LIST,DIMID,FOUND)
    IMPLICIT NONE
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    INTEGER,                       INTENT(IN)   :: DIMID
    LOGICAL,                       INTENT(OUT)  :: FOUND
    LOGICAL                                     :: F1
    TYPE(NCDIMP)                        ,POINTER :: CURRENT, PREVIOUS
    TYPE(NCVARP)                        ,POINTER :: CURRENT_VAR
    
    CURRENT_VAR => LIST%VARS%NEXT    
    DO 
       IF(.NOT. ASSOCIATED(CURRENT_VAR)) RETURN !END OF LIST

       IF(.NOT. ASSOCIATED(CURRENT_VAR%VAR)) THEN
          CALL FATAL_ERROR("DELETE_NCF_DIMP_BYDIMID: NULL VAR POINTER IN FILE LIST?")
       END IF

      CALL DELETE_DIM_LINK(CURRENT_VAR%VAR,DIMID,F1)

      CURRENT_VAR => CURRENT_VAR%NEXT
   END DO



    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( DIMID == CURRENT%DIM%DIMID ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT
    CALL KILL_DIM(CURRENT%DIM)
    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_NCF_DIMP_BYDIMID
!====================================================================
!====================================================================
  SUBROUTINE DELETE_NCF_DIMP_LIST(LIST)
    IMPLICIT NONE
    TYPE(NCFILE),                  INTENT(INOUT):: LIST
    TYPE(NCDIMP)                        ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST

       ! IF FOUND DELETE IT
       PREVIOUS%NEXT => CURRENT%NEXT
       CALL KILL_DIM(CURRENT%DIM)
       DEALLOCATE(CURRENT)
       CURRENT => PREVIOUS%NEXT

    END DO
    
  END SUBROUTINE DELETE_NCF_DIMP_LIST
!====================================================================
!====================================================================
  SUBROUTINE DELETE_VAR_DIMP_BYNAME(LIST,NAME,FOUND)
    IMPLICIT NONE
    TYPE(NCVAR),                   INTENT(INOUT):: LIST
    CHARACTER(LEN=*),INTENT(IN)   :: NAME
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCDIMP)                        ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST

       IF(.NOT. ASSOCIATED(CURRENT%DIM))THEN 
          !       CALL PRINT_VAR(LIST)
          CALL FATAL_ERROR("DELETE_VAR_DIMP_BYNAME: VARIABLE HAS UNASSOCIATED DIMENSION IN LIST?")
       END IF
       
       IF( NAME == CURRENT%DIM%DIMNAME ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT

    IF(CURRENT%DIM%DIMID==-1) THEN
       call KILL_DIM(CURRENT%DIM) ! IT BELONGS TO THE VAR ONLY, NOT PART
       ! OF A FILE
    ELSE
       NULLIFY(CURRENT%DIM) ! DO NOT DEALLOCATE POINTER TO ACTUAL
       ! DIMENSION, IT BELONGS TO A FILE
    END IF

    DEALLOCATE(CURRENT)
   
  END SUBROUTINE DELETE_VAR_DIMP_BYNAME
!====================================================================
!====================================================================
! THIS METHOD WILL NOT WORK FOR VARIABLES WHICH HAVE NOT BEEN
  ! ASSIGNED TO A FILE BECAUSE THE ID WILL STILL BE -1. ONLY AFTER
  ! DIMS ARE ASSIGNED TO A VARIABLE ARE THEY GIVEN AND ID NUMBER
  SUBROUTINE DELETE_VAR_DIMP_BYDIMID(LIST,DIMID,FOUND)
    IMPLICIT NONE
    TYPE(NCVAR),                   INTENT(INOUT):: LIST
    INTEGER,                       INTENT(IN)   :: DIMID
    LOGICAL,                       INTENT(OUT)  :: FOUND
    TYPE(NCDIMP)                        ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST

       IF(.NOT. ASSOCIATED(CURRENT%DIM))THEN 
          !       CALL PRINT_VAR(LIST)
          CALL FATAL_ERROR("DELETE_VAR_DIMP_BYDIMID: VARIABLE HAS UNASSOCIATED DIMENSION IN LIST?")
       END IF
    
       IF( DIMID == CURRENT%DIM%DIMID ) THEN ! FOUND THE FILE LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

    ! IF FOUND DELETE IT
    PREVIOUS%NEXT => CURRENT%NEXT

    IF(CURRENT%DIM%DIMID==-1) THEN
       call KILL_DIM(CURRENT%DIM) ! IT BELONGS TO THE VAR ONLY, NOT PART
       ! OF A FILE
    ELSE
       NULLIFY(CURRENT%DIM) ! DO NOT DEALLOCATE POINTER TO ACTUAL
       ! DIMENSION, IT BELONGS TO A FILE
    END IF
    DEALLOCATE(CURRENT)
  END SUBROUTINE DELETE_VAR_DIMP_BYDIMID
!====================================================================
!====================================================================
  SUBROUTINE DELETE_VAR_DIMP_LIST(LIST)
    IMPLICIT NONE
    TYPE(NCVAR),                   INTENT(INOUT):: LIST
    TYPE(NCDIMP)                        ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       PREVIOUS%NEXT => CURRENT%NEXT
       IF(.NOT. ASSOCIATED(CURRENT%DIM))THEN 
!          CALL PRINT_VAR(LIST)
          CALL FATAL_ERROR("VARIABLE HAS UNASSOCIATED DIMENSION IN LIST?")
       END IF

       IF(CURRENT%DIM%DIMID==-1) THEN
          call KILL_DIM(CURRENT%DIM) ! IT BELONGS TO THE VAR ONLY, NOT PART
          ! OF A FILE
       ELSE
          NULLIFY(CURRENT%DIM) ! DO NOT DEALLOCATE POINTER, IT
          ! BELONGS TO A FILE
       END IF
       
       DEALLOCATE(CURRENT)
       CURRENT => PREVIOUS%NEXT

    END DO

  END SUBROUTINE DELETE_VAR_DIMP_LIST
!====================================================================
!====================================================================
 SUBROUTINE KILL_DIM(DIM)
   TYPE(NCDIM),POINTER :: DIM
   INTEGER STATUS
   DEALLOCATE(DIM,STAT=STATUS)
   IF(STATUS /= 0) CALL FATAL_ERROR("KILL_DIM: COULD NOT DEALLOCATE")
   NULLIFY(DIM)

 END SUBROUTINE KILL_DIM
!====================================================================
!====================================================================

! INSERT_NCF_DIMP_BYNAME IS DANGEROUS AND DOES NOT EXIST

!====================================================================
!====================================================================
  SUBROUTINE INSERT_NCF_DIMP_BYDIM(LIST,DIM,FOUND)
    ! IF FOUND DO NOT INSERT DUPLICATE, RETURN FOUND
    ! INSERT UNLIMDIM AT THE END OF THE LIST
    ! INSERT NONE UNLIMDIM IN THE ORDER ADDED  
    IMPLICIT NONE
    TYPE(NCDIM), POINTER              :: DIM
    LOGICAL,             INTENT(OUT)  :: FOUND
    TYPE(NCFILE),        INTENT(INOUT):: LIST
    TYPE(NCDIMP),POINTER              :: CURRENT, PREVIOUS
    INTEGER CNT
    IF(.NOT.ASSOCIATED(DIM))&
         & CALL FATAL_ERROR("INSERT_NCF_DIMP_BYDIM: DIM NOT ASSOCIATED!")

    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.
    CNT = 1

!!$    IF (DIM%UNLIMITED) THEN ! ADD AT END OF LIST!
    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST, ADD DIM
       
       IF( DIM%DIMNAME == CURRENT%DIM%DIMNAME) THEN
          ! PROTECT AGAINST USER ERROR - MISMATCHED DIMENSIONS!
          IF(DIM%DIM .NE. CURRENT%DIM%DIM) &
               & CALL FATAL_ERROR("ATEMPTED TO ADD DIMENSION NAMED:"//TRIM(DIM%DIMNAME),&
               & "BUT THAT DIMENSION NAME ALREADY EXISTS WITH A DIFFERENT SIZE")
          
          ! PROTECT AGAINST USER ERROR - MISMATCHED UNLIMITED DIMENSIONS!
          IF(DIM%UNLIMITED .AND. .NOT. CURRENT%DIM%UNLIMITED)&
               & CALL FATAL_ERROR("ATEMPTED TO ADD DIMENSION NAMED:&
               &"//TRIM(DIM%DIMNAME)//"; AS UNLIMITED",&
               & "BUT THAT DIMENSION NAME ALREADY EXISTS AS NOT UNLIMITED")

          IF(.NOT. DIM%UNLIMITED .AND. CURRENT%DIM%UNLIMITED)&
               & CALL FATAL_ERROR("ATEMPTED TO ADD DIMENSION NAMED:&
               &"//TRIM(DIM%DIMNAME)//"; AS NOT UNLIMITED",&
               & "BUT THAT DIMENSION NAME ALREADY EXISTS AS UNLIMITED")
          
          FOUND = .TRUE.
          RETURN ! DIMENSION ALREADY EXISTS
       ELSE IF(DIM%UNLIMITED .AND. CURRENT%DIM%UNLIMITED) THEN 
          CALL FATAL_ERROR("ATTEMPT TO PUT A SECOND UNLIMITED DIMENSIO&
               &N IN THE FILE OBJECT","DIMENSION NAME: "//TRIM(DIM%DIMNAME))
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
          CNT = CNT +1
       END IF
    END DO
    ! NOT FOUND - ADD NEW DIM TO END OF LIST
    
    PREVIOUS%NEXT => NEW_DIMP()
    PREVIOUS%NEXT%DIM => DIM
    PREVIOUS%NEXT%NEXT => CURRENT
    
    DIM%DIMID = CNT ! SET THE DIMID OF THE NEW DIMENSION

    IF(DIM%UNLIMITED) LIST%UNLIMDIMID = DIM%DIMID
    
!!$ ELSE ! NOT AN UNLIMITED DIMENSION - ADD BEFORE ANY UNLIMITED DIMENSION
!!$    ! DO NOT MAKE DUPLICATE ENTRIES IN THE LIST
!!$    DO
!!$          IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST, ADD DIM TO END OF LIST
!!$          
!!$          IF( DIM%DIMNAME == CURRENT%DIM%DIMNAME) THEN
!!$             ! PROTECT AGAINST USER ERROR - MISMATCHED DIMENSIONS!
!!$             IF(DIM%DIM .NE. CURRENT%DIM%DIM) &
!!$                  & CALL FATAL_ERROR("ATEMPTED TO ADD DIMENSION NAMED:"//TRIM(DIM%DIMNAME),&
!!$                  & "BUT THAT DIMENSION NAME ALREADY EXISTS WITH A DIFFERENT SIZE")
!!$
!!$             ! PROTECT AGAINST USER ERROR - MISMATCHED UNLIMITED DIMENSIONS!
!!$             IF(CURRENT%DIM%UNLIMITED)&
!!$                  & CALL FATAL_ERROR("ATEMPTED TO ADD DIMENSION NAMED:&
!!$                  &"//TRIM(DIM%DIMNAME)//"; AS NOT UNLIMITED",&
!!$                  & "BUT THAT DIMENSION NAME ALREADY EXISTS AS UNLIMITED")
!!$
!!$             FOUND = .TRUE.
!!$             RETURN ! DIMENSION ALREADY EXISTS
!!$          ELSE IF(CURRENT%DIM%UNLIMITED) THEN ! FOUND THE UNLIMITED DIMENSION
!!$             IF(ASSOCIATED(CURRENT%NEXT)) &
!!$                  & CALL FATAL_ERROR("FOUND EXTRA DIMENSION LINK AFTER&
!!$                  & THE UNLIMITED DIMENSION", "WHILE ADDING NEW DIMENSION: "&
!!$                  &//TRIM(DIM%DIMNAME) )
!!$             
!!$             CURRENT%DIM%DIMID = CNT+1 ! INCRIMENT THE UNLIMITED DIMENSIONS DIMID
!!$             EXIT              ! Add the new dimension
!!$
!!$          ELSE
!!$             PREVIOUS => PREVIOUS%NEXT
!!$             CURRENT  => CURRENT%NEXT
!!$             CNT = CNT+1
!!$          END IF
!!$       END DO
!!$       ! NOT FOUND - INSERT DIM INTO LIST
!!$
!!$       PREVIOUS%NEXT => NEW_DIMP()
!!$       PREVIOUS%NEXT%DIM => DIM
!!$       DIM%DIMID = CNT  ! SET THE DIMID OF THE NEW DIMENSION
!!$
!!$       PREVIOUS%NEXT%NEXT => CURRENT
!!$    END IF

  END SUBROUTINE INSERT_NCF_DIMP_BYDIM
!====================================================================
!====================================================================

! INSERT_VAR_DIMP_BYNAME IS DANGEROUS AND DOES NOT EXIST

!====================================================================
!====================================================================
  SUBROUTINE INSERT_VAR_DIMP_BYDIM(LIST,DIM,FOUND)
    ! ALLOW NETCDF TO HAVE DUPLICATE DIMENSIONS IN A VARIABLE
    ! INSERT DIMS IN THE ORDER THEY ARE ADDED
    !
    ! DO NOT TOUCH DIMID - IT IS ONLY SET WHEN THE VARIABLE'S
    ! DIMENSIONS ARE POINTED TO THE FILES DIMENSIONS

    IMPLICIT NONE
    TYPE(NCDIM), POINTER              :: DIM
    LOGICAL,             INTENT(OUT)  :: FOUND
    TYPE(NCVAR),         INTENT(INOUT):: LIST
    TYPE(NCDIMP),POINTER              :: CURRENT, PREVIOUS

    IF(.NOT.ASSOCIATED(DIM))&
         & CALL FATAL_ERROR("INSERT_NCF_DIMP_BYDIM: DIM NOT ASSOCIATED!")

    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) THEN
          EXIT !END OF LIST, ADD DIM
       ELSE
          
          IF( DIM%DIMNAME == CURRENT%DIM%DIMNAME) FOUND = .TRUE.
          
          IF(DIM%UNLIMITED .AND. CURRENT%DIM%UNLIMITED) &
               & CALL FATAL_ERROR("ATTEMPT TO PUT A SECOND UNLIMITED DIMENSIO&
               &N IN THE VARIALBE NAME:"//TRIM(LIST%VARNAME),&
               & "DIMENSION NAME: "//TRIM(DIM%DIMNAME))

          ! MAKE AN EXCEPTION FOR CHARACTER DATA?!?!
          IF(CURRENT%DIM%UNLIMITED .AND. LIST%XTYPE .NE. NF90_CHAR) &
               & CALL FATAL_ERROR("ATTEMPT TO PUT A DIMENSION AFTER THE UNLIMITED DIMENSIO&
               &N IN THE VARIALBE NAME:"//TRIM(LIST%VARNAME),&
               & "DIMENSION NAME: "//TRIM(DIM%DIMNAME),"THE USER MUST &
               &ADD DIMENSION IN THE CORRECT ORDER! (UNLIMITEDS GO LAS&
               &T IN FORTRAN ORDER)")


          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO
    ! NOT FOUND - ADD NEW DIM TO END OF LIST, EVEN IF IT ALREADY
    ! EXISTS!
    PREVIOUS%NEXT => NEW_DIMP()
    PREVIOUS%NEXT%DIM => DIM
    PREVIOUS%NEXT%NEXT => CURRENT

  END SUBROUTINE INSERT_VAR_DIMP_BYDIM
!====================================================================
!====================================================================
  FUNCTION FIND_NCF_DIM_BYNAME(LIST,NAME,FOUND) RESULT(DIM)
    IMPLICIT NONE
    TYPE(NCDIM), POINTER               :: DIM
    TYPE(NCFILE),         INTENT(IN)   :: LIST
    CHARACTER(LEN=*),INTENT(IN)    :: NAME
    LOGICAL,              INTENT(OUT)  :: FOUND
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(DIM)
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( trim(NAME) == trim(CURRENT%DIM%DIMNAME) ) THEN ! FOUND THE FILE LINK
          DIM => CURRENT%DIM
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_NCF_DIM_BYNAME
!====================================================================
!====================================================================
  FUNCTION FIND_NCF_DIM_BYDIMID(LIST,DIMID,FOUND) RESULT(DIM)
    IMPLICIT NONE
    TYPE(NCDIM), POINTER               :: DIM
    INTEGER,              INTENT(IN)   :: DIMID
    TYPE(NCFILE),         INTENT(IN):: LIST
    LOGICAL,              INTENT(OUT)  :: FOUND
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(DIM)
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( DIMID == CURRENT%DIM%DIMID ) THEN ! FOUND THE DIM LINK
          DIM => CURRENT%DIM
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_NCF_DIM_BYDIMID
!====================================================================
!====================================================================
  FUNCTION FIND_NCF_DIM_UNLIMITED(LIST,FOUND) RESULT(DIM)
    IMPLICIT NONE
    TYPE(NCDIM), POINTER               :: DIM
    TYPE(NCFILE),         INTENT(IN):: LIST
    LOGICAL,              INTENT(OUT)  :: FOUND
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(DIM)
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF(CURRENT%DIM%UNLIMITED ) THEN ! FOUND THE DIM LINK
          DIM => CURRENT%DIM
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_NCF_DIM_UNLIMITED
!====================================================================
!====================================================================
  FUNCTION HAS_UNLIMITED_NCF(LIST) RESULT(FOUND)
    IMPLICIT NONE
    TYPE(NCFILE),         INTENT(IN):: LIST
    LOGICAL                        :: FOUND
    TYPE(NCDIMP)           ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF(CURRENT%DIM%UNLIMITED ) THEN ! FOUND THE DIM LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION HAS_UNLIMITED_NCF
!====================================================================
!====================================================================
  FUNCTION FIND_VAR_DIM_BYNAME(LIST,NAME,FOUND) RESULT(DIM)
    IMPLICIT NONE
    TYPE(NCDIM), POINTER               :: DIM
    TYPE(NCVAR),          INTENT(IN)   :: LIST
    CHARACTER(LEN=*),INTENT(IN)    :: NAME
    LOGICAL,              INTENT(OUT)  :: FOUND
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(DIM)
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( TRIM(NAME) == trim(CURRENT%DIM%DIMNAME) ) THEN ! FOUND THE FILE LINK
          DIM => CURRENT%DIM
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_VAR_DIM_BYNAME
!====================================================================
!====================================================================
  FUNCTION FIND_VAR_DIM_BYDIMID(LIST,DIMID,FOUND) RESULT(DIM)
    IMPLICIT NONE
    TYPE(NCDIM), POINTER               :: DIM
    INTEGER,              INTENT(IN)   :: DIMID
    TYPE(NCVAR),          INTENT(IN)   :: LIST
    LOGICAL,              INTENT(OUT)  :: FOUND
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(DIM)
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF( DIMID == CURRENT%DIM%DIMID ) THEN ! FOUND THE DIM LINK
          DIM => CURRENT%DIM
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_VAR_DIM_BYDIMID
!====================================================================
!====================================================================
  FUNCTION FIND_VAR_DIM_UNLIMITED(LIST,FOUND) RESULT(DIM)
    IMPLICIT NONE
    TYPE(NCDIM), POINTER               :: DIM
    TYPE(NCVAR),          INTENT(IN)   :: LIST
    LOGICAL,              INTENT(OUT)  :: FOUND
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    NULLIFY(DIM)
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF(CURRENT%DIM%UNLIMITED ) THEN ! FOUND THE DIM LINK
          DIM => CURRENT%DIM
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION FIND_VAR_DIM_UNLIMITED
!====================================================================
!====================================================================
  FUNCTION HAS_UNLIMITED_VAR(LIST) RESULT(FOUND)
    IMPLICIT NONE
    TYPE(NCVAR),          INTENT(IN)   :: LIST
    LOGICAL                            :: FOUND
    TYPE(NCDIMP)              ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    FOUND = .FALSE.

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       
       IF(CURRENT%DIM%UNLIMITED ) THEN ! FOUND THE DIM LINK
          FOUND = .TRUE.
          EXIT
       ELSE
          PREVIOUS => PREVIOUS%NEXT
          CURRENT  => CURRENT%NEXT
       END IF
    END DO

  END FUNCTION HAS_UNLIMITED_VAR
!====================================================================
!====================================================================
  FUNCTION VAR_DIMIDS(LIST) RESULT(DIMIDS)
    TYPE(NCVAR), INTENT(IN) :: LIST
    INTEGER, POINTER        :: DIMIDS(:)
    TYPE(NCDIMP),POINTER    :: CURRENT, PREVIOUS
    INTEGER                 :: CNT, status, SZ

    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    ALLOCATE(DIMIDS(COUNT_DIM_LIST(LIST)),stat=status)
    if(status /= 0) CALL FATAL_ERROR("VAR_DIMIDS: Can not allocate DIMIDS")

    IF(SIZE(DIMIDS)==0) RETURN

    DIMIDS=0

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       CNT = CNT + 1       

       DIMIDS(CNT) = CURRENT%DIM%DIMID
       
       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       
    END DO

    if(CNT .NE. Count_dim_list(LIST)) CALL FATAL_ERROR&
         &("VAR_DIMIDS: THE NUMBER OF DIMENSION OBJECTS IN THE VARIABL&
         &ES LIST EXCEEDS THE VARIABLES NDIMS PROPERTY")


  END FUNCTION VAR_DIMIDS
!====================================================================
!====================================================================
  FUNCTION VAR_DIMS(LIST) RESULT(DIMS)
    TYPE(NCVAR), INTENT(IN) :: LIST
    INTEGER, POINTER        :: DIMS(:)
    TYPE(NCDIMP),POINTER    :: CURRENT, PREVIOUS
    INTEGER                 :: CNT, status

    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    ALLOCATE(DIMS(COUNT_DIM_LIST(LIST)),stat=status)
    if(status /= 0) CALL FATAL_ERROR("VAR_DIMS: Can not allocate DIMS")

    IF(SIZE(DIMS)==0) RETURN

    DIMS=0

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       CNT= CNT+1

       DIMS(CNT) = CURRENT%DIM%DIM
       
       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       
    END DO

    if(CNT .NE. Count_dim_list(LIST)) CALL FATAL_ERROR&
         &("VAR_DIMS: THE NUMBER OF DIMENSION OBJECTS IN THE VARIABL&
         &ES LIST EXCEEDS THE VARIABLES NDIMS PROPERTY")


  END FUNCTION VAR_DIMS
!====================================================================
!====================================================================
  FUNCTION MEM_DIMS(LIST) RESULT(DIMS)
    TYPE(NCVAR), INTENT(IN) :: LIST
    INTEGER, POINTER        :: DIMS(:)
    TYPE(NCDIMP),POINTER    :: CURRENT, PREVIOUS
    INTEGER                 :: CNT, status

    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
!    ALLOCATE(DIMS(COUNT_DIM_LIST(LIST)),stat=status)

    ALLOCATE(DIMS(COUNT_NONSINGLETON_DIM_LIST(LIST)),stat=status)
    if(status /= 0) CALL FATAL_ERROR("VAR_DIMS: Can not allocate DIMS")

    IF(SIZE(DIMS)==0) RETURN

    DIMS=0

    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       CNT= CNT+1

       DIMS(CNT) = CURRENT%DIM%DIM
       
       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       
    END DO

    if(CNT .NE. Count_dim_list(LIST)) CALL FATAL_ERROR&
         &("VAR_DIMS: THE NUMBER OF DIMENSION OBJECTS IN THE VARIABL&
         &ES LIST EXCEEDS THE VARIABLES NDIMS PROPERTY")


  END FUNCTION MEM_DIMS
!====================================================================
!====================================================================
  FUNCTION COUNT_NCF_NS_DIM_LIST(LIST) RESULT(CNT)
    IMPLICIT NONE
    INTEGER                            :: CNT
    TYPE(NCFILE),    INTENT(IN)        :: LIST
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST

       ! DO NOT COUNT SINGLETON DIMENSIONS
       IF(CURRENT%DIM%DIM .GT. 1) CNT = CNT + 1

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       
    END DO
    
  END FUNCTION COUNT_NCF_NS_DIM_LIST
!====================================================================
!====================================================================
  FUNCTION COUNT_VAR_NS_DIM_LIST(LIST) RESULT(CNT)
    IMPLICIT NONE
    INTEGER                            :: CNT
    TYPE(NCVAR),     INTENT(IN)        :: LIST
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST

       ! DO NOT COUNT SINGLETON DIMENSIONS
       IF(CURRENT%DIM%DIM .GT. 1) CNT = CNT + 1

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT

    END DO
    
  END FUNCTION COUNT_VAR_NS_DIM_LIST
!====================================================================
!====================================================================
  SUBROUTINE ALLOC_VAR(VAR,MYDIMS)
    IMPLICIT NONE
    TYPE(NCVAR),INTENT(INOUT), POINTER :: VAR
    INTEGER,OPTIONAL,INTENT(OUT) :: MYDIMS
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR_TMP

    INTEGER :: Ndims
    INTEGER, POINTER :: DIMS(:)
    LOGICAL :: FOUND
    INTEGER :: STATUS

    NULLIFY(DIM,DIMS)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: ALLOC_VAR"

    IF(.not. ASSOCIATED(VAR))CALL FATAL_ERROR&
            & ("ALLOC_VAR: Unassociated Var!")
    
    ! GET RID OF TIME AND SINGLETON DIMENSIONS WHEN ALLOCATING SPACE
    VAR_TMP => COPY_VAR(VAR)

    DIM => FIND_UNLIMITED(VAR_TMP,FOUND)
    IF (FOUND) THEN
       CALL DELETE_DIM_LINK(VAR_TMP,DIM%DIMID,FOUND)
    END IF

    DIMS => MEM_DIMS(VAR_TMP)

    IF(.not. Associated(DIMS)) CALL FATAL_ERROR("ALLOC_VAR: Could not allocate Dims?")
    Ndims=size(DIMS)
    IF(PRESENT(MYDIMS)) MYDIMS = NDIMS

    CALL KILL_VAR(VAR_TMP)
    

    IF(DBG_SET(DBG_SBRIO)) WRITE(IPT,*) "ALLOCATING SPACE FOR:"

    select case(VAR%XTYPE)
    case(NF90_CHAR)
       IF(DBG_SET(DBG_SBRIO)) THEN
          WRITE(IPT,*) "XYTPE       :: CHAR"
          WRITE(IPT,*) "DIMS       ::",dims
       END IF

       SELECT CASE(NDIMS)
       CASE(2)
          ALLOCATE(VAR%VEC_CHR(DIMS(2)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("VAR_DIMS: Can not allocate VEC_CHR")
          VAR%VEC_CHR=""
       CASE(1)
          ALLOCATE(VAR%SCL_CHR,stat=status)
          if(status /= 0) CALL FATAL_ERROR("VAR_DIMS: Can not allocate SCL_CHR")
          VAR%SCL_CHR=""
       CASE(0)
          CALL FATAL_ERROR("Unsupported Character data dimension: 0")
       CASE DEFAULT
          CALL FATAL_ERROR("Unsupported Character data dimension")
       END SELECT

    case(NF90_BYTE)
       WRITE(IPT,*) "XYTPE       :: BYTE"
       CALL FATAL_ERROR("No Byte Type Available")

    case(NF90_SHORT)
       WRITE(IPT,*) "XYTPE       :: SHORT"
       CALL FATAL_ERROR("No Short Type Available")

    case(NF90_INT)
       IF(DBG_SET(DBG_SBRIO)) THEN
          WRITE(IPT,*) "XYTPE       :: INT"
          WRITE(IPT,*) "dims       ::",dims
       END IF


       SELECT CASE(NDIMS)
       CASE(4)
          ALLOCATE(VAR%FDA_INT(0:DIMS(1),DIMS(2),DIMS(3),DIMS(4)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate INT_FDA")
          VAR%FDA_INT=0
       CASE(3)
          ALLOCATE(VAR%CUB_INT(0:DIMS(1),DIMS(2),DIMS(3)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate INT_CUB")
          VAR%CUB_INT=0
       CASE(2)
          ALLOCATE(VAR%ARR_INT(0:DIMS(1),DIMS(2)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate INT_ARR")
          VAR%ARR_INT=0
       CASE(1)
          ALLOCATE(VAR%VEC_INT(0:DIMS(1)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate INT_VEC")
          VAR%VEC_INT=0
       CASE(0)
          ALLOCATE(VAR%SCL_INT,stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate INT_SCL")
          VAR%SCL_INT=0
       CASE DEFAULT
          CALL FATAL_ERROR("Unsupported Integer data dimension")
       END SELECT


    case(NF90_FLOAT)
       IF(DBG_SET(DBG_SBRIO)) THEN
          WRITE(IPT,*) "XYTPE       :: FLOAT"
          WRITE(IPT,*) "dims       ::",dims
       END IF

       SELECT CASE(NDIMS)
       CASE(4)
          ALLOCATE(VAR%FDA_FLT(0:DIMS(1),DIMS(2),DIMS(3),DIMS(4)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate FLT_FDA")
          VAR%FDA_FLT=0.0_SPA
       CASE(3)
          ALLOCATE(VAR%CUB_FLT(0:DIMS(1),DIMS(2),DIMS(3)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate FLT_CUB")
          VAR%CUB_FLT=0.0_SPA
       CASE(2)
          ALLOCATE(VAR%ARR_FLT(0:DIMS(1),DIMS(2)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate FLT_ARR")
          VAR%ARR_FLT=0.0_SPA
       CASE(1)
          ALLOCATE(VAR%VEC_FLT(0:DIMS(1)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate FLT_VEC")
          VAR%VEC_FLT=0.0_SPA
       CASE(0)
          ALLOCATE(VAR%SCL_FLT,stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate FLT_SCL")
          VAR%SCL_FLT=0.0_SPA
       CASE DEFAULT
          CALL FATAL_ERROR("Unsupported Integer data dimension")
       END SELECT
       
    case(NF90_DOUBLE)
       IF(DBG_SET(DBG_SBRIO)) THEN
          WRITE(IPT,*) "XYTPE       :: DOUBLE"
          WRITE(IPT,*) "dims       ::",dims
       END IF

       SELECT CASE(NDIMS)
       CASE(4)
          ALLOCATE(VAR%FDA_DBL(0:DIMS(1),DIMS(2),DIMS(3),DIMS(4)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate DBL_FDA")
          VAR%FDA_DBL=0.0_DP
       CASE(3)
          ALLOCATE(VAR%CUB_DBL(0:DIMS(1),DIMS(2),DIMS(3)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate DBL_CUB")
          VAR%CUB_DBL=0.0_DP
       CASE(2)
          ALLOCATE(VAR%ARR_DBL(0:DIMS(1),DIMS(2)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate DBL_ARR")
          VAR%ARR_DBL=0.0_DP
       CASE(1)
          ALLOCATE(VAR%VEC_DBL(0:DIMS(1)),stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate DBL_VEC")
          VAR%VEC_DBL=0.0_DP
       CASE(0)
          ALLOCATE(VAR%SCL_DBL,stat=status)
          if(status /= 0) CALL FATAL_ERROR("ALLOC_VAR: Can not allocate DBL_SCL")
          VAR%SCL_DBL=0.0_DP
       CASE DEFAULT
          CALL FATAL_ERROR("Unsupported Integer data dimension")
       END SELECT
       
    END select
    

    DEALLOCATE(DIMS)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: ALLOC_VAR"
  END SUBROUTINE ALLOC_VAR
!====================================================================
!====================================================================
  SUBROUTINE ALLOCATE_ASSOCIATED_VARS(LIST)
    ! THIS IS A SHADY THING TO DO - I AM INTENTIONALLY ALLOCATING AN
    ! ASSOCIATED VARIABLE ASSUMING THAT IT WAS ASSOCIATED NOT
    ! ALLOCATED TO BEGIN WITH, THUS NOT CAUSING A MEMORY LEAK

    IMPLICIT NONE
    TYPE(NCFILE),    INTENT(IN)        :: LIST
    TYPE(NCVARP),POINTER :: CURRENT
    TYPE(NCVAR),POINTER :: VAR

    CURRENT => LIST%VARS%NEXT
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) THEN
          EXIT ! END OF VAR LIST
       END IF

       IF(.NOT. ASSOCIATED(CURRENT%VAR)) THEN
          CALL FATAL_ERROR("ALLOCATE_ASSOCIATED_VARS: FOUND NULL VAR POINTER IN THE LIST")
       END IF

       VAR => CURRENT%VAR

       ! INTEGERS
       IF(Associated(VAR%SCL_INT)) THEN
          ALLOCATE(VAR%SCL_INT)
          VAR%SCL_INT=0
       END IF

       IF(Associated(VAR%VEC_INT))THEN
          ALLOCATE(VAR%VEC_INT(lbound(VAR%VEC_INT,1):ubound(VAR%VEC_INT,1)))
          VAR%VEC_INT=0
       END IF

       IF(Associated(VAR%ARR_INT))THEN
          ALLOCATE(VAR%ARR_INT(lbound(VAR%ARR_INT,1):ubound(VAR%ARR_INT,1),&
               & lbound(VAR%ARR_INT,2):ubound(VAR%ARR_INT,2)) )
          VAR%ARR_INT=0
       END IF

       IF(Associated(VAR%CUB_INT))THEN
          ALLOCATE(VAR%CUB_INT(lbound(VAR%CUB_INT,1):ubound(VAR%CUB_INT,1),&
               & lbound(VAR%CUB_INT,2):ubound(VAR%CUB_INT,2),&
               & lbound(VAR%CUB_INT,3):ubound(VAR%CUB_INT,3)))
          VAR%CUB_INT=0
       END IF

       IF(Associated(VAR%FDA_INT))THEN
          ALLOCATE(VAR%FDA_INT(lbound(VAR%FDA_INT,1):ubound(VAR%FDA_INT,1),&
               & lbound(VAR%FDA_INT,2):ubound(VAR%FDA_INT,2),&
               & lbound(VAR%FDA_INT,3):ubound(VAR%FDA_INT,3),&
	       & lbound(VAR%FDA_INT,4):ubound(VAR%FDA_INT,4)))
          VAR%FDA_INT=0
       END IF

       
       !FLOATING POINT VALUES
       IF(Associated(VAR%SCL_FLT)) THEN
          ALLOCATE(VAR%SCL_FLT)
          VAR%SCL_FLT=0.0_SPA
       END IF
       IF(Associated(VAR%VEC_FLT))THEN
          ALLOCATE(VAR%VEC_FLT(lbound(VAR%VEC_FLT,1):ubound(VAR%VEC_FLT,1)))
          VAR%VEC_FLT=0.0_SPA
       END IF

       IF(Associated(VAR%ARR_FLT))THEN
          ALLOCATE(VAR%ARR_FLT(lbound(VAR%ARR_FLT,1):ubound(VAR%ARR_FLT,1),&
               & lbound(VAR%ARR_FLT,2):ubound(VAR%ARR_FLT,2)) )
          VAR%ARR_FLT=0.0_SPA
       END IF

       IF(Associated(VAR%CUB_FLT))THEN
          ALLOCATE(VAR%CUB_FLT(lbound(VAR%CUB_FLT,1):ubound(VAR%CUB_FLT,1),&
               & lbound(VAR%CUB_FLT,2):ubound(VAR%CUB_FLT,2),&
               & lbound(VAR%CUB_FLT,3):ubound(VAR%CUB_FLT,3)))
          VAR%CUB_FLT=0.0_SPA
       END IF
       
       IF(Associated(VAR%FDA_FLT))THEN
          ALLOCATE(VAR%FDA_FLT(lbound(VAR%FDA_FLT,1):ubound(VAR%FDA_FLT,1),&
               & lbound(VAR%FDA_FLT,2):ubound(VAR%FDA_FLT,2),&
               & lbound(VAR%FDA_FLT,3):ubound(VAR%FDA_FLT,3),&
               & lbound(VAR%FDA_FLT,4):ubound(VAR%FDA_FLT,4)))
          VAR%FDA_FLT=0.0_SPA
       END IF
       
       !DOUBLE PRECISION
       IF(Associated(VAR%SCL_DBL)) THEN
          ALLOCATE(VAR%SCL_DBL)
          VAR%SCL_DBL=0.0_DP
       END IF
       IF(Associated(VAR%VEC_DBL))THEN
          ALLOCATE(VAR%VEC_DBL(lbound(VAR%VEC_DBL,1):ubound(VAR%VEC_DBL,1)))
          VAR%VEC_DBL=0.0_DP
       END IF
       IF(Associated(VAR%ARR_DBL))THEN
          ALLOCATE(VAR%ARR_DBL(lbound(VAR%ARR_DBL,1):ubound(VAR%ARR_DBL,1),&
               & lbound(VAR%ARR_DBL,2):ubound(VAR%ARR_DBL,2)) )
          VAR%ARR_DBL=0.0_DP
       END IF

       IF(Associated(VAR%CUB_DBL))THEN
          ALLOCATE(VAR%CUB_DBL(lbound(VAR%CUB_DBL,1):ubound(VAR%CUB_DBL,1),&
               & lbound(VAR%CUB_DBL,2):ubound(VAR%CUB_DBL,2),&
               & lbound(VAR%CUB_DBL,3):ubound(VAR%CUB_DBL,3)))
          VAR%CUB_DBL=0.0_DP
       END IF
       
       IF(Associated(VAR%FDA_DBL))THEN
          ALLOCATE(VAR%FDA_DBL(lbound(VAR%FDA_DBL,1):ubound(VAR%FDA_DBL,1),&
               & lbound(VAR%FDA_DBL,2):ubound(VAR%FDA_DBL,2),&
               & lbound(VAR%FDA_DBL,3):ubound(VAR%FDA_DBL,3),&
               & lbound(VAR%FDA_DBL,4):ubound(VAR%FDA_DBL,4)))
          VAR%FDA_DBL=0.0_DP
       END IF
       
       ! CHARACTER ARRAY DATA
       IF(Associated(VAR%SCL_CHR)) THEN
          ALLOCATE(VAR%SCL_CHR)
          VAR%SCL_CHR=""
       END IF
       IF(Associated(VAR%VEC_CHR))THEN
          ALLOCATE(VAR%VEC_CHR(lbound(VAR%VEC_CHR,1):ubound(VAR%VEC_CHR,1)))
          VAR%VEC_CHR=""
       END IF

       CURRENT => CURRENT%NEXT

    END DO



  END SUBROUTINE ALLOCATE_ASSOCIATED_VARS
!====================================================================
!====================================================================
  FUNCTION COUNT_NCF_DIM_LIST(LIST) RESULT(CNT)
    IMPLICIT NONE
    INTEGER                            :: CNT
    TYPE(NCFILE),    INTENT(IN)        :: LIST
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       CNT = CNT + 1
    END DO
    
  END FUNCTION COUNT_NCF_DIM_LIST
!====================================================================
!====================================================================
  FUNCTION COUNT_VAR_DIM_LIST(LIST) RESULT(CNT)
    IMPLICIT NONE
    INTEGER                            :: CNT
    TYPE(NCVAR),     INTENT(IN)        :: LIST
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    
    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) RETURN !END OF LIST
       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
       CNT = CNT + 1
    END DO
    
  END FUNCTION COUNT_VAR_DIM_LIST
!====================================================================
!====================================================================
  SUBROUTINE PRINT_NCF_DIM_LIST(LIST)
    IMPLICIT NONE
    type(NCFILE), intent(IN) :: LIST
    TYPE(NCDIMP), POINTER    :: CURRENT, PREVIOUS
    INTEGER                  :: CNT
    Character(len=4)         :: chr

    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT

    IF(.NOT. ASSOCIATED(CURRENT)) THEN ! EMPTY LIST
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%%%%%% FILE DIMENSION LIST IS EMPTY %%%%%%%%%%%%%"
       RETURN
    ELSE
       if(DBG_SET(DBG_LOG)) then
          write(IPT,*)"%%%%%%% PRINTING FILE DIMENSION LIST %%%%%%%%%"
          write(IPT,*)"%%%%%%% FILE NAME: "//TRIM(LIST%FNAME)//" %%%%%%%%%"
       end if
    END IF
    
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       CNT = CNT + 1
       write(chr,'(I4.4)')CNT
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"! PRINTING DIMENSION LIST ENTRY #"//CHR
       CALL PRINT_DIM(CURRENT%DIM)

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
    END DO
    if(DBG_SET(DBG_LOG)) &
         & write(IPT,*)"%%%%%%%%%%% END OF DIMENSION LIST %%%%%%%%%%%%%"
  END SUBROUTINE PRINT_NCF_DIM_LIST
!====================================================================
!====================================================================
  SUBROUTINE PRINT_VAR_DIM_LIST(LIST)
    IMPLICIT NONE
    type(NCVAR), intent(IN)           :: LIST
    TYPE(NCDIMP)               ,POINTER :: CURRENT, PREVIOUS
    INTEGER                            :: CNT
    Character(len=4)                   :: chr

    PREVIOUS => LIST%DIMS
    CURRENT  => PREVIOUS%NEXT

    IF(.NOT. ASSOCIATED(CURRENT)) THEN ! EMPTY LIST
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%%%%%% VARIABLE DIMENSION LIST IS EMPTY %%%%%%%%%%%%%"
       RETURN
    ELSE
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"%%%%%%% PRINTING VARIABLE: "//TRIM(LIST&
            &%VARNAME)//"; DIMENSION LIST %%%%%%%%%"
    END IF
    
    CNT = 0
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT !END OF LIST
       CNT = CNT + 1
       write(chr,'(I4.4)')CNT
       if(DBG_SET(DBG_LOG)) &
            & write(IPT,*)"! PRINTING DIMENSION LIST ENTRY #"//CHR
       CALL PRINT_DIM(CURRENT%DIM)

       PREVIOUS => PREVIOUS%NEXT
       CURRENT  => CURRENT%NEXT
    END DO
    if(DBG_SET(DBG_LOG)) &
         & write(IPT,*)"%%%%%%%%%%% END OF DIMENSION LIST %%%%%%%%%%%%%"
  END SUBROUTINE PRINT_VAR_DIM_LIST
!====================================================================
!====================================================================
  SUBROUTINE PRINT_DIM(DIM)
    implicit none
    type(NCDIM),pointer, intent(IN) :: DIM
    
    if(DBG_SET(dbg_log)) then
       WRITE(IPT,*) "======== PRINT NCDIM TYPE ======="
       if(.not. associated(DIM)) then
          WRITE(IPT,*) "THIS NCDIM HAS NOT BEEN ASSOCIATED"
          WRITE(IPT,*) "======= PRINTED NCDIM TYPE ======"
          return
       end if
       WRITE(IPT,*) "DIMNAME  ::"//TRIM(DIM%DIMNAME)
       WRITE(IPT,*) "DIMID    ::",DIM%DIMID
       WRITE(IPT,*) "DIMLEN   ::",DIM%DIM
       WRITE(IPT,*) "UNLIMITED::",DIM%UNLIMITED
       WRITE(IPT,*) "======= PRINTED NCDIM TYPE ======"
    end if
  END SUBROUTINE PRINT_DIM
!====================================================================
!====================================================================
!====================== MIXED OPERATORS =============================
!====================================================================
!====================================================================
!
! THESE OPERATATIONS ADD TWO MIXED NC OBJECTS, NULLIFY THE ARGUMENTS
! AND RETURN THE LIST OF THE SUM. IN GENERAL THEY DO NOT DELETE OR 
! ALLOCATE ANY MEMORY - ONLY POINT THE LIST OF THE FIRST OBJECT TO
! THE NEW DATA IN THE MEMORY OF THE SECOND.
!     EXCEPTION: IF A DIMENSION OF A VARIABLE ALREADY EXISTS IN A
!     FILE WHEN YOU ADD IT, THAT DIMENSION IS DELETED FROM THE
!     VARIABLE AND THE VARIABLES DIMENSION IS POINTED TO THE FILE'S
!     DIMENSION  
! 
! THESE METHODS ARE 'SAFE' AND INTENDED FOR END USERS TO WORK WITH
! THE PRINT_* METHODS ARE ALSO INTENDED FOR END USERS
!
  FUNCTION VAR_PLUS_DIM(VAR,DIM) RESULT(RET)
    ! ADD THE DIMENSION TO THE VARIABLE DIM LIST IF IT DOES NOT EXIST
    IMPLICIT NONE
    type(NCVAR), POINTER        :: VAR
    type(NCVAR), POINTER        :: RET
    type(NCDIM), POINTER        :: DIM
    LOGICAL FOUND
    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START VAR_PLUS_DIM"

    IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
         & ("VAR_PLUS_DIM: THE VAR OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")

    IF(.NOT. ASSOCIATED(DIM)) CALL FATAL_ERROR &
         & ("VAR_PLUS_DIM: THE DIM OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")

    IF(VAR%CONNECTED) CALL FATAL_ERROR&
         ("CAN NOT ADD DIMENSIONS TO A VARIABLE OBJECT ALREADY IN A FILE")

    CALL INSERT_DIM_LINK(VAR,DIM,FOUND)
    

!     IF(FOUND) CALL FATAL_ERROR("ERROR ADDIND DIMENSION TO VARIABLE",&
!         & "THE DIMENSION: "//TRIM(DIM%DIMNAME)//"; ALREADY EXISTS",&
!         & "IN THE VARIABLE: "//TRIM(VAR%VARNAME))


    !ALL MEMORY NOW BELONGS TO RET
    RET => VAR
    NULLIFY(VAR)
    NULLIFY(DIM) 
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "END VAR_PLUS_DIM"
  END FUNCTION VAR_PLUS_DIM
!====================================================================
!====================================================================
  FUNCTION VAR_PLUS_ATT(VAR,ATT) RESULT(RET)
    ! ADD THE ATTRIBUTE TO THE VARIABLE ATT LIST IF IT DOES NOT EXIST
    IMPLICIT NONE
    type(NCVAR), POINTER        :: VAR
    type(NCATT), POINTER        :: ATT
    type(NCVAR), POINTER        :: RET
    LOGICAL FOUND
    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START VAR_PLUS_ATT"

    IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
         & ("VAR_PLUS_ATT: THE VAR OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")

    IF(.NOT. ASSOCIATED(ATT)) CALL FATAL_ERROR &
         & ("VAR_PLUS_ATT: THE ATT OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")

    IF(VAR%CONNECTED) CALL FATAL_ERROR&
         ("CAN NOT ADD ATTRIBUTES TO A VARIABLE OBJECT ALREADY IN A FILE")

    CALL INSERT_ATT_LINK(VAR,ATT,FOUND)
    ! ADD ERROR HANDLING LATER FOR MERGING ATTS
    IF(FOUND) CALL FATAL_ERROR("ERROR ADDIND ATTRIBUTE TO VARIABLE",&
         & "THE ATTRIBUTE: "//TRIM(ATT%ATTNAME)//"; ALREADY EXISTS",&
         & "IN THE VARIABLE: "//TRIM(VAR%VARNAME))


    !ALL MEMORY NOW BELONGS TO RET
    RET => VAR
    NULLIFY(VAR)
    NULLIFY(ATT)
    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END VAR_PLUS_ATT"
  END FUNCTION VAR_PLUS_ATT
!====================================================================
!====================================================================
  FUNCTION NCF_PLUS_VAR(NCF,VAR) RESULT(RET)
    ! ADD THE VARIABLE TO THE FILE LIST IF IT DOES NOT ALREADY EXIST
    !
    ! ADD THE VARIABLES DIMENSIONS IF THEY DO NOT EXIST IN THE FILE LIST
    ! OR POINT TO THE FILES DIMENSIONS AND USE THE FILES DIMIDS IF THE
    ! DIMENSION ALREADY EXISTS IN THE FILE'S LIST
    !
    ! ATTRIBUTES BELONG TO THE VARIABLE, JUST ADD THEM WITH THE VARIABLE
    IMPLICIT NONE
    type(NCVAR),  POINTER       :: VAR
    type(NCFILE), POINTER       :: RET
    type(NCFILE), POINTER       :: NCF
    type(NCDIM),  POINTER       :: VDIM, GDIM
    type(NCDIMP), POINTER       :: CURRENT
    
    LOGICAL :: FOUND, F
    CHARACTER(LEN=NF90_MAX_NAME+1):: NAME
    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NCF_PLUS_VAR"

    IF(.NOT. ASSOCIATED(NCF)) CALL FATAL_ERROR &
         & ("NCF_PLUS_VAR: THE FILE OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")

    IF(.NOT. ASSOCIATED(VAR)) CALL FATAL_ERROR &
         & ("NCF_PLUS_VAR: THE VAR OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")


    IF(NCF%CONNECTED) CALL FATAL_ERROR&
         ("CAN NOT ADD VARIABLE TO FILE OBJECTS ALREADY CONNECTED TO A NETCDF FILE")

    IF(VAR%CONNECTED) CALL FATAL_ERROR&
         ("CAN NOT ADD VARIABLE TO A FILE OBJECT WHEN IT IS ALREADY CONNECTED")


    ! LOOK TO SEE IF WE NEED TO ADD ANY NEW DIMENSIONS TO THE FILE
    CURRENT=>VAR%DIMS%NEXT
    DO
       IF(.NOT. ASSOCIATED(CURRENT)) EXIT

       IF (.NOT. ASSOCIATED(CURRENT%DIM)) &
            & CALL FATAL_ERROR("NCF_PLUS_VAR: ATTEMPT TO ADD VARIABLE &
            &THAT HAS UNASSOCIATED DIM POINTERS IN ITS LINK LIST!")

       VDIM => CURRENT%DIM

       CALL INSERT_DIM_LINK(NCF,VDIM,FOUND)
       ! THIS WILL ADD THE VARIABLES DIMENSION IF IT IS NOT THERE ALREADY
       ! IN THE FILE. THERE IS SOME BASIC CONSISTANCY CHECKING.

       ! IF WE JUST ADDED IT, GREAT! NOW THE FILE'S DIM LINK LIST AND
       ! THE VARIABLE'S DIM LINK LIST ALREADY POINT TO THE SAME MEMORY
       
       IF(FOUND) THEN ! THE DIMENSION IS ALREADY THERE. 

          NAME = VDIM%DIMNAME
          GDIM => FIND_DIM(NCF,NAME,F)
          IF(.NOT. F) CALL FATAL_ERROR("NCF_PLUS_VAR: CAN'T FIND D&
               &IMENSION BUT IT WAS THERE A MINUTE AGO?")

          IF (.not. associated(VDIM,target = GDIM)) THEN ! CHECKS TO SEE IF IT IS
             ! POINTING TO THE SAME MEMORY OR TRUELY A DUPLICATE
             
             ! DELETE THE DUPLICATE AND POINT TO THE FILES DIM LINK LIST
             
             ! THIS ALLOWS THE FILE TO SET THE DIMIDS. DIMIDS SET IN THE
             ! VARIABLES DIM LIST ARE IGNORED
             CALL KILL_DIM(VDIM)

             ! MOVE NOW EMPTY VARIABLE'S DIM LINK LIST POINTER TO FILE'S DIM LINK LIST
             CURRENT%DIM => GDIM
          END IF
       END IF
       CURRENT => CURRENT%NEXT
    END DO

    ! NOW ADD THE VARIABLE TO THE FILE

    CALL INSERT_VAR_LINK(NCF,VAR,FOUND)
    IF(FOUND) CALL FATAL_ERROR("NCF_PLUS_VAR: THIS VARIABLE ALREADY EX&
         &ISTS IN THE FILE. YOU CAN'T ADD IT AGAIN!", &
         & "VARIABLE NAME: "//TRIM(VAR%VARNAME))


!    IF(ASSOCIATED(VAR%NCID)) THEN
      ! IF(VAR%NCID /= NCF%NCID) &
!            CALL FATAL_ERROR&
!            &("CAN'T ADD VARIABLE TO A FILE WITH A DIFFERENT NCID")
!    ELSE
    
!    END IF

    VAR%NCID => NCF%NCID
    VAR%CONNECTED=.TRUE.

    !ALL MEMORY NOW BELONGS TO RET
    RET => NCF
    NULLIFY(NCF)
    NULLIFY(VAR)
    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NCF_PLUS_VAR"

  END FUNCTION NCF_PLUS_VAR
!====================================================================
!====================================================================
  FUNCTION NCF_PLUS_ATT(NCF,ATT) RESULT(RET)
    IMPLICIT NONE
    type(NCATT),  POINTER       :: ATT
    type(NCFILE), POINTER       :: RET
    type(NCFILE), POINTER       :: NCF
    LOGICAL FOUND
    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NCF_PLUS_ATT"

    IF(.NOT. ASSOCIATED(NCF)) CALL FATAL_ERROR &
         & ("NCF_PLUS_ATT: THE FILE OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")

    IF(.NOT. ASSOCIATED(ATT)) CALL FATAL_ERROR &
         & ("NCF_PLUS_ATT: THE ATT OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")

    IF(NCF%CONNECTED) CALL FATAL_ERROR&
         ("CAN NOT ADD ATTRIBUTE TO FILE OBJECTS ALREADY CONNECTED TO A NETCDF FILE")

    CALL INSERT_ATT_LINK(NCF,ATT,FOUND)
    ! ADD ERROR HANDLING LATER FOR MERGING ATTS
    IF(FOUND) CALL FATAL_ERROR("ERROR ADDIND ATTRIBUTE TO FILE",&
         & "THE ATTRIBUTE: "//TRIM(ATT%ATTNAME)//"; ALREADY EXISTS",&
         & "IN THE FILE: "//TRIM(NCF%FNAME))


    !ALL MEMORY NOW BELONGS TO RET
    RET => NCF
    NULLIFY(NCF)
    NULLIFY(ATT)
    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NCF_PLUS_ATT"
  END FUNCTION NCF_PLUS_ATT
!====================================================================
!====================================================================
  FUNCTION NCF_PLUS_DIM(NCF,DIM) RESULT(RET)
    IMPLICIT NONE
    type(NCDIM),  POINTER       :: DIM
    type(NCFILE), POINTER       :: RET
    type(NCFILE), POINTER       :: NCF
    LOGICAL FOUND
    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NCF_PLUS_DIM"

    IF(.NOT. ASSOCIATED(NCF)) CALL FATAL_ERROR &
         & ("NCF_PLUS_DIM: THE FILE OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")

    IF(.NOT. ASSOCIATED(DIM)) CALL FATAL_ERROR &
         & ("NCF_PLUS_DIM: THE DIM OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")


    IF(NCF%CONNECTED) CALL FATAL_ERROR&
         ("CAN NOT ADD DIMENSION TO FILE OBJECTS ALREADY CONNECTED TO A NETCDF FILE")

    CALL INSERT_DIM_LINK(NCF,DIM,FOUND)

    ! ADDING DUPLICATE DIMENSIONS TO A FILE IS NOT IS OKAY
    IF(FOUND) CALL FATAL_ERROR("ERROR ADDIND DIMENSION TO FILE",&
         & "THE DIMENSION: "//TRIM(DIM%DIMNAME)//"; ALREADY EXISTS",&
         & "IN THE FILE: "//TRIM(NCF%FNAME))


    !ALL MEMORY NOW BELONGS TO RET
    RET => NCF
    NULLIFY(NCF)
    NULLIFY(DIM)

    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NCF_PLUS_DIM"
  END FUNCTION NCF_PLUS_DIM
!====================================================================
!====================================================================
  FUNCTION NCF_PLUS_NCF(NCF1,NCF2) RESULT(RET)
    ! ADD EACH VARIABLE FROM NCF2 TO NCF1
    ! THE VARIABLES WILL TAKE THEIR DIMENSIONS WITH THEM AS NEEDED
    !
    ! ADD ANY ATTRIBUTES THAT DO NOT ALREADY EXIST IN NCF1

    IMPLICIT NONE
    type(NCFILE), POINTER       :: NCF1
    type(NCFILE), POINTER       :: NCF2
    type(NCFILE), POINTER       :: RET

    type(NCVARP), POINTER       :: CURRENT_VAR
    type(NCVAR), POINTER        :: VAR

    type(NCDIMP), POINTER       :: CURRENT_DIM
    type(NCDIM), POINTER        :: DIM1
    type(NCDIM), POINTER        :: DIM2

    TYPE(NCATTP), POINTER       :: CURRENT_ATT
!!$
    TYPE(NCATT),  POINTER       :: ATTT
!!$    

    CHARACTER(LEN=NF90_MAX_NAME+1):: NAME
    LOGICAL FOUND1, FOUND2

    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "START NCF_PLUS_NCF"


    IF(.NOT. ASSOCIATED(NCF1)) CALL FATAL_ERROR &
         & ("NCF_PLUS_NCF: THE FIRST FILE OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")

    IF(.NOT. ASSOCIATED(NCF2)) CALL FATAL_ERROR &
         & ("NCF_PLUS_NCF: THE SECOND FILE OBJECT ARGUMENT IS NOT ASSOC&
         &IATED. THIS IS ILLEGAL.","THE DEPARTMENT OF FVCOM SECURITY HAS BEEN CONTACTED")
    
    IF(NCF1%CONNECTED .OR. NCF2%CONNECTED) CALL FATAL_ERROR&
         ("CAN NOT COMBINE FILE OBJECTS ALREADY CONNECTED TO A NETCDF FILE")
    
    CURRENT_VAR => NCF2%VARS%NEXT

    DO
       IF(.NOT. ASSOCIATED(CURRENT_VAR)) EXIT
       
       IF (.NOT. ASSOCIATED(CURRENT_VAR%VAR)) &
            & CALL FATAL_ERROR("NCF_PLUS_NCF: ATTEMPT TO ADD A FILE &
            &THAT HAS UNASSOCIATED VAR POINTER IN ITS LINK LIST!")
       
       ! Variable from the second file
       VAR =>Copy_Var(CURRENT_VAR%VAR)
       !DISCONECT THE VARIABLE MOMENTARILY
       VAR%CONNECTED=.false.


       CURRENT_DIM=>VAR%DIMS%NEXT
       DO
          IF(.NOT. ASSOCIATED(CURRENT_DIM)) EXIT
          
          IF (.NOT. ASSOCIATED(CURRENT_DIM%DIM)) &
               & CALL FATAL_ERROR("NCF_PLUS_NCF: ATTEMPT TO ADD VARIABLE &
               &THAT HAS UNASSOCIATED DIM POINTERS IN ITS LINK LIST!")
          
          !DIM from the variable of the second file
          DIM2 => CURRENT_DIM%DIM
          
          CALL INSERT_DIM_LINK(NCF1,DIM2,FOUND1)
          ! THIS WILL ADD THE VARIABLES DIMENSION IF IT IS NOT THERE ALREADY
          ! IN THE FILE. THERE IS SOME BASIC CONSISTANCY CHECKING.
          
          ! IF WE JUST ADDED IT, GREAT! NOW THE FILE'S DIM LINK LIST AND
          ! THE VARIABLE'S DIM LINK LIST ALREADY POINT TO THE SAME MEMORY
          
          ! IF WE FOUND IT IN THE FILE TO WHICH THE VARIABLE IS BEING
          ! ADDED MOVE THE POINTER OVER
          IF (FOUND1) THEN
             CURRENT_DIM%DIM => FIND_DIM(NCF1,DIM2%DIMNAME,FOUND2)
             IF (.NOT. FOUND2) CALL FATAL_ERROR("NCF_PLUS_NCF: CAN'T FIN&
                  &D DIMENSION IN FILE BUT IT WAS HERE A SECOND AGO?")
             CALL KILL_DIM(DIM2)

          END IF
          
          CURRENT_DIM => CURRENT_DIM%NEXT
       END DO

       ! NOW ADD THE VARIABLE TO THE FILE
       
       CALL INSERT_VAR_LINK(NCF1,VAR,FOUND1)
       IF(FOUND1) CALL FATAL_ERROR("NCF_PLUS_VAR: THIS VARIABLE ALREADY EX&
            &ISTS IN THE FILE. YOU CAN'T ADD IT AGAIN!", &
            & "VARIABLE NAME: "//TRIM(VAR%VARNAME))
       
       
       VAR%NCID => NCF1%NCID
       VAR%CONNECTED=.TRUE.
       
       CURRENT_VAR => CURRENT_VAR%NEXT
    END DO
    

!!$    ! NOW GO THROUGH DIMENSION LIST IN BOTH FILES AND DELETE DUPLICATES!
!!$    CURRENT_DIM=>NCF2%DIMS%NEXT
!!$    DO
!!$       IF(.NOT. ASSOCIATED(CURRENT_DIM)) EXIT
!!$       
!!$       IF (.NOT. ASSOCIATED(CURRENT_DIM%DIM)) &
!!$            & CALL FATAL_ERROR("NCF_PLUS_NCF: ATTEMPT TO ADD VARIABLE &
!!$            &THAT HAS UNASSOCIATED DIM POINTERS IN ITS LINK LIST!")
!!$              
!!$       !DIM from the second file
!!$       DIM2 => CURRENT_DIM%DIM
!!$       
!!$       DIM1 => FIND_DIM(NCF1,DIM2%DIMNAME,FOUND1)
!!$       IF (.NOT. FOUND1) THEN
!!$          CALL FATAL_ERROR("NCF_PLUS_NCF: ALL DIMENSIO&
!!$               &NS FROM FILE TWO SHOULD ALREADY BE ADDED TO FILE ONE:",&
!!$               &"FOUND DIMNAME: "//DIM2%DIMNAME//" IN FILE 2 BUT NOT IN 1&
!!$               & ?")
!!$       ELSE
!!$          IF(.not. associated(DIM1,target = DIM2)) CALL KILL_DIM(DIM2)
!!$       END IF
!!$
!!$       CURRENT_DIM => CURRENT_DIM%NEXT
!!$
!!$    END DO


   CURRENT_ATT => NCF2%ATTS%NEXT
    DO
       IF(.NOT. ASSOCIATED(CURRENT_ATT)) EXIT
 
       IF (.NOT. ASSOCIATED(CURRENT_ATT%ATT)) &
            & CALL FATAL_ERROR("NCF_PLUS_NCF: ATTEMPT TO ADD A FILE &
            &THAT HAS UNASSOCIATED ATT POINTER IN ITS LINK LIST!")

       ATTT => Copy_Att(CURRENT_ATT%ATT)
       NCF1 => NCF_PLUS_ATT(NCF1,ATTT)
!!$       NCF1 => NCF_PLUS_ATT(NCF1,Copy_Att(CURRENT_ATT%ATT))
       !ALREADY INCRIMENTED NATTS

       CURRENT_ATT => CURRENT_ATT%NEXT
    END DO



    RET => NCF1
    NULLIFY(NCF1)
    CALL KILL_FILE(NCF2)
    NULLIFY(NCF2)

    if(DBG_SET(dbg_sbr)) &
       & write(IPT,*) "END NCF_PLUS_NCF"

  END FUNCTION NCF_PLUS_NCF
!====================================================================
!====================================================================
  FUNCTION NCFLIST_PLUS_NCF(NCFLIST,NCF) RESULT(RET)
    IMPLICIT NONE
    type(NCFILELIST), POINTER   :: NCFLIST
    type(NCFILELIST), POINTER   :: RET
    type(NCFILE),     POINTER   :: NCF
    LOGICAL FOUND

    IF (.NOT. ASSOCIATED(NCFLIST))THEN
       ALLOCATE(NCFLIST)
       NCFLIST%FIRST => NEW_FILEP()
    END IF

    CALL INSERT_FILE_LINK(NCFLIST,NCF,FOUND)
    IF(FOUND) CALL FATAL_ERROR("NCFLIST_PLUS_NCF: THIS FILE ALREADY EX&
         &ISTS IN THE FILELIST. YOU CAN'T ADD IT AGAIN!", &
         & "FILE NAME: "//TRIM(NCF%FNAME))

    !ALL MEMORY IN NCFLIST NOW BELONGS TO RET
    RET => NCFLIST
    NULLIFY(NCFLIST)
    NULLIFY(NCF)
  END FUNCTION NCFLIST_PLUS_NCF



  
END MODULE MOD_NCLL
