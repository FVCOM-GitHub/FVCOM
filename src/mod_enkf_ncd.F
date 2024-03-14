!!$/=============================================================================/
!!$Copyright (c) 2007, The University of Massachusetts Dartmouth 
!!$Produced at the School of Marine Science & Technology 
!!$Marine Ecosystem Dynamics Modeling group
!!$All rights reserved.
!!$
!!$The FVCOM Offline Lagrangian Model has been developed by the joint UMASSD-WHOI
!!$research team.   For details of authorship and attribution of credit please see
!!$the FVCOM technical manual or contact the MEDM group.
!!$
!!$ 
!!$This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu/ The
!!$full copyright notice is contained in the file COPYRIGHT located in the root
!!$directory of the FVCOM code. This original header must be maintained in all
!!$distributed versions.
!!$
!!$THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!!$ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
!!$IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
!!$ARE  DISCLAIMED.  
!!$
!!$/-----------------------------------------------------------------------------/
!!$CVS VERSION INFORMATION
!!$$Id$
!!$$Name$
!!$$Revision$
!!$/=============================================================================/
!==============================================================================|
!   NCD UTILITIES                                                              !
!==============================================================================|

MODULE MOD_NCD   

  USE NETCDF
  USE MOD_PREC
  IMPLICIT NONE
  SAVE

  INTEGER :: NC_FID


CONTAINS
  !==============================================================================!
  INTEGER FUNCTION GETDIM(FID,SSIZE,DIMNAME)
    !==============================================================================!
    !  Read dimensions
    !==============================================================================!

    IMPLICIT NONE
    INTEGER, INTENT(IN)              :: FID
    INTEGER, INTENT(IN)              :: SSIZE
    CHARACTER(LEN=SSIZE), INTENT(IN) :: DIMNAME
    INTEGER                          :: LENGTH
    INTEGER                          :: IERR
    INTEGER                          :: DIMID
    CHARACTER(LEN=NF90_MAX_NAME)     :: TEMPNAME

    IERR = NF90_INQ_DIMID(FID,TRIM(DIMNAME),DIMID)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'Error getting dimension id: ',TRIM(DIMNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    IERR = NF90_INQUIRE_DIMENSION(FID,DIMID,TEMPNAME,LENGTH)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'Error getting dimension: ',TRIM(DIMNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    GETDIM = LENGTH

  END FUNCTION GETDIM
  !==============================================================================!

  !==============================================================================!
  SUBROUTINE GETSVAR(FID,NLEN,VARNAME,I1,I2,TEMP)
    !============================================================================!
    !  Read Static variables
    !==============================================================================!

    IMPLICIT NONE
    INTEGER, INTENT(IN)             :: FID
    INTEGER, INTENT(IN)             :: NLEN
    CHARACTER(LEN=NLEN), INTENT(IN) :: VARNAME
    INTEGER, INTENT(IN)             :: I1,I2
    REAL(SP),INTENT(OUT)            :: TEMP(I1,I2)
    INTEGER                         :: IERR
    INTEGER                         :: VARID
    INTEGER, ALLOCATABLE            :: DIMS(:)

    IF(I2 == 1)THEN
       ALLOCATE(DIMS(1))
       DIMS(1) = 1
    ELSE
       ALLOCATE(DIMS(2))
       DIMS(1) = 1
       DIMS(2) = 1
    END IF

    IERR = NF90_INQ_VARID(FID,TRIM(VARNAME),VARID)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'error getting variable id: ',TRIM(VARNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    IERR = NF90_GET_VAR(FID,VARID,TEMP,DIMS)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'1 error getting variable: ',TRIM(VARNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    RETURN
  END SUBROUTINE GETSVAR
  !==============================================================================!

  !==============================================================================!
  SUBROUTINE GETSVAR_D(FID,NLEN,VARNAME,I1,I2,TEMP)
    !============================================================================!
    !  Read Static variables
    !==============================================================================!

    IMPLICIT NONE
    INTEGER, INTENT(IN)             :: FID
    INTEGER, INTENT(IN)             :: NLEN
    CHARACTER(LEN=NLEN), INTENT(IN) :: VARNAME
    INTEGER, INTENT(IN)             :: I1,I2
    REAL(DP),INTENT(OUT)            :: TEMP(I1,I2)
    INTEGER                         :: IERR
    INTEGER                         :: VARID
    INTEGER, ALLOCATABLE            :: DIMS(:)

    IF(I2 == 1)THEN
       ALLOCATE(DIMS(1))
       DIMS(1) = 1
    ELSE
       ALLOCATE(DIMS(2))
       DIMS(1) = 1
       DIMS(2) = 1
    END IF

    IERR = NF90_INQ_VARID(FID,TRIM(VARNAME),VARID)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'error getting variable id: ',TRIM(VARNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    IERR = NF90_GET_VAR(FID,VARID,TEMP,DIMS)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'1 error getting variable: ',TRIM(VARNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    RETURN
  END SUBROUTINE GETSVAR_D
  !==============================================================================!

  !==============================================================================!
  SUBROUTINE GETDVAR(FID,NLEN,VARNAME,I1,I2,TEMP,NT)
    !============================================================================!
    !  Read time dynamic variables
    !==============================================================================!

    IMPLICIT NONE
    INTEGER, INTENT(IN)             :: FID
    INTEGER, INTENT(IN)             :: NLEN
    CHARACTER(LEN=NLEN), INTENT(IN) :: VARNAME
    INTEGER, INTENT(IN)             :: I1,I2
    REAL(SP),  INTENT(OUT)          :: TEMP(I1,I2)
    INTEGER                         :: IERR
    INTEGER                         :: VARID
    INTEGER                         :: NT
    INTEGER, ALLOCATABLE            :: DIMS(:)

    IF(I2 == 1)THEN
       ALLOCATE(DIMS(2))
       DIMS(1) = 1
       DIMS(2) = NT 
    ELSE
       ALLOCATE(DIMS(3))
       DIMS(1) = 1
       DIMS(2) = 1
       DIMS(3) = NT 
    END IF

    IERR = NF90_INQ_VARID(FID,TRIM(VARNAME),VARID)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'error getting variable id: ',TRIM(VARNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    IERR = NF90_GET_VAR(FID,VARID,TEMP,DIMS)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'2 error getting variable: ',TRIM(VARNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    RETURN
  END SUBROUTINE GETDVAR
  !==============================================================================!

  !==============================================================================!
  SUBROUTINE PUTDVAR(FID,NLEN,VARNAME,I1,TEMP,NT)
    !============================================================================!
    !  Write dynamic time variables
    !==============================================================================!

    IMPLICIT NONE
    INTEGER, INTENT(IN)             :: FID
    INTEGER, INTENT(IN)             :: NLEN
    CHARACTER(LEN=NLEN), INTENT(IN) :: VARNAME
    INTEGER, INTENT(IN)             :: I1
    REAL(SP),INTENT(IN)             :: TEMP(I1)
    INTEGER                         :: IERR
    INTEGER                         :: VARID
    INTEGER, ALLOCATABLE            :: DIMS(:)
    INTEGER                         :: NT

    IF(I1 == 1)THEN
       ALLOCATE(DIMS(1))
       DIMS(1) = NT 
    ELSE
       ALLOCATE(DIMS(2))
       DIMS(1) = 1
       DIMS(2) = NT 
    END IF

    IERR = NF90_INQ_VARID(FID,TRIM(VARNAME),VARID)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'error getting variable id: ',TRIM(VARNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    IERR = NF90_PUT_VAR(FID,VARID,TEMP,DIMS)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'3 error getting variable: ',TRIM(VARNAME)
       print *, temp,dims ! should be mark
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    RETURN
  END SUBROUTINE PUTDVAR
  !==============================================================================!

  !==============================================================================!
  SUBROUTINE PUTSVAR(FID,NLEN,VARNAME,I1,TEMP)
    !============================================================================!
    !  Write static variables
    !==============================================================================!

    IMPLICIT NONE
    INTEGER, INTENT(IN)             :: FID
    INTEGER, INTENT(IN)             :: NLEN
    CHARACTER(LEN=NLEN), INTENT(IN) :: VARNAME
    INTEGER, INTENT(IN)             :: I1
    REAL(SP),  INTENT(IN)           :: TEMP(I1)
    INTEGER                         :: IERR
    INTEGER                         :: VARID
    INTEGER, DIMENSION(1)           :: DIMS

    DIMS(1)=1

    IERR = NF90_INQ_VARID(FID,TRIM(VARNAME),VARID)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'error getting variable id: ',TRIM(VARNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    IERR = NF90_PUT_VAR(FID,VARID,TEMP,DIMS)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'4 error getting variable: ',TRIM(VARNAME)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    RETURN
  END SUBROUTINE PUTSVAR
  !==============================================================================!

END MODULE MOD_NCD
