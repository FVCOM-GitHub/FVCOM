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

module mod_utils
  USE MOD_PREC
  use CONTROL, only: IPT,TPI, DEG2RAD
  implicit none
  save
  
!
!--Communication Parameters FOR MPI IO MODE
!
   INTEGER, PARAMETER ::        SYNC_TAG  =3601
   INTEGER, PARAMETER ::        EXT_CODE  =3602
   INTEGER, PARAMETER ::        WAIT_CODE =3603
   INTEGER ::                   RESTART_CODE
   INTEGER ::                   NC_CODE 
   INTEGER ::                   NCSF_CODE 
   INTEGER ::                   NCAV_CODE 
   INTEGER ::                   INIT_CODE
   INTEGER ::                   NESTING_CODE 
   INTEGER ::                   INITNEST_CODE 

#  if defined (WAVE_CURRENT_INTERACTION)
   INTEGER ::                   NESTING_CODE_WAVE 
   INTEGER ::                   INITNEST_CODE_WAVE 
#  endif

!
!--Enumerate debugging levels
!
  integer,parameter::dbg_nbr=9 ! [nbr] Number of different debugging levels

  integer,parameter::dbg_log=0    ! [enm] Production mode. Debugging is
  integer,parameter::dbg_io=1     ! [enm] IO Filenames
  integer,parameter::dbg_scl=2    ! [enm] Scalars
  integer,parameter::dbg_mpi=3    ! [enm] Mpi Communication 
  integer,parameter::dbg_sbr=4    ! [enm] Subroutine names on entry and exit
  integer,parameter::dbg_sbrio=5  ! [enm] Subroutine I/O
  integer,parameter::dbg_vec=6    ! [enm] Entire vectors
  integer,parameter::dbg_vrb=7    ! [enm] Everything

  CHARACTER(LEN=80) :: prg_nm ! used in mod_sng and mod_input
  integer :: dbg_lvl ! [enm] Debugging level, initialized in input
  logical :: dbg_par=.false.

  INTERFACE SHUTDOWN_CHECK
     MODULE PROCEDURE SHUTDOWN_CHECK_1D
     MODULE PROCEDURE SHUTDOWN_CHECK_2D
  END INTERFACE



  ! THE PROJ BASED SUBROUTINE HAVE A DUMMY DEFINITION WHICH EXISTS
  ! EVEN IF PROJ IS NOT DEFINED TO REDUCE THE NUMBER OF IF DEFS IN
  ! THE CODE
  INTERFACE DEGREES2METERS
     MODULE PROCEDURE DEGREES2METERS_SCL_FLT
     MODULE PROCEDURE DEGREES2METERS_VEC_FLT
     MODULE PROCEDURE DEGREES2METERS_ARR_FLT
     MODULE PROCEDURE DEGREES2METERS_SCL_DBL
     MODULE PROCEDURE DEGREES2METERS_VEC_DBL
     MODULE PROCEDURE DEGREES2METERS_ARR_DBL
  END INTERFACE

  INTERFACE METERS2DEGREES
     MODULE PROCEDURE METERS2DEGREES_SCL_FLT
     MODULE PROCEDURE METERS2DEGREES_VEC_FLT
     MODULE PROCEDURE METERS2DEGREES_ARR_FLT
     MODULE PROCEDURE METERS2DEGREES_SCL_DBL
     MODULE PROCEDURE METERS2DEGREES_VEC_DBL
     MODULE PROCEDURE METERS2DEGREES_ARR_DBL
  END INTERFACE


  INTERFACE INTERP_ANODAL
     MODULE PROCEDURE INTERP_ANODAL_2D_FLT
     MODULE PROCEDURE INTERP_ANODAL_3D_FLT
     MODULE PROCEDURE INTERP_ANODAL_2D_DBL
     MODULE PROCEDURE INTERP_ANODAL_3D_DBL
  END INTERFACE

  INTERFACE INTERP_PNODAL
     MODULE PROCEDURE INTERP_PNODAL_2D_FLT
     MODULE PROCEDURE INTERP_PNODAL_3D_FLT
     MODULE PROCEDURE INTERP_PNODAL_2D_DBL
     MODULE PROCEDURE INTERP_PNODAL_3D_DBL
  END INTERFACE

  INTERFACE INTERP_AZONAL
     MODULE PROCEDURE INTERP_AZONAL_2D_FLT
     MODULE PROCEDURE INTERP_AZONAL_3D_FLT
     MODULE PROCEDURE INTERP_AZONAL_2D_DBL
     MODULE PROCEDURE INTERP_AZONAL_3D_DBL
  END INTERFACE

  INTERFACE INTERP_PZONAL
     MODULE PROCEDURE INTERP_PZONAL_2D_FLT
     MODULE PROCEDURE INTERP_PZONAL_3D_FLT
     MODULE PROCEDURE INTERP_PZONAL_2D_DBL
     MODULE PROCEDURE INTERP_PZONAL_3D_DBL
  END INTERFACE

  ! ELEMENTID = FIND_ELEMENT_CONTAINING(XLOC,YLOC, GUESS)
  ! GUESS IS OPTIONAL ARGUMENT

CONTAINS

  SUBROUTINE INITIALIZE_CONTROL(NAME)
    USE LIMS
    USE CONTROL
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    
!==============================================================================!
!  FVCOM VERSION                                                               !
!==============================================================================!

   FVCOM_VERSION     = 'FVCOM_4.3'
   FVCOM_WEBSITE     = 'http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu'
   INSTITUTION       = 'School for Marine Science and Technology'

   ! Set the IO UNIT value to screen output for now:
   IPT = 6


!==============================================================================!
!   SETUP PARALLEL ENVIRONMENT                                                 !
!==============================================================================!
   ! DEFAULT SETTINGS
                 SERIAL = .TRUE. 
                    PAR = .FALSE. 
                    MSR = .TRUE.
         IN_MPI_IO_LOOP = .FALSE.
        USE_MPI_IO_MODE = .FALSE.
                   MYID = 1
                  MSRID = 1
                 NPROCS = 1
               PRG_NAME = NAME
            MX_NBR_ELEM = 0
           NPROCS_TOTAL => NPROCS
               IOPROCID =-1

         ZEROTIME%MUSOD = 0
           ZEROTIME%MJD = 0
# if defined(MULTIPROCESSOR)
        MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
# else
        MPI_FVCOM_GROUP = -1
# endif
  END SUBROUTINE INITIALIZE_CONTROL


  logical  Function  dbg_set(vrb)
    USE CONTROL, only: MSR
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: vrb
    
    ! ONLY RETURN TRUE IF LEVEL PASSED IS .LE. THE RUN TIME LEVEL 
    dbg_set =.false.
    if (vrb <= dbg_lvl) then
       if ( MSR .OR. dbg_par) dbg_set =.true.
    end if
       
    return
  End Function DBG_SET


  Subroutine dbg_init(IPT_BASE,outtofile)
    use control, only: INFOFILE, PRG_NAME, MSR
    use lims, only: myid
    implicit none
    integer, intent(in):: IPT_BASE
    logical, intent(in):: outtofile
    character(LEN=3) :: ch3
    character(len=100) :: debugname

    if (outtofile .AND. MSR) then
       WRITE(IPT,*)"========================================================================"
       WRITE(IPT,*)"=== All further standard output goes to the user specified log file ===="
       WRITE(IPT,*)"=== Any further standard error messages will still print to screen ====="
       WRITE(IPT,*)"=== LOG FILE NAME: "//trim(INFOFILE)
       WRITE(IPT,*)"========================================================================"
       IPT = IPT_BASE+1
       CALL FOPEN(IPT, TRIM(INFOFILE) ,"ofr")
       
    end if

    if (dbg_par .AND. .NOT. MSR) then       
       IPT = IPT_BASE + MYID
       write(ch3,'(i3.3)') myid
       debugname=trim(PRG_NAME)//"_DEBUG."&
            & // trim(adjustl(ch3)) // ".log"
 
       CALL FOPEN(IPT, TRIM(debugname) ,"ofr")
       
    end if

  End Subroutine dbg_init


  Subroutine Fatal_Error(ER1,ER2,ER3,ER4)
    USE CONTROL, only: PRG_NAME
    implicit none
    character(Len=*) :: ER1
    CHARACTER(LEN=*), OPTIONAL   :: ER2
    CHARACTER(LEN=*), OPTIONAL   :: ER3
    CHARACTER(LEN=*), OPTIONAL   :: ER4

    if(DBG_SET(dbg_log)) then
       write(IPT,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       write(IPT,*) TRIM(PRG_NAME)//" Fatal Error!"
       write(IPT,*) ER1
       IF(PRESENT(ER2)) WRITE(IPT,*) ER2
       IF(PRESENT(ER3)) WRITE(IPT,*) ER3
       IF(PRESENT(ER4)) WRITE(IPT,*) ER4
       write(IPT,*) "Stopping "//TRIM(PRG_NAME)
       write(IPT,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    end if
    call pstop
  End Subroutine Fatal_Error

  Subroutine Warning(ER1,ER2,ER3,ER4)
    USE CONTROL, only: PRG_NAME
    implicit none
    character(Len=*) :: ER1
    CHARACTER(LEN=*), OPTIONAL   :: ER2
    CHARACTER(LEN=*), OPTIONAL   :: ER3
    CHARACTER(LEN=*), OPTIONAL   :: ER4

    if(DBG_SET(dbg_log)) then
       write(IPT,*) "+++++++++++++++++++++++++++++++++++++++++++++++++"
       write(IPT,*) TRIM(PRG_NAME)//" WARNING!"
       write(IPT,*) ER1
       IF(PRESENT(ER2)) WRITE(IPT,*) ER2
       IF(PRESENT(ER3)) WRITE(IPT,*) ER3
       IF(PRESENT(ER4)) WRITE(IPT,*) ER4
       write(IPT,*) TRIM(PRG_NAME)//" CONTINUEING"
       write(IPT,*) "+++++++++++++++++++++++++++++++++++++++++++++++++"
    end if

  End Subroutine Warning

  !==============================================================================|
  SUBROUTINE PSTOP               
    !==============================================================================|
# if defined (MULTIPROCESSOR)
    USE LIMS, ONLY: IOPROCID
    use CONTROL, only: IN_MPI_IO_LOOP, MSR

    INTEGER :: IERR, RECV, Ecode

!    if(IN_MPI_IO_LOOP)then
!       RECV = IOPROCID -1
!       CALL MPI_SEND(EXT_CODE,1,MPI_INTEGER,RECV,SYNC_TAG &
!            &,MPI_COMM_WORLD,IERR)
!    end if

    ! MPI_ABORT SHOULD CLEAR ALL COMMUNICATION
    ! I RECOMMEND USING '-kill' in mpiexec to force PBS to kill jobs
    ! where one processor hits a mpi_abort or mpi_finalize.

    ecode = -1
    CALL MPI_ABORT(MPI_COMM_WORLD,ecode,IERR)

    CALL MPI_FINALIZE(IERR)
# endif
    STOP
  END SUBROUTINE PSTOP
  

  !==============================================================================|
  SUBROUTINE PSHUTDOWN               
    !==============================================================================|
# if defined (MULTIPROCESSOR)
    USE LIMS, ONLY: IOPROCID
    use CONTROL, only: IN_MPI_IO_LOOP, MSR
    INTEGER IERR, RECV

    if(IN_MPI_IO_LOOP .AND. MSR)then
       RECV = IOPROCID -1
       CALL MPI_SEND(EXT_CODE,1,MPI_INTEGER,RECV,SYNC_TAG &
            &,MPI_COMM_WORLD,IERR)
    end if

    CALL MPI_FINALIZE(IERR)
# endif
    STOP
  END SUBROUTINE PSHUTDOWN
  
!==============================================================================|
!  CHECK DEPTH ARRAY FOR NAN.  SHUTDOWN IF FOUND                               |
!==============================================================================|

  SUBROUTINE SHUTDOWN_CHECK_1D(VAR,MSG)

    !==============================================================================|
    USE CONTROL
    IMPLICIT NONE
    REAL(DP) :: SBUF,RBUF  
    INTEGER  :: IERR,I,J
    REAL(SP), DIMENSION(:), INTENT(IN) :: VAR
    CHARACTER(LEN=*), OPTIONAL :: MSG
    !==============================================================================|

    IF (DBG_SET(DBG_SBR)) THEN
       IF (PRESENT(MSG)) THEN
          WRITE(IPT,*) "START: SHUTDOWN CHECK: "//MSG
       ELSE
          WRITE(IPT,*) "START: SHUTDOWN CHECK: no msg" 
       END IF
    END IF

    !Collect Depth Average to Master Processor
    SBUF = SUM(DBLE(VAR(1:ubound(VAR,1))))
    RBUF = SBUF
#  if defined (MULTIPROCESSOR)
    IF(PAR)CALL MPI_ALLREDUCE(SBUF,RBUF,1,MPI_DP,MPI_SUM,MPI_FVCOM_GROUP,IERR)
#  endif

    !Halt FVCOM if Depth Average = NaN          
#  if !defined (GFORTRAN) && !defined (ABSOFT)
    IF(ISNAN(RBUF))THEN 
       IF (PRESENT(MSG)) THEN
          CALL FATAL_ERROR("SHUTDOWN_CHECK FOUND NON FINITE VALUE:",&
               & MSG )
       ELSE
          CALL FATAL_ERROR('NON FINITE VALUE (DEPTH?) FOUND',&
               & 'MODEL HAS BECOME UNSTABLE')
       END IF

    END IF
#  else
    IF(RBUF /= RBUF)THEN
       IF (PRESENT(MSG)) THEN
          CALL FATAL_ERROR("SHUTDOWN_CHECK FOUND NON FINITE VALUE:",&
               & MSG )
       ELSE
          CALL FATAL_ERROR('NON FINITE VALUE (DEPTH?) FOUND',&
               & 'MODEL HAS BECOME UNSTABLE')
       END IF

    END IF
#  endif


    IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SHUTDOWN CHECK"

    RETURN
  END SUBROUTINE SHUTDOWN_CHECK_1D
!==============================================================================|
  SUBROUTINE SHUTDOWN_CHECK_2D(VAR,MSG)

    !==============================================================================|
    USE CONTROL
    IMPLICIT NONE
    REAL(DP) :: SBUF,RBUF  
    INTEGER  :: IERR,I,J
    REAL(SP), DIMENSION(:,:),INTENT(IN) :: VAR
    CHARACTER(LEN=*), OPTIONAL :: MSG
    !==============================================================================|

    IF (DBG_SET(DBG_SBR)) THEN
       IF (PRESENT(MSG)) THEN
          WRITE(IPT,*) "START: SHUTDOWN CHECK: "//MSG
       ELSE
          WRITE(IPT,*) "START: SHUTDOWN CHECK: no msg" 
       END IF
    END IF

    !Collect Depth Average to Master Processor
    SBUF = SUM(SUM(DBLE(VAR(1:ubound(VAR,1),:)),1),1)
    RBUF = SBUF
#  if defined (MULTIPROCESSOR)
    IF(PAR)CALL MPI_ALLREDUCE(SBUF,RBUF,1,MPI_DP,MPI_SUM,MPI_FVCOM_GROUP,IERR)
#  endif

    !Halt FVCOM if Depth Average = NaN          
#  if !defined (GFORTRAN) && !defined (ABSOFT)
    IF(ISNAN(RBUF))THEN 
       IF (PRESENT(MSG)) THEN
          CALL FATAL_ERROR("SHUTDOWN_CHECK FOUND NON FINITE VALUE:",&
               & MSG )
       ELSE
          CALL FATAL_ERROR('NON FINITE VALUE (DEPTH?) FOUND',&
               & 'MODEL HAS BECOME UNSTABLE')
       END IF

    END IF
#  else
    IF(RBUF /= RBUF)THEN
       IF (PRESENT(MSG)) THEN
          CALL FATAL_ERROR("SHUTDOWN_CHECK FOUND NON FINITE VALUE:",&
               & MSG )
       ELSE
          CALL FATAL_ERROR('NON FINITE VALUE (DEPTH?) FOUND',&
               & 'MODEL HAS BECOME UNSTABLE')
       END IF

    END IF
#  endif


    IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SHUTDOWN CHECK"

    RETURN
  END SUBROUTINE SHUTDOWN_CHECK_2D

!
! -- PROJECTION TOOL BOX STUFF
!
# if defined (PROJ)
  FUNCTION HAVE_PROJ(proj_ref) RESULT(RES)
    USE PROJ4
    IMPLICIT NONE
    LOGICAL :: RES
    character(LEN=*) :: proj_ref
    integer :: status, IDX, I
    type(prj90_projection) :: proj
    character(len=40), pointer :: PROJ4_PARAMS(:)

    RES = .FALSE.
  
    PROJ4_PARAMS => SET_PROJ_PARAMS(PROJ_REF)
    IF(.not. ASSOCIATED(PROJ4_PARAMS)) THEN
       CALL FATAL_ERROR("ERROR IN SET_PROJ_PARAMS: RESULT NOT ASSOCIAT&
            &ED")
    ELSE
       IDX = size(PROJ4_PARAMS)
       IF (IDX .LT. 1) CALL FATAL_ERROR&
            ("ERROR IN SET_PROJ_PARAMS: RESULT IS LENGTH 0")
    END IF

    IF (PROJ4_PARAMS(1) == 'none') RETURN
         
    if (dbg_set(dbg_log)) Then
       write(ipt,*) "! =========== Set Projection Parameters ============="
       DO I = 1,IDX
          write(ipt,*) "! "//TRIM(PROJ4_PARAMS(I))
       END DO
    write(ipt,*) "! ====================================================="
    end if


    status=prj90_init(proj,PROJ4_PARAMS)
    if (status.ne.PRJ90_NOERR)then
       CALL WARNING("DEGREES2METERS: prj90 returned an error;",&
            & prj90_strerrno(status))
    else
       RES=.true.
    end if
       
    status = prj90_free(proj)
    if (status.ne.PRJ90_NOERR)&
         & CALL FATAL_ERROR("DEGREES2METERS: prj90 returned an error which could not recover::",&
         & prj90_strerrno(status))

    DEALLOCATE(PROJ4_PARAMS)

    return
  END FUNCTION HAVE_PROJ

  FUNCTION set_proj_params(LINE_IN) result(STRINGVAL)
    IMPLICIT NONE
    Character(len=*):: LINE_IN
    Character(len=300):: TEXT_LINE
    Character(len=40), pointer :: STRINGVAL(:)
    Character(len=40) :: Tstring
    
    LOGICAL NEED_UNITS
    INTEGER EQCNT
    INTEGER LENGTH,EQLOC,blank1,blank2,CNT
    INTEGER I,unitsloc, IDX

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START set_proj_params"

    TEXT_LINE=TRIM(LINE_IN)

    LENGTH = LEN_TRIM(TEXT_LINE)

    EQCNT=0
    CNT =0

    EQLOC=1
    blank1=1
    blank2=1


    ! CHECK FOR ESCAPE VALUES
    IF('none' .EQ. TRIM(TEXT_LINE(1:4)) .OR. &
         & 'None' .EQ. TRIM(TEXT_LINE(1:4)) .OR.&
         & 'NONE' .EQ. TRIM(TEXT_LINE(1:4))) THEN

       allocate(STRINGVAL(1))
       STRINGVAL='none'
       RETURN
    END IF
       
!
!-----------------------CHANGE COMMAS TO BLANKS--------------------------------!
!
    DO I=1,LENGTH
       IF(TEXT_LINE(I:I) == ",") TEXT_LINE(I:I) = " "
    END DO
!
!-----------------------CHANGE 'PLUS' TO BLANKS--------------------------------!
!
    DO I=1,LENGTH
       IF(TEXT_LINE(I:I) == "+") TEXT_LINE(I:I) = " "
    END DO

!
!-------------- COUNT THE NUMBER OF PARAMETERS!
!

    DO
       IDX = INDEX(TEXT_LINE(EQLOC:LENGTH),"=")
       IF (IDX==0) THEN
          exit
       ELSE
          EQLOC=EQLOC+IDX
          EQCNT=EQCNT+1
       END IF
       
    END DO


   IF(EQCNT == 0) then
      CALL FATAL_ERROR(& 
        &'Could not find correct parameters for PROJ4.',&
        &'The null value is "none" otherwise the PROJECTION_REFERENCE variable',&
        &'in the NML_GRID_COORDINATES name list must contain a valid PROJ4 string.')
      RETURN
   END IF


   ! LOOK FOR UNITS STRING, ADD IF NEEDED
   UNITSLOC= INDEX(TEXT_LINE,"units=m")
   IF(UNITSLOC==0) THEN
      ALLOCATE(STRINGVAL(EQCNT+1))
      NEED_UNITS = .TRUE.
   ELSE
      ALLOCATE(STRINGVAL(EQCNT))
      NEED_UNITS=.FALSE.
   END IF


   ! PARSE STRING TO MAKE PARAMETER STRING
    EQLOC=1
    DO
       IDX = INDEX(TEXT_LINE(blank1:LENGTH),"=")
       IF (IDX ==0) EXIT
       
       EQLOC=BLANK1+IDX

!       write(IPT,*)"TEXT_LINE='"//text_line(BLANK1:LENGTH)//'"'
!       WRITE(IPT,*) "EQLOC=",EQLOC,"; IDX=",IDX
       

       IDX = INDEX(TEXT_LINE(EQLOC:LENGTH)," ")
       IF (IDX ==0) THEN
          BLANK2 = LENGTH
       ELSE
          BLANK2 = EQLOC+IDX-2
       END IF

!       write(IPT,*)"TEXT_LINE='"//text_line(EQLOC:LENGTH)//'"'
!       WRITE(IPT,*) "BLANK1=",BLANK1,"; BLANK2=",BLANK2,"; IDX=",IDX

!       write(IPT,*)"RESULT_LINE='"//text_line(BLANK1:BLANK2)//'"'


       CNT= CNT+1      
       STRINGVAL(CNT)=TRIM(ADJUSTL(TEXT_LINE(BLANK1:BLANK2)))
       
       BLANK1 = BLANK2+2

!       WRITE(IPT,*) "==========================================="

    END DO
    

!    WRITE(IPT,*)"HERE4"
   
   IF (CNT /= EQCNT) CALL FATAL_ERROR &
           &("The number of parameters found does not match the number of '=' in the proj4 string?")
   
   IF (NEED_UNITS) STRINGVAL(EQCNT+1)="units=m"

    IF(DBG_SET(DBG_SBRIO)) THEN
       DO I=1,SIZE(STRINGVAL)
          write(ipt,*) STRINGVAL(I)
       END DO
    END IF

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END set_proj_params"

 END FUNCTION set_proj_params

!=================================================================
  SUBROUTINE DEGREES2METERS_SCL_FLT(LON,LAT,proj_ref,X,Y)
    implicit none
    REAL(SPA), INTENT(IN)  :: LON, LAT
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(out):: X, Y
    
    REAL(DP):: XD, YD
    REAL(DP):: LOND, LATD

    LATD=LAT
    LOND=LON

    CALL DEGREES2METERS_SCL_DBL(LOND,LATD,proj_ref,XD,YD)

    X = XD
    Y = YD

  END SUBROUTINE DEGREES2METERS_SCL_FLT
!=================================================================
  SUBROUTINE DEGREES2METERS_VEC_FLT(LON,LAT,proj_ref,X,Y,nsze)
    implicit none
   integer, intent(in) :: nsze
    REAL(SPA), INTENT(IN)  :: LON(nsze), LAT(nsze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(out):: X(nsze), Y(nsze)

    REAL(DP):: XD(nsze), YD(nsze)
    REAL(DP):: LOND(nsze), LATD(nsze)

    LATD=LAT
    LOND=LON

    CALL DEGREES2METERS_VEC_DBL(LOND,LATD,proj_ref,XD,YD,nsze)

    X = XD
    Y = YD

  END SUBROUTINE DEGREES2METERS_VEC_FLT
!=================================================================
  SUBROUTINE DEGREES2METERS_ARR_FLT(LON,LAT,proj_ref,X,Y,nsze,msze)
    implicit none
    integer, intent(in) :: nsze, msze
    REAL(SPA), INTENT(IN)  :: LON(nsze,msze), LAT(nsze,msze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(out):: X(nsze,msze), Y(nsze,msze)
    REAL(DP) :: XD(nsze*msze),YD(nsze*msze)
    REAL(DP) :: LATD(nsze*msze),LOND(nsze*msze)
    integer :: i,lb,ub

    DO i = 1,msze
       lb = (i-1)*nsze+1
       ub = i*nsze
       LOND(lb:ub)=LON(1:nsze,i)
       LATD(lb:ub)=LAT(1:nsze,i)
    END DO

    CALL DEGREES2METERS_VEC_DBL(LOND,LATD,proj_ref,XD,YD,nsze*msze)

    DO i = 1,msze
       lb = (i-1)*nsze+1
       ub = i*nsze
       X(1:nsze,i)=XD(lb:ub)
       Y(1:nsze,i)=YD(lb:ub)
    END DO

  END SUBROUTINE DEGREES2METERS_ARR_FLT
!=================================================================
  SUBROUTINE DEGREES2METERS_SCL_DBL(LON,LAT,proj_ref,X,Y)
    USE PROJ4
    implicit none
    REAL(DP), INTENT(IN)  :: LON, LAT
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(out):: X, Y
    character(len=40), pointer :: PROJ4_PARAMS(:)
    integer :: status 
    type(prj90_projection) :: proj


    PROJ4_PARAMS => SET_PROJ_PARAMS(PROJ_REF)

    status=prj90_init(proj,PROJ4_PARAMS)
    if (status.ne.PRJ90_NOERR)&
         &CALL FATAL_ERROR("DEGREES2METERS: prj90 returned an error;",&
         & prj90_strerrno(status))

    status = prj90_fwd(proj,LON,LAT,x,y)
    if (status.ne.PRJ90_NOERR)&
         &CALL FATAL_ERROR("DEGREES2METERS: prj90 returned an error;",&
         & prj90_strerrno(status))

  status = prj90_free(proj)

  DEALLOCATE(PROJ4_PARAMS)

  END SUBROUTINE DEGREES2METERS_SCL_DBL
!=================================================================
  SUBROUTINE DEGREES2METERS_VEC_DBL(LON,LAT,proj_ref,X,Y,nsze)
    USE PROJ4
    implicit none
    integer, intent(in) :: nsze
    REAL(DP), INTENT(IN)  :: LON(nsze), LAT(nsze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(out):: X(nsze), Y(nsze)
    character(len=40), pointer :: PROJ4_PARAMS(:)
    integer :: status,I
    type(prj90_projection) :: proj


    PROJ4_PARAMS => SET_PROJ_PARAMS(PROJ_REF)

    status=prj90_init(proj,PROJ4_PARAMS)
    if (status.ne.PRJ90_NOERR)&
         &CALL FATAL_ERROR("DEGREES2METERS: prj90 returned an error;",&
         & prj90_strerrno(status))

    status = prj90_fwd(proj,LON,LAT,x,y)
    if (status.ne.PRJ90_NOERR)&
         &CALL FATAL_ERROR("DEGREES2METERS: prj90 returned an error;",&
         & prj90_strerrno(status))
   
  status = prj90_free(proj)

  DEALLOCATE(PROJ4_PARAMS)
  
  END SUBROUTINE DEGREES2METERS_VEC_DBL
!=================================================================
  SUBROUTINE DEGREES2METERS_ARR_DBL(LON,LAT,proj_ref,X,Y,nsze,msze)
    implicit none
    integer, intent(in) :: nsze, msze
    REAL(DP), INTENT(IN)  :: LON(nsze,msze), LAT(nsze,msze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(out):: X(nsze,msze), Y(nsze,msze)
    REAL(DP) :: Xvec(nsze*msze),Yvec(nsze*msze)
    REAL(DP) :: LATvec(nsze*msze),LONvec(nsze*msze)
    integer :: i,lb,ub

    DO i = 1,msze
       lb = (i-1)*nsze+1
       ub = i*nsze
       LONvec(lb:ub)=LON(1:nsze,i)
       LATvec(lb:ub)=LAT(1:nsze,i)
    END DO

    CALL DEGREES2METERS_VEC_DBL(LONvec,LATvec,proj_ref,Xvec,Yvec,nsze*msze)

    DO i = 1,msze
       lb = (i-1)*nsze+1
       ub = i*nsze
       X(1:nsze,i)=Xvec(lb:ub)
       Y(1:nsze,i)=Yvec(lb:ub)
    END DO


  END SUBROUTINE DEGREES2METERS_ARR_DBL
!=================================================================

!/////////////////////////////////////////////////////////////////

!=================================================================
  SUBROUTINE METERS2DEGREES_SCL_FLT(X,Y,proj_ref,LON,LAT)
    implicit none
    REAL(SPA), INTENT(OUT)  :: LON, LAT
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(IN):: X, Y
    
    REAL(DP):: XD, YD
    REAL(DP):: LOND, LATD


    XD = X
    YD = Y

    CALL METERS2DEGREES_SCL_DBL(XD,YD,proj_ref,LOND,LATD)

    LAT=LATD
    LON=LOND

  END SUBROUTINE METERS2DEGREES_SCL_FLT
!=================================================================
  SUBROUTINE METERS2DEGREES_VEC_FLT(X,Y,proj_ref,LON,LAT,nsze)
    implicit none
   integer, intent(in) :: nsze
    REAL(SPA), INTENT(OUT)  :: LON(nsze), LAT(nsze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(IN):: X(nsze), Y(nsze)

    REAL(DP):: XD(nsze), YD(nsze)
    REAL(DP):: LOND(nsze), LATD(nsze)

    XD = X
    YD = Y

    CALL METERS2DEGREES_VEC_DBL(XD,YD,proj_ref,LOND,LATD,nsze)

    LAT=LATD
    LON=LOND


  END SUBROUTINE METERS2DEGREES_VEC_FLT
!=================================================================
  SUBROUTINE METERS2DEGREES_ARR_FLT(X,Y,proj_ref,LON,LAT,nsze,msze)
    implicit none
    integer, intent(in) :: nsze, msze
    REAL(SPA), INTENT(OUT)  :: LON(nsze,msze), LAT(nsze,msze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(IN):: X(nsze,msze), Y(nsze,msze)
    REAL(DP) :: XD(nsze*msze),YD(nsze*msze)
    REAL(DP) :: LATD(nsze*msze),LOND(nsze*msze)
    integer :: i,lb,ub


    DO i = 1,msze
       lb = (i-1)*nsze+1
       ub = i*nsze
       XD(lb:ub)=X(1:nsze,i)
       YD(lb:ub)=Y(1:nsze,i)
    END DO

    CALL METERS2DEGREES_VEC_DBL(XD,YD,proj_ref,LOND,LATD,nsze*msze)

    DO i = 1,msze
       lb = (i-1)*nsze+1
       ub = i*nsze
       LON(1:nsze,i)=LOND(lb:ub)
       LAT(1:nsze,i)=LATD(lb:ub)
    END DO


  END SUBROUTINE METERS2DEGREES_ARR_FLT
!=================================================================
  SUBROUTINE METERS2DEGREES_SCL_DBL(X,Y,proj_ref,LON,LAT)
    USE PROJ4
    implicit none
    REAL(DP), INTENT(OUT)  :: LON, LAT
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(IN):: X, Y
    character(len=40), pointer :: PROJ4_PARAMS(:)
    integer :: status 
    type(prj90_projection) :: proj

    PROJ4_PARAMS => SET_PROJ_PARAMS(PROJ_REF)

    status=prj90_init(proj,PROJ4_PARAMS)
    if (status.ne.PRJ90_NOERR)&
         &CALL FATAL_ERROR("METERS2DEGREES: prj90 returned an error;",&
         & prj90_strerrno(status))

    status = prj90_inv(proj,x,y,LON,LAT)
    if (status.ne.PRJ90_NOERR)&
         &CALL FATAL_ERROR("METERS2DEGREES: prj90 returned an error;",&
         & prj90_strerrno(status))

  status = prj90_free(proj)
  
  DEALLOCATE(PROJ4_PARAMS)

  END SUBROUTINE METERS2DEGREES_SCL_DBL
!=================================================================
  SUBROUTINE METERS2DEGREES_VEC_DBL(X,Y,proj_ref,LON,LAT,nsze)
    USE PROJ4
    implicit none
    integer, intent(in) :: nsze
    REAL(DP), INTENT(OUT)  :: LON(nsze), LAT(nsze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(IN):: X(nsze), Y(nsze)
    character(len=40), pointer :: PROJ4_PARAMS(:)
    integer :: status
    type(prj90_projection) :: proj

    PROJ4_PARAMS => SET_PROJ_PARAMS(PROJ_REF)

    status=prj90_init(proj,PROJ4_PARAMS)
    if (status.ne.PRJ90_NOERR)&
         &CALL FATAL_ERROR("METERS2DEGREES: prj90 returned an error;",&
         & prj90_strerrno(status))

    status = prj90_inv(proj,x,y,LON,LAT)
    if (status.ne.PRJ90_NOERR)&
         &CALL FATAL_ERROR("METERS2DEGREES: prj90 returned an error;",&
         & prj90_strerrno(status))
   
  status = prj90_free(proj)

  DEALLOCATE(PROJ4_PARAMS)

  END SUBROUTINE METERS2DEGREES_VEC_DBL
!=================================================================
  SUBROUTINE METERS2DEGREES_ARR_DBL(X,Y,proj_ref,LON,LAT,nsze,msze)
    implicit none
    integer, intent(in) :: nsze, msze
    REAL(DP), INTENT(OUT)  :: LON(nsze,msze), LAT(nsze,msze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(IN):: X(nsze,msze), Y(nsze,msze)
    REAL(DP) :: Xvec(nsze*msze),Yvec(nsze*msze)
    REAL(DP) :: LATvec(nsze*msze),LONvec(nsze*msze)
    integer :: i,lb,ub

    DO i = 1,msze
       lb = (i-1)*nsze+1
       ub = i*nsze
       XVEC(lb:ub)=X(1:nsze,i)
       YVEC(lb:ub)=Y(1:nsze,i)
    END DO

    CALL METERS2DEGREES_VEC_DBL(Xvec,Yvec,proj_ref,LONvec,LATvec,nsze*msze)


    DO i = 1,msze
       lb = (i-1)*nsze+1
       ub = i*nsze
       LON(1:nsze,i)=LONvec(lb:ub)
       LAT(1:nsze,i)=LATvec(lb:ub)
    END DO



  END SUBROUTINE METERS2DEGREES_ARR_DBL
# else
! IF PROJ IS NOT DEFINED!
!%%%%%%%%%%%%%%%%%%%%%%%% DUMMY ROUTINES TO AVOID IF DEFS! %%%%%%%%%%%%%%%%%%%%
  FUNCTION HAVE_PROJ(proj_ref) RESULT(RES)
     LOGICAL :: RES
     character(LEN=*), INTENT(IN)    :: proj_ref
     RES=.FALSE.
   END FUNCTION HAVE_PROJ
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE DEGREES2METERS_SCL_FLT(LON,LAT,proj_ref,X,Y)
    implicit none
    REAL(SPA), INTENT(IN)  :: LON, LAT
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(out):: X, Y
    X=-1.0_SPA
    Y=-1.0_SPA
  END SUBROUTINE DEGREES2METERS_SCL_FLT
!=================================================================
  SUBROUTINE DEGREES2METERS_VEC_FLT(LON,LAT,proj_ref,X,Y,nsze)
    implicit none
   integer, intent(in) :: nsze
    REAL(SPA), INTENT(IN)  :: LON(nsze), LAT(nsze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(out):: X(nsze), Y(nsze)

    X=-1.0_SPA
    Y=-1.0_SPA
  END SUBROUTINE DEGREES2METERS_VEC_FLT
!=================================================================
  SUBROUTINE DEGREES2METERS_ARR_FLT(LON,LAT,proj_ref,X,Y,nsze,msze)
    implicit none
    integer, intent(in) :: nsze, msze
    REAL(SPA), INTENT(IN)  :: LON(nsze,msze), LAT(nsze,msze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(out):: X(nsze,msze), Y(nsze,msze)

    X=-1.0_SPA
    Y=-1.0_SPA
  END SUBROUTINE DEGREES2METERS_ARR_FLT
!=================================================================
  SUBROUTINE DEGREES2METERS_SCL_DBL(LON,LAT,proj_ref,X,Y)
!    USE PROJ4
    implicit none
    REAL(DP), INTENT(IN)  :: LON, LAT
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(out):: X, Y
    
    X=-1.0_DP
    Y=-1.0_DP
  END SUBROUTINE DEGREES2METERS_SCL_DBL
!=================================================================
  SUBROUTINE DEGREES2METERS_VEC_DBL(LON,LAT,proj_ref,X,Y,nsze)
!    USE PROJ4
    implicit none
    integer, intent(in) :: nsze
    REAL(DP), INTENT(IN)  :: LON(nsze), LAT(nsze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), INTENT(OUT) :: X(nsze), Y(nsze)
    X=-1.0_DP
    Y=-1.0_DP
  END SUBROUTINE DEGREES2METERS_VEC_DBL
!=================================================================
  SUBROUTINE DEGREES2METERS_ARR_DBL(LON,LAT,proj_ref,X,Y,nsze,msze)
    implicit none
    integer, intent(in) :: nsze, msze
    REAL(DP), INTENT(IN)  :: LON(nsze,msze), LAT(nsze,msze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), INTENT(OUT) :: X(nsze,msze), Y(nsze,msze)
    X=-1.0_DP
    Y=-1.0_DP
  END SUBROUTINE DEGREES2METERS_ARR_DBL
!=================================================================

!/////////////////////////////////////////////////////////////////

!=================================================================
  SUBROUTINE METERS2DEGREES_SCL_FLT(X,Y,proj_ref,LON,LAT)
    implicit none
    REAL(SPA), INTENT(OUT)  :: LON, LAT
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(IN):: X, Y
    
    LON=-1.0_SPA
    LAT=-1.0_SPA
  END SUBROUTINE METERS2DEGREES_SCL_FLT
!=================================================================
  SUBROUTINE METERS2DEGREES_VEC_FLT(X,Y,proj_ref,LON,LAT,nsze)
    implicit none
   integer, intent(in) :: nsze
    REAL(SPA), INTENT(OUT)  :: LON(nsze), LAT(nsze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(IN):: X(nsze), Y(nsze)

    LON=-1.0_SPA
    LAT=-1.0_SPA
  END SUBROUTINE METERS2DEGREES_VEC_FLT
!=================================================================
  SUBROUTINE METERS2DEGREES_ARR_FLT(X,Y,proj_ref,LON,LAT,nsze,msze)
    implicit none
    integer, intent(in) :: nsze, msze
    REAL(SPA), INTENT(OUT)  :: LON(nsze,msze), LAT(nsze,msze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(SPA), intent(IN):: X(nsze,msze), Y(nsze,msze)
        
    LON=-1.0_SPA
    LAT=-1.0_SPA
  END SUBROUTINE METERS2DEGREES_ARR_FLT
!=================================================================
  SUBROUTINE METERS2DEGREES_SCL_DBL(X,Y,proj_ref,LON,LAT)
!    USE PROJ4
    implicit none
    REAL(DP), INTENT(OUT)  :: LON, LAT
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(IN):: X, Y

    LON=-1.0_DP
    LAT=-1.0_DP
  END SUBROUTINE METERS2DEGREES_SCL_DBL
!=================================================================
  SUBROUTINE METERS2DEGREES_VEC_DBL(X,Y,proj_ref,LON,LAT,nsze)
!    USE PROJ4
    implicit none
    integer, intent(in) :: nsze
    REAL(DP), INTENT(OUT)  :: LON(nsze), LAT(nsze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(IN):: X(nsze), Y(nsze)

    LON=-1.0_DP
    LAT=-1.0_DP
  END SUBROUTINE METERS2DEGREES_VEC_DBL
!=================================================================
  SUBROUTINE METERS2DEGREES_ARR_DBL(X,Y,proj_ref,LON,LAT,nsze,msze)
    implicit none
    integer, intent(in) :: nsze, msze
    REAL(DP), INTENT(OUT)  :: LON(nsze,msze), LAT(nsze,msze)
    character(LEN=*), INTENT(IN)    :: proj_ref
    REAL(DP), intent(IN):: X(nsze,msze), Y(nsze,msze)

    LON=-1.0_DP
    LAT=-1.0_DP
  END SUBROUTINE METERS2DEGREES_ARR_DBL

# endif
! END OF PROJ4 SUBROUTINES
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
   FUNCTION INTERP_ANODAL_2D_FLT(xloc,yloc,i,Field) RESULT(FPT)
     IMPLICIT NONE
     REAL(SPA) :: FPT
     INTEGER, INTENT(IN) :: i ! Cell Number
     REAL(SPA), INTENT(IN):: xloc,yloc
     REAL(SPA), ALLOCATABLE, TARGET, INTENT(IN) ::FIELD(:) 
     REAL(SPA), POINTER :: PFIELD(:)
     
     PFIELD => FIELD
     FPT = INTERP_PNODAL_2D_FLT(xloc,yloc,i,PField)

   END FUNCTION INTERP_ANODAL_2D_FLT
!==============================================================================|
   FUNCTION INTERP_ANODAL_2D_DBL(xloc,yloc,i,Field) RESULT(FPT)
     IMPLICIT NONE
     REAL(DP) :: FPT
     INTEGER, INTENT(IN) :: i ! Cell Number
     REAL(DP), INTENT(IN):: xloc,yloc
     REAL(DP), ALLOCATABLE, TARGET, INTENT(IN) ::FIELD(:) 
     REAL(DP), POINTER :: PFIELD(:)
     
     PFIELD => FIELD
     FPT = INTERP_PNODAL_2D_DBL(xloc,yloc,i,PField)

   END FUNCTION INTERP_ANODAL_2D_DBL
!==============================================================================|
   FUNCTION INTERP_PNODAL_2D_FLT(xloc,yloc,i,Field) RESULT(FPT)
     USE ALL_VARS, only : aw0,awx,awy,nv,xc,yc
     IMPLICIT NONE
     REAL(SPA) :: FPT
     INTEGER, INTENT(IN) :: i  ! Cell Number
     REAL(SPA), INTENT(IN):: xloc,yloc
     REAL(SPA), POINTER, INTENT(IN) :: FIELD(:) 

     REAL(SPA):: X0c, Y0c,F0,Fx,Fy
     INTEGER :: n1,n2,n3
     REAL(SPA) :: dx_sph
     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PNODAL_2D"

     !element location (i) and surrounding nodes (n1,n2,n3)
     n1  = nv(i,1)
     n2  = nv(i,2)
     n3  = nv(i,3)

!offset from element center
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
     y0c = TPI * (Yloc - YC(I))
     dx_sph = Xloc - XC(I)
     if(dx_sph > 180.0_SPA) then
        dx_sph = dx_sph -360.0_SPA
     elseif(dx_sph < -180.0_SPA) then
        dx_sph = dx_sph + 360.0_SPA
     end if

     x0c = TPI * COS(DEG2RAD*(Yloc + YC(I))*0.5_SPA) * dx_sph
# else

     x0c = xloc - xc(i)
     y0c = yloc - yc(i)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

     !linear interpolation of Field
     F0 = aw0(i,1)*Field(n1)+aw0(i,2)*Field(n2)+aw0(i,3)*Field(n3)
     Fx = awx(i,1)*Field(n1)+awx(i,2)*Field(n2)+awx(i,3)*Field(n3)
     Fy = awy(i,1)*Field(n1)+awy(i,2)*Field(n2)+awy(i,3)*Field(n3)
     FPT = F0 + Fx*x0c + Fy*y0c


     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"END: INTERP_PNODAL_2D"
   END FUNCTION INTERP_PNODAL_2D_FLT
!==============================================================================|
!==============================================================================|
   FUNCTION INTERP_PNODAL_2D_DBL(xloc,yloc,i,Field) RESULT(FPT)
     USE ALL_VARS, only : aw0,awx,awy,nv,xc,yc
     IMPLICIT NONE
     REAL(DP) :: FPT
     INTEGER, INTENT(IN) :: i  ! Cell Number
     REAL(DP), INTENT(IN):: xloc,yloc
     REAL(DP), POINTER, INTENT(IN) :: FIELD(:) 

     REAL(DP):: X0c, Y0c,F0,Fx,Fy
     INTEGER :: n1,n2,n3
     REAL(DP) :: dx_sph
     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PNODAL_2D"

     !element location (i) and surrounding nodes (n1,n2,n3)
     n1  = nv(i,1)
     n2  = nv(i,2)
     n3  = nv(i,3)

!offset from element center
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------

     y0c = TPI * (Yloc - YC(I))
     dx_sph = Xloc - XC(I)
     if(dx_sph > 180.0_DP) then
        dx_sph = dx_sph -360.0_DP
     elseif(dx_sph < -180.0_DP) then
        dx_sph = dx_sph + 360.0_DP
     end if

     x0c = TPI * DCOS(DEG2RAD*(Yloc + YC(I))*0.5_DP) * dx_sph

# else
     !offset from element center
     x0c = xloc - xc(i)
     y0c = yloc - yc(i)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

     !linear interpolation of Field
     F0 = aw0(i,1)*Field(n1)+aw0(i,2)*Field(n2)+aw0(i,3)*Field(n3)
     Fx = awx(i,1)*Field(n1)+awx(i,2)*Field(n2)+awx(i,3)*Field(n3)
     Fy = awy(i,1)*Field(n1)+awy(i,2)*Field(n2)+awy(i,3)*Field(n3)
     FPT = F0 + Fx*x0c + Fy*y0c


     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"END: INTERP_PNODAL_2D"
   END FUNCTION INTERP_PNODAL_2D_DBL
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
   FUNCTION INTERP_ANODAL_3D_FLT(xloc,yloc,sigloc,lvls,i,Field) RESULT(FPT)
     IMPLICIT NONE
     REAL(SPA) :: FPT
     INTEGER, INTENT(IN) :: i,lvls ! Cell Number, Number of Levels in data (kb or kbm1)
     REAL(SPA), INTENT(IN):: xloc,yloc,sigloc
     REAL(SPA), ALLOCATABLE, TARGET, INTENT(IN) ::FIELD(:,:)
     REAL(SPA), POINTER :: PFIELD(:,:)

     PFIELD => FIELD
     FPT = INTERP_PNODAL_3D_FLT(xloc,yloc,sigloc,lvls,i,PField)

   END FUNCTION INTERP_ANODAL_3D_FLT
!==============================================================================|
   FUNCTION INTERP_ANODAL_3D_DBL(xloc,yloc,sigloc,lvls,i,Field) RESULT(FPT)
     IMPLICIT NONE
     REAL(DP) :: FPT
     INTEGER, INTENT(IN) :: i,lvls ! Cell Number, Number of Levels in data (kb or kbm1)
     REAL(DP), INTENT(IN):: xloc,yloc,sigloc
     REAL(DP), ALLOCATABLE, TARGET, INTENT(IN) ::FIELD(:,:)
     REAL(DP), POINTER :: PFIELD(:,:)

     PFIELD => FIELD
     FPT = INTERP_PNODAL_3D_DBL(xloc,yloc,sigloc,lvls,i,PField)

   END FUNCTION INTERP_ANODAL_3D_DBL
!==============================================================================|
   FUNCTION INTERP_PNODAL_3D_FLT(xloc,yloc,sigloc,lvls,i,Field) RESULT(fpt)
     USE ALL_VARS, only : aw0,awx,awy,nv,xc,yc,kbm2,kbm1,kb,z1,zz1,dz1,dzz1
     IMPLICIT NONE
     REAL(SPA) :: FPT
     INTEGER, INTENT(IN) :: i, lvls ! Cell Number, Number of Levels in data (kb or kbm1)
     REAL(SPA), INTENT(IN):: xloc,yloc,sigloc
     REAL(SPA), POINTER, INTENT(IN) :: FIELD(:,:) 

     REAL(SPA):: X0c,Y0c,F0,Fx,Fy,F_LOWER,F_UPPER, alpha,dsig
     INTEGER :: n1,n2,n3,k1,k2,k
     REAL(SPA) :: dx_sph

     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PNODAL_3D"
     
     !element location (i) and surrounding nodes (n1,n2,n3)
     n1  = nv(i,1)
     n2  = nv(i,2)
     n3  = nv(i,3)

!offset from element center
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
     y0c = TPI * (Yloc - YC(I))
     dx_sph = Xloc - XC(I)
     if(dx_sph > 180.0_SPA) then
        dx_sph = dx_sph -360.0_SPA
     elseif(dx_sph < -180.0_SPA) then
        dx_sph = dx_sph + 360.0_SPA
     end if

     x0c = TPI * COS(DEG2RAD*(Yloc + YC(I))*0.5_SPA) * dx_sph
# else

     x0c = xloc - xc(i)
     y0c = yloc - yc(i)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

     ! Determine the layer in which the point resides
     IF(LVLS == KBM1) THEN

        !top
        if(sigloc >= zz1(i,1))then
           k1 = 1
           k2 = 1
           alpha = -1
        elseif(sigloc > zz1(i,KBM1)) then!intermediate 
           do k=1,kbm2
              if(sigloc  < zz1(i,k) .and. sigloc >= zz1(i,k+1) )then
                 k1 = k
                 k2 = k+1 
                 alpha = (zz1(i,k)-sigloc)/dzz1(i,k)
                 exit
              endif
           end do
        else
           ! TOP
           k1 = KBM1
           k2 = KBM1
           alpha = -1
        endif
     ELSE IF(LVLS == KB) THEN
        
        !top
        if(sigloc >= z1(i,1))then
           k1 = 1
           k2 = 1
           alpha = -1
        elseif(sigloc > z1(i,KB)) then !intermediate 
           do k=1,kbm1
              if(sigloc  < z1(i,k) .and. sigloc >= z1(i,k+1) )then
                 k1 = k
                 k2 = k+1 
                 alpha = (z1(i,k)-sigloc)/dz1(i,k)
              endif
           end do
        else
           !bottom
           k1 = KBM1
           k2 = KBM1
           alpha = -1 
        endif

     ELSE
        CALL FATAL_ERROR("INTERP_PNODAL_3D: Invalid number of levels passed",&
             & "(Must be equal to either KB or KBM1")
     END IF
        
     !linear interpolation of Field
     F0 = aw0(i,1)*Field(n1,k1)+aw0(i,2)*Field(n2,k1)+aw0(i,3)*Field(n3,k1)
     Fx = awx(i,1)*Field(n1,k1)+awx(i,2)*Field(n2,k1)+awx(i,3)*Field(n3,k1)
     Fy = awy(i,1)*Field(n1,k1)+awy(i,2)*Field(n2,k1)+awy(i,3)*Field(n3,k1)
     F_UPPER = F0 + Fx*x0c + Fy*y0c

     IF(K1 == K2) THEN
        FPT = F_UPPER

     ELSE
        
        F0 = aw0(i,1)*Field(n1,k2)+aw0(i,2)*Field(n2,k2)+aw0(i,3)*Field(n3,k2)
        Fx = awx(i,1)*Field(n1,k2)+awx(i,2)*Field(n2,k2)+awx(i,3)*Field(n3,k2)
        Fy = awy(i,1)*Field(n1,k2)+awy(i,2)*Field(n2,k2)+awy(i,3)*Field(n3,k2)
        F_LOWER = F0 + Fx*x0c + Fy*y0c

        FPT = (alpha)*F_lower + (1.0_SPA-alpha)*F_upper
     END IF


     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"END: INTERP_PNODAL_3D"
   END FUNCTION INTERP_PNODAL_3D_FLT
!==============================================================================|
   FUNCTION INTERP_PNODAL_3D_DBL(xloc,yloc,sigloc,lvls,i,Field) RESULT(fpt)
     USE ALL_VARS, only : aw0,awx,awy,nv,xc,yc,kbm2,kbm1,kb,z1,zz1,dz1,dzz1
     IMPLICIT NONE
     REAL(DP) :: FPT
     INTEGER, INTENT(IN) :: i, lvls ! Cell Number, Number of Levels in data (kb or kbm1)
     REAL(DP), INTENT(IN):: xloc,yloc,sigloc
     REAL(DP), POINTER, INTENT(IN) :: FIELD(:,:) 

     REAL(DP):: X0c,Y0c,F0,Fx,Fy,F_LOWER,F_UPPER, alpha,dsig
     INTEGER :: n1,n2,n3,k1,k2,k
     REAL(DP) :: dx_sph

     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PNODAL_3D"
     
     !element location (i) and surrounding nodes (n1,n2,n3)
     n1  = nv(i,1)
     n2  = nv(i,2)
     n3  = nv(i,3)

!offset from element center
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------

     y0c = TPI * (Yloc - YC(I))
     dx_sph = Xloc - XC(I)
     if(dx_sph > 180.0_DP) then
        dx_sph = dx_sph -360.0_DP
     elseif(dx_sph < -180.0_DP) then
        dx_sph = dx_sph + 360.0_DP
     end if
     
     x0c = TPI * DCOS(DEG2RAD*(Yloc + YC(I))*0.5_DP) * dx_sph

# else
     
     x0c = xloc - xc(i)
     y0c = yloc - yc(i)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------


     ! Determine the layer in which the point resides
     IF(LVLS == KBM1) THEN

        !top
        if(sigloc >= zz1(i,1))then
           k1 = 1
           k2 = 1
           alpha = -1
        elseif(sigloc > zz1(i,KBM1)) then!intermediate 
           do k=1,kbm2
              if(sigloc  < zz1(i,k) .and. sigloc >= zz1(i,k+1) )then
                 k1 = k
                 k2 = k+1 
                 alpha = (zz1(i,k)-sigloc)/dzz1(i,k)
                 exit
              endif
           end do
        else
           ! TOP
           k1 = KBM1
           k2 = KBM1
           alpha = -1
        endif
     ELSE IF(LVLS == KB) THEN
        
        !top
        if(sigloc >= z1(i,1))then
           k1 = 1
           k2 = 1
           alpha = -1
        elseif(sigloc > z1(i,KB)) then !intermediate 
           do k=1,kbm1
              if(sigloc  < z1(i,k) .and. sigloc >= z1(i,k+1) )then
                 k1 = k
                 k2 = k+1 
                 alpha = (z1(i,k)-sigloc)/dz1(i,k)
              endif
           end do
        else
           !bottom
           k1 = KBM1
           k2 = KBM1
           alpha = -1 
        endif

     ELSE
        CALL FATAL_ERROR("INTERP_PNODAL_3D: Invalid number of levels passed",&
             & "(Must be equal to either KB or KBM1")
     END IF
        
     !linear interpolation of Field
     F0 = aw0(i,1)*Field(n1,k1)+aw0(i,2)*Field(n2,k1)+aw0(i,3)*Field(n3,k1)
     Fx = awx(i,1)*Field(n1,k1)+awx(i,2)*Field(n2,k1)+awx(i,3)*Field(n3,k1)
     Fy = awy(i,1)*Field(n1,k1)+awy(i,2)*Field(n2,k1)+awy(i,3)*Field(n3,k1)
     F_UPPER = F0 + Fx*x0c + Fy*y0c

     IF(K1 == K2) THEN
        FPT = F_UPPER

     ELSE
        
        F0 = aw0(i,1)*Field(n1,k2)+aw0(i,2)*Field(n2,k2)+aw0(i,3)*Field(n3,k2)
        Fx = awx(i,1)*Field(n1,k2)+awx(i,2)*Field(n2,k2)+awx(i,3)*Field(n3,k2)
        Fy = awy(i,1)*Field(n1,k2)+awy(i,2)*Field(n2,k2)+awy(i,3)*Field(n3,k2)
        F_LOWER = F0 + Fx*x0c + Fy*y0c

        FPT = (alpha)*F_lower + (1.0_DP-alpha)*F_upper
     END IF


     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"END: INTERP_PNODAL_3D"
   END FUNCTION INTERP_PNODAL_3D_DBL
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   FUNCTION INTERP_AZONAL_2D_FLT(xloc,yloc,i,Field) RESULT(FPT)
     IMPLICIT NONE
     REAL(SPA) :: FPT
     INTEGER, INTENT(IN) :: i ! Cell Number
     REAL(SPA), INTENT(IN):: xloc,yloc
     REAL(SPA), ALLOCATABLE, TARGET, INTENT(IN) ::FIELD(:) 
     REAL(SPA), POINTER :: PFIELD(:)
     
     PFIELD => FIELD
     FPT = INTERP_PZONAL_2D_FLT(xloc,yloc,i,PField)

   END FUNCTION INTERP_AZONAL_2D_FLT
!==============================================================================|
   FUNCTION INTERP_AZONAL_2D_DBL(xloc,yloc,i,Field) RESULT(FPT)
     IMPLICIT NONE
     REAL(DP) :: FPT
     INTEGER, INTENT(IN) :: i ! Cell Number
     REAL(DP), INTENT(IN):: xloc,yloc
     REAL(DP), ALLOCATABLE, TARGET, INTENT(IN) ::FIELD(:) 
     REAL(DP), POINTER :: PFIELD(:)
     
     PFIELD => FIELD
     FPT = INTERP_PZONAL_2D_DBL(xloc,yloc,i,PField)

   END FUNCTION INTERP_AZONAL_2D_DBL
!==============================================================================|
!  NOTE: ZONAL INTERP REQUIRES FIELD TO BE ALLOCATED 0:NT
!==============================================================================|
   FUNCTION INTERP_PZONAL_2D_FLT(xloc,yloc,i,Field) RESULT(FPT)
     USE ALL_VARS, only : a1u,a2u,xc,yc,nbe
     IMPLICIT NONE
     REAL(SPA) :: FPT
     INTEGER, INTENT(IN) :: i  ! Cell Number
     REAL(SPA), INTENT(IN):: xloc,yloc
     REAL(SPA), POINTER, INTENT(IN) :: FIELD(:) 

     REAL(SPA):: X0c, Y0c,F0,Fx,Fy
     INTEGER :: e1,e2,e3
     REAL(SPA) :: dx_sph
     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PZONAL_2D"

!offset from element center
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
     y0c = TPI * (Yloc - YC(I))
     dx_sph = Xloc - XC(I)
     if(dx_sph > 180.0_SPA) then
        dx_sph = dx_sph -360.0_SPA
     elseif(dx_sph < -180.0_SPA) then
        dx_sph = dx_sph + 360.0_SPA
     end if

     x0c = TPI * COS(DEG2RAD*(Yloc + YC(I))*0.5_SPA) * dx_sph
# else

     x0c = xloc - xc(i)
     y0c = yloc - yc(i)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------
     !Surrounding Element IDs
     e1  = nbe(i,1)
     e2  = nbe(i,2)
     e3  = nbe(i,3)

     !interpolate Field to the location
     Fx = a1u(i,1)*Field(i)+a1u(i,2)*Field(e1)+a1u(i,3)*Field(e2)+a1u(i,4)*Field(e3)
     Fy = a2u(i,1)*Field(i)+a2u(i,2)*Field(e1)+a2u(i,3)*Field(e2)+a2u(i,4)*Field(e3)
     FPT = Field(i) + Fx*x0c + Fy*y0c
     
     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PZONAL_2D"

   END FUNCTION INTERP_PZONAL_2D_FLT
!==============================================================================|
   FUNCTION INTERP_PZONAL_2D_DBL(xloc,yloc,i,Field) RESULT(FPT)
     USE ALL_VARS, only : a1u,a2u,xc,yc,nbe
     IMPLICIT NONE
     REAL(DP) :: FPT
     INTEGER, INTENT(IN) :: i  ! Cell Number
     REAL(DP), INTENT(IN):: xloc,yloc
     REAL(DP), POINTER, INTENT(IN) :: FIELD(:) 

     REAL(DP):: X0c, Y0c,F0,Fx,Fy
     INTEGER :: e1,e2,e3
     REAL(DP) :: dx_sph
     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PZONAL_2D"

!offset from element center
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------

     y0c = TPI * (Yloc - YC(I))
     dx_sph = Xloc - XC(I)
     if(dx_sph > 180.0_DP) then
        dx_sph = dx_sph -360.0_DP
     elseif(dx_sph < -180.0_DP) then
        dx_sph = dx_sph + 360.0_DP
     end if

     x0c = TPI * DCOS(DEG2RAD*(Yloc + YC(I))*0.5_DP) * dx_sph

# else
     !offset from element center
     x0c = xloc - xc(i)
     y0c = yloc - yc(i)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

     !Surrounding Element IDs
     e1  = nbe(i,1)
     e2  = nbe(i,2)
     e3  = nbe(i,3)

     !interpolate Field to the location
     Fx = a1u(i,1)*Field(i)+a1u(i,2)*Field(e1)+a1u(i,3)*Field(e2)+a1u(i,4)*Field(e3)
     Fy = a2u(i,1)*Field(i)+a2u(i,2)*Field(e1)+a2u(i,3)*Field(e2)+a2u(i,4)*Field(e3)
     FPT = Field(i) + Fx*x0c + Fy*y0c
     
     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PZONAL_2D"

   END FUNCTION INTERP_PZONAL_2D_DBL
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   FUNCTION INTERP_AZONAL_3D_FLT(xloc,yloc,sigloc,lvls,i,Field) RESULT(FPT)
     IMPLICIT NONE
     REAL(SPA) :: FPT
     INTEGER, INTENT(IN) :: i, lvls ! Cell Number, Number of levels(kb,or kbm1)
     REAL(SPA), INTENT(IN):: xloc,yloc,sigloc
     REAL(SPA), ALLOCATABLE, TARGET, INTENT(IN) ::FIELD(:,:) 
     REAL(SPA), POINTER :: PFIELD(:,:)
     
     PFIELD => FIELD
     FPT = INTERP_PZONAL_3D_FLT(xloc,yloc,sigloc,lvls,i,PField)

   END FUNCTION INTERP_AZONAL_3D_FLT
!==============================================================================|
   FUNCTION INTERP_AZONAL_3D_DBL(xloc,yloc,sigloc,lvls,i,Field) RESULT(FPT)
     IMPLICIT NONE
     REAL(DP) :: FPT
     INTEGER, INTENT(IN) :: i, lvls ! Cell Number, Number of levels(kb,or kbm1)
     REAL(DP), INTENT(IN):: xloc,yloc,sigloc
     REAL(DP), ALLOCATABLE, TARGET, INTENT(IN) ::FIELD(:,:) 
     REAL(DP), POINTER :: PFIELD(:,:)
     
     PFIELD => FIELD
     FPT = INTERP_PZONAL_3D_DBL(xloc,yloc,sigloc,lvls,i,PField)

   END FUNCTION INTERP_AZONAL_3D_DBL
!==============================================================================|
!  NOTE: ZONAL INTERP REQUIRES FIELD TO BE ALLOCATED 0:NT
!==============================================================================|
   FUNCTION INTERP_PZONAL_3D_FLT(xloc,yloc,sigloc,lvls,i,Field) RESULT(FPT)
     USE ALL_VARS, only : a1u,a2u,xc,yc,nbe,kbm2,kbm1,kb,z1,dz1,zz1,dzz1
     IMPLICIT NONE
     REAL(SPA) :: FPT
     INTEGER, INTENT(IN) :: i,lvls  ! Cell Number, Number of levels(kb,or kbm1)
     REAL(SPA), INTENT(IN):: xloc,yloc,sigloc
     REAL(SPA), POINTER, INTENT(IN) :: FIELD(:,:) 

     REAL(SPA):: X0c, Y0c,F0,Fx,Fy,alpha,F_upper,F_lower
     INTEGER :: e1,e2,e3,k1,k2,k
     REAL(SPA) :: dx_sph
     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PZONAL_3D"

!offset from element center
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
     y0c = TPI * (Yloc - YC(I))
     dx_sph = Xloc - XC(I)
     if(dx_sph > 180.0_SPA) then
        dx_sph = dx_sph -360.0_SPA
     elseif(dx_sph < -180.0_SPA) then
        dx_sph = dx_sph + 360.0_SPA
     end if

     x0c = TPI * COS(DEG2RAD*(Yloc + YC(I))*0.5_spa) * (dx_sph)

# else
     !offset from element center
     x0c = xloc - xc(i)
     y0c = yloc - yc(i)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

     !Surrounding Element IDs
     e1  = nbe(i,1)
     e2  = nbe(i,2)
     e3  = nbe(i,3)

!----determine sigma layers above and below sigloc

     IF(LVLS == KBM1) THEN

        k1  = 1
        k2  = 1
        alpha = -1
        if(sigloc < zz1(i,kbm1)) then !!particle near bottom
           k1 = kbm1
           k2 = kbm1
           alpha = -1
        else
           do k=1,kbm2
              if(sigloc  < zz1(i,k) .and. sigloc >= zz1(i,k+1) )then
                 k1 = k
                 k2 = k+1 
                 alpha = (zz1(i,k)-sigloc)/dzz1(i,k)
                 exit
              endif
           end do
        end if
        
     ELSE IF(LVLS == KB) THEN
        
        !surface (default)
        k1 = 1
        k2 = 1
        alpha = -1        
        !bottom
        if(sigloc < z1(i,kb))then
           k1 = kb
           k2 = kb
           alpha = -1
        else !intermediate 
           do k=1,kbm1
              if(sigloc  < z1(i,k) .and. sigloc >= z1(i,k+1) )then
                 k1 = k
                 k2 = k+1 
                 alpha = (z1(i,k)-sigloc)/dz1(i,k)
              endif
           end do
        endif

     ELSE
        CALL FATAL_ERROR("INTERP_PZONAL_3D: Invalid number of levels passed",&
             & "(Must be equal to either KB or KBM1")
     END IF

     !interpolate Field to the location
     Fx = a1u(i,1)*Field(i,k1)+a1u(i,2)*Field(e1,k1)+a1u(i,3)*Field(e2,k1)+a1u(i,4)*Field(e3,k1)
     Fy = a2u(i,1)*Field(i,k1)+a2u(i,2)*Field(e1,k1)+a2u(i,3)*Field(e2,k1)+a2u(i,4)*Field(e3,k1)
     F_upper = Field(i,k1) + Fx*x0c + Fy*y0c
     
     IF(K1 == K2) THEN
        FPT = F_UPPER
     ELSE
        
        Fx = a1u(i,1)*Field(i,k2)+a1u(i,2)*Field(e1,k2)+a1u(i,3)*Field(e2,k2)+a1u(i,4)*Field(e3,k2)
        Fy = a2u(i,1)*Field(i,k2)+a2u(i,2)*Field(e1,k2)+a2u(i,3)*Field(e2,k2)+a2u(i,4)*Field(e3,k2)
        F_lower = Field(i,k2) + Fx*x0c + Fy*y0c
        
        FPT = (alpha)*F_lower + (1.0-alpha)*F_upper
     END IF


     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PZONAL_3D"

   END FUNCTION INTERP_PZONAL_3D_FLT
!==============================================================================|
   FUNCTION INTERP_PZONAL_3D_DBL(xloc,yloc,sigloc,lvls,i,Field) RESULT(FPT)
     USE ALL_VARS, only : a1u,a2u,xc,yc,nbe,kbm2,kbm1,kb,z1,dz1,zz1,dzz1
     IMPLICIT NONE
     REAL(DP) :: FPT
     INTEGER, INTENT(IN) :: i,lvls  ! Cell Number, Number of levels(kb,or kbm1)
     REAL(DP), INTENT(IN):: xloc,yloc,sigloc
     REAL(DP), POINTER, INTENT(IN) :: FIELD(:,:) 

     REAL(DP):: X0c, Y0c,F0,Fx,Fy,alpha,F_upper,F_lower
     INTEGER :: e1,e2,e3,k1,k2,k

     REAL(DP) :: dx_sph
     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PZONAL_3D"


!offset from element center
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------

     y0c = TPI * (Yloc - YC(I))
     dx_sph = Xloc - XC(I)
     if(dx_sph > 180.0_DP) then
        dx_sph = dx_sph -360.0_DP
     elseif(dx_sph < -180.0_DP) then
        dx_sph = dx_sph + 360.0_DP
     end if

     x0c = TPI * DCOS(DEG2RAD*(Yloc + YC(I))*0.5_DP) * dx_sph

# else
     !offset from element center
     x0c = xloc - xc(i)
     y0c = yloc - yc(i)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

     !Surrounding Element IDs
     e1  = nbe(i,1)
     e2  = nbe(i,2)
     e3  = nbe(i,3)

!----determine sigma layers above and below sigloc

     IF(LVLS == KBM1) THEN

        k1  = 1
        k2  = 1
        alpha = -1
        if(sigloc < zz1(i,kbm1)) then !!particle near bottom
           k1 = kbm1
           k2 = kbm1
           alpha = -1
        else
           do k=1,kbm2
              if(sigloc  < zz1(i,k) .and. sigloc >= zz1(i,k+1) )then
                 k1 = k
                 k2 = k+1 
                 alpha = (zz1(i,k)-sigloc)/dzz1(i,k)
                 exit
              endif
           end do
        end if
        
     ELSE IF(LVLS == KB) THEN
        
        !surface (default)
        k1 = 1
        k2 = 1
        alpha = -1        
        !bottom
        if(sigloc < z1(i,kb))then
           k1 = kb
           k2 = kb
           alpha = -1
        else !intermediate 
           do k=1,kbm1
              if(sigloc  < z1(i,k) .and. sigloc >= z1(i,k+1) )then
                 k1 = k
                 k2 = k+1 
                 alpha = (z1(i,k)-sigloc)/dz1(i,k)
              endif
           end do
        endif

     ELSE
        CALL FATAL_ERROR("INTERP_PZONAL_3D: Invalid number of levels passed",&
             & "(Must be equal to either KB or KBM1")
     END IF

     !interpolate Field to the location
     Fx = a1u(i,1)*Field(i,k1)+a1u(i,2)*Field(e1,k1)+a1u(i,3)*Field(e2,k1)+a1u(i,4)*Field(e3,k1)
     Fy = a2u(i,1)*Field(i,k1)+a2u(i,2)*Field(e1,k1)+a2u(i,3)*Field(e2,k1)+a2u(i,4)*Field(e3,k1)
     F_upper = Field(i,k1) + Fx*x0c + Fy*y0c
     
     IF(K1 == K2) THEN
        FPT = F_UPPER
     ELSE
        
        Fx = a1u(i,1)*Field(i,k2)+a1u(i,2)*Field(e1,k2)+a1u(i,3)*Field(e2,k2)+a1u(i,4)*Field(e3,k2)
        Fy = a2u(i,1)*Field(i,k2)+a2u(i,2)*Field(e1,k2)+a2u(i,3)*Field(e2,k2)+a2u(i,4)*Field(e3,k2)
        F_lower = Field(i,k2) + Fx*x0c + Fy*y0c
        
        FPT = (alpha)*F_lower + (1.0_DP-alpha)*F_upper
     END IF


     IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: INTERP_PZONAL_3D"

   END FUNCTION INTERP_PZONAL_3D_DBL
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   Function Find_Element_Containing(xloc,yloc,GUESS) RESULT(EID)
     !==============================================================================|
     !  find home element for points (x,y)                                          |
     !  search nearest element to progressively further elements.
     !==============================================================================|

     !------------------------------------------------------------------------------|

     use all_vars
     implicit none
     INTEGER :: EID

     !------------------------------------------------------------------------------|
     real(sp),INTENT(IN) :: xloc,yloc
     INTEGER, INTENT(IN),OPTIONAL :: GUESS

     IF(PRESENT(GUESS)) THEN
        EID  = FIND_ELEMENT_CONTAINING_quick(xloc,yloc,Guess)
        IF (EID /= 0) RETURN
     END IF

     EID  = FIND_ELEMENT_CONTAINING_robust(xloc,yloc)

   END Function Find_Element_Containing
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   Function Find_Element_Containing_robust(xloc,yloc) RESULT(EID)
     !==============================================================================|
     !  find home element for points (x,y)                                          |
     !  search nearest element to progressively further elements.
     !==============================================================================|

     !------------------------------------------------------------------------------|

     use all_vars
     USE MOD_SPHERICAL
     implicit none
     INTEGER :: EID

     !------------------------------------------------------------------------------|
     real(sp),INTENT(IN) :: xloc,yloc

     integer i,min_loc
     real(sp), dimension(1:nt,1) :: radlist
     real(sp), dimension(3) :: xtri,ytri
     real(sp) :: radlast
     integer  :: locij(2), cnt

     !==============================================================================|

     EID = 0

     cnt = 0
!     radlist(1:nt,1) = sqrt((xc(1:nt)-xloc)**2 + (yc(1:nt)-yloc)**2)

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
     radlist=10.0

     Do I=1,NT
        if(abs(yloc-yc(i))<1.0) &
        CALL ARC(xloc,yloc,XC(I),YC(I),radlist(i,1))
        if(radlist(i,1)<1.0E-6) goto 10
     END Do
10 continue
# else
     radlist(1:nt,1) = abs(xc(1:nt)-xloc) + abs(yc(1:nt)-yloc)
# endif

     radlast = -1.0_sp
     in:  do while(cnt < 50)
        cnt = cnt+1
        locij   = minloc(radlist,radlist>radlast)
        min_loc = locij(1)
        if(min_loc == 0) then
           exit in
        end if

        if(isintriangle(min_loc,xloc,yloc))then
           EID = min_loc
           exit in 
        end if
        radlast = radlist(min_loc,1)
     end do in

     return
   end Function Find_Element_Containing_robust
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   FUNCTION FIND_ELEMENT_CONTAINING_quick(xloc,yloc,Guess) RESULT(EID)
     !==============================================================================|
     !  determine which element a location reside in by searching neighboring 
     !  elements.  
     !==============================================================================|

     !------------------------------------------------------------------------------|

     use all_vars, only: nv,nbve,ntve,nv
     implicit none
     INTEGER :: EID

     REAL(SP), INTENT(IN) :: Xloc,Yloc
     INTEGER, INTENT(IN)  :: GUESS

     integer i,j,k,iney,ncheck
     real(sp), dimension(3) :: xlast,ylast,xney,yney

     !==============================================================================|
     EID = 0 
     IF (GUESS == 0) RETURN

     if(isintriangle(GUESS,xloc,yloc))then       !!particle remains in element
        EID = GUESS
     else                                             !!check neighbors
        outer: do j=1,3
           ncheck = nv(GUESS,j)
           do k=1,ntve(ncheck)
              iney = nbve(ncheck,k) 
              if(isintriangle(iney,xloc,yloc))then
                 EID   = iney 
                 exit outer
              end if
           end do
        end do outer
     end if
     
     return
   end FUNCTION FIND_ELEMENT_CONTAINING_quick
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   LOGICAL FUNCTION ISINTRI(X0,Y0,Xt,Yt)
     use mod_prec
     implicit none
     real(sp), intent(in) :: x0,Y0
     real(sp), intent(in) :: xt(3),yt(3)
     real(sp) :: f1,f2,f3

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------

     real(sp) :: lx0, ly0
     real(sp) :: lxt(3),lyt(3)


     isintri = .false. 
     
     ! MAKE SURE ALL PTS IN TRIANGLE ARE NOT SEPERATED BY MORE THAN 180
     LXT = XT
     WHERE(XT -XT(1) < -180.0_SP)
        LXT = XT + 360.0_SP
     ELSEWHERE(XT-XT(1) > 180.0_SP)
        LXT = XT - 360.0_SP
     END WHERE

     ! MAKE SURE THE PNT IS NOT MORE THAN 180 FROM THE PTS OF THE TRIANGLE
     LX0 = X0
     IF(ANY(LXT-X0 < -180.0_SP)) THEN
        LX0 = X0 - 360.0_SP
     ELSEIF(ANY(LXT-X0 > 180.0_SP)) THEN
        LX0 = X0 + 360.0_SP
     END IF

     LY0 = Y0
     LYT = YT

     ! GET LOCAL PROJECTED COORDINATES
!     ly0 = DEG2RAD*y0
!     lx0 = cos(ly0)*DEG2RAD*lx0
     
!     lyt = DEG2RAD*yt
!     lxt = cos(lyt)*DEG2RAD*lxt

!     WRITE(IPT,*) "====== ISINTRI==========="
!          write(ipt,*) "LXT",LXT
!          write(ipt,*) "LYT",LYT
!          write(ipt,*) "LX0,LY0",LX0,LY0



     if(ly0 < minval(lyt) .or. ly0 > maxval(lyt)) then
        isintri = .false.
        return
     endif
     if(lx0 < minval(lxt) .or. lx0 > maxval(lxt)) then
        isintri = .false.
        return
     endif
     
     f1 = (ly0-lyt(1))*(lxt(2)-lxt(1)) - (lx0-lxt(1))*(lyt(2)-lyt(1))
     f2 = (ly0-lyt(3))*(lxt(1)-lxt(3)) - (lx0-lxt(3))*(lyt(1)-lyt(3))
     f3 = (ly0-lyt(2))*(lxt(3)-lxt(2)) - (lx0-lxt(2))*(lyt(3)-lyt(2))
     if(f1*f3 >= 0.0_sp .and. f3*f2 >= 0.0_sp) isintri = .true.

!---------------------------------------------------------------
# else
!---------------------------------------------------------------

     isintri = .false. 
     
     if(y0 < minval(yt) .or. y0 > maxval(yt)) then
     isintri = .false.
     return
   endif
   if(x0 < minval(xt) .or. x0 > maxval(xt)) then
     isintri = .false.
     return
   endif

   f1 = (y0-yt(1))*(xt(2)-xt(1)) - (x0-xt(1))*(yt(2)-yt(1))
   f2 = (y0-yt(3))*(xt(1)-xt(3)) - (x0-xt(3))*(yt(1)-yt(3))
   f3 = (y0-yt(2))*(xt(3)-xt(2)) - (x0-xt(2))*(yt(3)-yt(2))
   if(f1*f3 >= 0.0_sp .and. f3*f2 >= 0.0_sp) isintri = .true.
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------


   return

 END FUNCTION ISINTRI


!==============================================================================|
   logical function isintriangle(i,x0,y0) 
!==============================================================================|
!  determine if point (x0,y0) is in triangle defined by nodes (xt(3),yt(3))    |
!  using algorithm used for scene rendering in computer graphics               |
!  algorithm works well unless particle happens to lie in a line parallel      |
!  to the edge of a triangle.                                                  |
!  This can cause problems if you use a regular grid, say for idealized        |
!  modelling and you happen to see particles right on edges or parallel to     |
!  edges.                                                                      |
!==============================================================================|

   use mod_prec
   use all_vars, only : nv, vx, vy
   implicit none
   real(sp), intent(in) :: x0,y0
   integer, intent(in)  :: i
   real(sp) :: xt(3),yt(3)
   real(sp) :: f1,f2,f3
   real(sp) :: x1(2)
   real(sp) :: x2(2)
   real(sp) :: x3(2)
   real(sp) :: p(2)

!------------------------------------------------------------------------------|

   isintriangle = .false. 

   xt = vx(nv(i,1:3))
   yt = vy(nv(i,1:3))

   isintriangle = ISINTRI(X0,Y0,Xt,Yt)

!
!   if(sameside(p,x1,x2,x3).and.sameside(p,x2,x1,x3).and. &
!      sameside(p,x3,x1,x2)) isintriangle = .true. 
!!$   if(y0 < minval(yt) .or. y0 > maxval(yt)) then
!!$     isintriangle = .false.
!!$     return
!!$   endif
!!$   if(x0 < minval(xt) .or. x0 > maxval(xt)) then
!!$     isintriangle = .false.
!!$     return
!!$   endif
!!$
!!$   f1 = (y0-yt(1))*(xt(2)-xt(1)) - (x0-xt(1))*(yt(2)-yt(1))
!!$   f2 = (y0-yt(3))*(xt(1)-xt(3)) - (x0-xt(3))*(yt(1)-yt(3))
!!$   f3 = (y0-yt(2))*(xt(3)-xt(2)) - (x0-xt(2))*(yt(3)-yt(2))
!!$   if(f1*f3 >= 0.0_sp .and. f3*f2 >= 0.0_sp) isintriangle = .true.

   return
 end function isintriangle
!==============================================================================|
  function sameside(p1,p2,a,b) result(value)
     real(sp), intent(in) :: p1(2)
     real(sp), intent(in) :: p2(2)
     real(sp), intent(in) :: a(2)
     real(sp), intent(in) :: b(2)
     logical value
     real(sp) :: cp1,cp2
  
     cp1 = (b(1)-a(1))*(p1(2)-a(2)) - (b(2)-a(2))*(p1(1)-a(1))
     cp2 = (b(1)-a(1))*(p2(2)-a(2)) - (b(2)-a(2))*(p2(1)-a(1))
  
     value = .false.
     if(cp1*cp2 >= 0) value = .true.

  end function sameside
!==============================================================================|

  ! SEQUENTIALLY FIND NEIGHBORING CELLS OR NODES WHERE DATA EXISTS
  ! AND MAKE A CREATE AN INDEX TO AVERAGE INTO CELL WHERE DATA IS MISSING

  ! ONLY WORKS FOR LOCAL DOMAIN DATA
  SUBROUTINE GRID_NEIGHBOR_INDEX(FOUND,IDEX,CNT,ORDER)
    USE ALL_VARS
    IMPLICIT NONE
    INTEGER, POINTER :: FOUND(:)
    INTEGER, POINTER :: IDEX(:,:)
    INTEGER, POINTER :: CNT(:)
    INTEGER, POINTER :: ORDER(:)
    
    INTEGER :: i,j
    INTEGER :: LOOP
    INTEGER :: ORD
    ! DO NO ALLOCATION HERE - ASSUME ALL VARIALBE ARE ALLOCATED
    
    ! LOOKING FOR NEIGHBORING NODES
    IF (ubound(FOUND,1) == MT) THEN
    
       LOOP = 0
       ORD = 0
       DO WHILE(ANY(FOUND==-1))
          
          LOOP = LOOP +1
          IF(LOOP>M) CALL FATAL_ERROR&
               &("LOOP COUNT EXCEEDED IN GRID_NEIGHBOR_INDEX")

          DO i=1,M
             
             IF(FOUND(i) == -1) THEN

                ! look to see if this node has neighbors which are set
                DO j=1,NTSN(i)
                   IF(FOUND(NBSN(i,j))>-1.and.FOUND(NBSN(i,j))<loop )THEN
                      ! INCREASE THE COUNT FOR AVERAGE
                      CNT(i)= CNT(i) +1
                      ! RECORD THE INDEX
                      IDEX(i,CNT(i)) = NBSN(i,j)
                      ! RECORD WHICH LOOP FOUND IT
                      FOUND(I) = LOOP

                   END IF
                   
                   
                END DO
                
                ! RECORD THE ORDER IN WHICH THIS NODE WAS FOUND
                IF (CNT(i) >0) THEN
                   ORD = ORD +1
                   ORDER(ORD) = i
                END IF
                   
             END IF

          END DO

       END DO

    ! LOOKING FOR NEIGHBORING ELEMENTS
    ELSEIF (ubound(FOUND,1) == NT) THEN

       LOOP = 0
       ORD = 0
       DO WHILE(ANY(FOUND==-1))
          
          LOOP = LOOP +1
          IF(LOOP>N) CALL FATAL_ERROR&
               &("LOOP COUNT EXCEEDED IN GRID_NEIGHBOR_INDEX")

          DO i=1,N
             
             IF(FOUND(i) == -1) THEN

                ! look to see if this element has neighbors which are found
                DO j=1,3
                   IF(FOUND(NBE(i,j))>-1 .and.FOUND(NBE(i,j))<loop )THEN
                      ! INCREASE THE COUNT FOR AVERAGE
                      CNT(i)= CNT(i) +1
                      ! RECORD THE INDEX
                      IDEX(i,CNT(i)) = NBE(i,j)
                      ! RECORD WHICH LOOP FOUND IT
                      FOUND(i) = LOOP

                   END IF
                   
                   
                END DO
                
                ! RECORD THE ORDER IN WHICH THIS NODE WAS FOUND
                IF (CNT(i) >0) THEN
                   ORD = ORD +1
                   ORDER(ORD) = i
                END IF
                   
             END IF

          END DO

       END DO


    ELSE

       CALL FATAL_ERROR("PASSED INVALID SIZE TO GRID_NEIGHBOR_INDEX ???")
    
    END IF
  END SUBROUTINE GRID_NEIGHBOR_INDEX

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
   SUBROUTINE WRITE_BANNER(PAR,NP,ID)
     IMPLICIT NONE
     LOGICAL, INTENT(IN) :: PAR
     INTEGER, INTENT(IN) :: NP,ID
     character(len=8) :: chr_np,chr_id

! CREATED USING: http://www.mudmagic.com/figlet-server/
!     146) rounded.flf :

WRITE(IPT,*)'!================================================================!'
WRITE(IPT,*)'  _______  _     _  _______  _______  _______  ______     _____  '
WRITE(IPT,*)' (_______)(_)   (_)(_______)(_______)(_______)(_____ \   (_____) '
WRITE(IPT,*)'  _____    _     _  _        _     _  _  _  _  _____) )  _  __ _ '
WRITE(IPT,*)' |  ___)  | |   | || |      | |   | || ||_|| |(_____ (  | |/ /| |'
WRITE(IPT,*)' | |       \ \ / / | |_____ | |___| || |   | | _____) )_|   /_| |'
WRITE(IPT,*)' |_|        \___/   \______) \_____/ |_|   |_|(______/(_)\_____/ '
WRITE(IPT,*)' -- Beta Release'
WRITE(IPT,*)'!================================================================!'
WRITE(IPT,*)'!                                                                !'
WRITE(IPT,*)'!========DOMAIN DECOMPOSITION USING: METIS 4.0.1 ================!'
WRITE(IPT,*)'!======Copyright 1998, Regents of University of Minnesota========!'
WRITE(IPT,*)'!                                                                !'
IF(PAR) THEN
   WRITE(chr_np,'(I3.3)') NP
   WRITE(chr_id,'(I3.3)') ID
WRITE(IPT,*)'!================================================================!'
WRITE(IPT,*)'!                                                                !'
WRITE(IPT,*)'! RUNNING IN PARALLEL: '//trim(chr_np)//' Processors                            !'
WRITE(IPT,*)'! MYID is '//trim(chr_id)//'                                                    !'
WRITE(IPT,*)'!================================================================!'
END IF
   RETURN
   END SUBROUTINE WRITE_BANNER
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==============================================================================|
SUBROUTINE FOPEN(IUNIT,INSTR,IOPT) 
!==============================================================================|
! FOPEN Utility to open non Netcdf files
!==============================================================================|
  
   IMPLICIT NONE
   INTEGER, INTENT(IN)             :: IUNIT
   CHARACTER(LEN=*)                :: INSTR 
   CHARACTER(LEN=3), INTENT(IN)    :: IOPT  
   CHARACTER(LEN=11) :: FORMSTR
   CHARACTER(LEN=7) :: STATSTR
   LOGICAL CHECK,FEXIST
   CHARACTER(LEN=2) :: cios
   integer :: ios

   IF(IOPT(1:1) == "c")THEN  
     STATSTR = "old"
     CHECK = .TRUE.
   ELSE IF(IOPT(1:1) == "o") THEN 
     STATSTR = "unknown"
     CHECK = .FALSE.
   ELSE
     CALL Fatal_Error("FIRST LETTER IN FOPEN OPTION STRING MUST BE 'c' OR 'o'")
   END IF

   IF(IOPT(2:2) == "f")THEN  
     FORMSTR = "formatted"
   ELSE IF(IOPT(2:2) == "u") THEN 
     FORMSTR = "unformatted"
   ELSE
     CALL FATAL_ERROR("ERROR PROCESSING FOPEN ON FILE",INSTR,"2ND LETTER IN FOPEN OPTION STRING MUST BE 'f' OR 'u'")
   END IF

   IF(CHECK)THEN
     INQUIRE(FILE=INSTR,EXIST=FEXIST)
     IF(.NOT. FEXIST)  CALL FATAL_ERROR("FILE "//INSTR//" NOT FOUND")
   END IF

   OPEN(IUNIT,FILE=INSTR,STATUS=TRIM(STATSTR),FORM=TRIM(FORMSTR),IOSTAT=ios) 


   write(cios,'(i2.2)') ios
   if (ios == 9) then
      CALL FATAL_ERROR("UNABLE TO OPEN THE FILE:",&
           & INSTR, "IOSTAT ERROR#"//cios//"; suggests bad permissions ?")
      
   elseif (ios ==29) then
      CALL FATAL_ERROR("UNABLE TO OPEN THE FILE:",&
           & INSTR, "IOSTAT ERROR#"//cios//"; suggests bad directory path ?")
      
      
   elseif (IOS /= 0)then
      Call FATAL_ERROR("UNABLE TO OPEN THE FILE:",INSTR,"IOSTAT ERROR# &
           &"//CIOS//"; UNNKOWN ERROR ?")
   END IF


   IF(IOPT(3:3) == "r")  REWIND(IUNIT)

   if(DBG_SET(dbg_io))  &
        & write(IPT,*) "Opend File: ",INSTR
   


END SUBROUTINE FOPEN


INTEGER FUNCTION OPEN_DAT(FNAME,UNIT,PATH)
  USE CONTROL, only : INPUT_DIR
  implicit none
  CHARACTER(LEN=*) :: FNAME
  INTEGER :: UNIT
  CHARACTER(LEN=*), OPTIONAL :: PATH
  
  CHARACTER(LEN=400) :: PATHNFILE
  
  OPEN_DAT = -1
  IF (LEN_TRIM(FNAME) ==0) return
  
  IF (PRESENT(PATH)) THEN
     IF (LEN_TRIM(PATH) ==0) return
     
     pathnfile = trim(PATH)//trim(FNAME)
  ELSE
     IF (LEN_TRIM(INPUT_DIR) ==0) return
     pathnfile = trim(INPUT_DIR)//trim(FNAME)
  END IF
  Call FOPEN(UNIT,trim(pathnfile),'cfr')
  OPEN_DAT = 0
  
  
END FUNCTION OPEN_DAT



!==============================================================================!
!  DECOMPOSE INPUT LINE INTO VARIABLE NAME AND VARIABLE VALUE(S)               !
!==============================================================================!

SUBROUTINE GET_VALUE(LNUM,NUMCHAR,TEXT_LINE,VARNAME,VARTYPE,LOGVAL,STRINGVAL,&
                    REALVAL,INTVAL,NVAL)

!==============================================================================!
  USE MOD_PREC
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: LNUM,NUMCHAR
  CHARACTER(LEN=NUMCHAR) :: TEXT_LINE
  CHARACTER(LEN=40), INTENT(OUT) :: VARNAME
  CHARACTER(LEN=7), INTENT(OUT) :: VARTYPE
  LOGICAL, INTENT(OUT) :: LOGVAL
  CHARACTER(LEN=80), INTENT(OUT) :: STRINGVAL(150)
  REAL(DP), INTENT(INOUT) :: REALVAL(150)
  INTEGER, INTENT(INOUT) :: INTVAL(150)
  INTEGER, INTENT(OUT) :: NVAL
!------------------------------------------------------------------------------!
  CHARACTER(LEN=NUMCHAR) :: VARVAL,TEMP,FRAG(200)
  CHARACTER(LEN=80) :: TSTRING
  CHARACTER(LEN=3) :: ERRSTRING
  CHARACTER(LEN=16) :: NUMCHARS 
  INTEGER LENGTH,EQLOC,LVARVAL,DOTLOC
  INTEGER I,J,LOCEX,NP
  LOGICAL ONFRAG

!==============================================================================!
  FRAG = " "
  NUMCHARS = "0123456789+-Ee. " 
  VARTYPE = "error"
  LOGVAL = .FALSE.
  LENGTH = LEN_TRIM(TEXT_LINE) 
  WRITE(ERRSTRING,"(I3)") LNUM
  LOCEX = INDEX(TEXT_LINE,"!")

!
!-----------------------CHECK FOR BLANK LINE OR COMMENT------------------------!
!
  IF(LENGTH == 0 .OR. LOCEX==1)THEN
    VARTYPE = "no data"
    VARNAME = "no data"
    RETURN
  END IF

!
!-----------------------CHANGE COMMAS TO BLANKS--------------------------------!
!
  DO I=1,LENGTH
    IF(TEXT_LINE(I:I) == ",") TEXT_LINE(I:I) = " "
  END DO
!
!-----------------------REMOVING TRAILING COMMENTS-----------------------------!
!
  IF(LOCEX /= 0)THEN
    TEMP = TEXT_LINE(1:LOCEX-1)
    TEXT_LINE = TEMP
   END IF
!
!--------------------ENSURE "=" EXISTS AND DETERMINE LOCATION------------------!
!
   EQLOC = INDEX(TEXT_LINE,"=")
   IF(EQLOC == 0) then
      CALL WARNING(& 
        &'Could not find correct variable name in the datafile header',&
        &'Header comment lines must start with "!", Data lines must contain "="',&
        &'DATA LINE '//ERRSTRING//': This often occurs if the header variable is missing!')
      VARTYPE = "no data"
      VARNAME = "no data"
      RETURN
   END IF

!
!--------------------SPLIT OFF VARNAME AND VARVAL STRINGS----------------------!
!
   VARNAME = TEXT_LINE(1:EQLOC-1)
   VARVAL  = ADJUSTL(TEXT_LINE(EQLOC+1:LENGTH))
   LVARVAL = LEN_TRIM(VARVAL)
   IF(LVARVAL == 0) then
      CALL WARNING('IN DATA PARAMETER FILE', &
           & 'VARIABLE'//VARNAME//'; LINE'//ERRSTRING//' HAS NO ASSOCIATED VALUE')
      VARTYPE = "no data"
      VARNAME = "no data"
      RETURN
   END IF

!-----------------DETERMINE TYPE OF VARVAL-------------------------------------!
!

!
!  CHECK FOR LOGICAL
!
   IF((VARVAL(1:1) == "T" .OR. VARVAL(1:1) == "F") .AND. LVARVAL == 1)THEN 
     VARTYPE = "logical"
     IF(VARVAL(1:1) == "T") LOGVAL = .TRUE.
     RETURN
   END IF

!
!  CHECK IF IT IS A STRING  (CONTAINS CHARACTERS OTHER THAN 0-9,+,-,e,E,.)
!
   DO I=1,LVARVAL
     IF(INDEX(NUMCHARS,VARVAL(I:I)) == 0) VARTYPE = "string" 
   END DO

!
!  PROCESS STRING (MAY BE MULTIPLE)
!
   IF(VARTYPE == "string") THEN
     TSTRING = VARVAL
     STRINGVAL(1) = TSTRING 
     NVAL = 1
     ONFRAG = .TRUE.
     DO I=1,LVARVAL
       IF(VARVAL(I:I) /= " ")THEN
         FRAG(NVAL) = TRIM(FRAG(NVAL))//VARVAL(I:I)
         ONFRAG = .TRUE.
       ELSE
         IF(ONFRAG) NVAL = NVAL + 1
         ONFRAG = .FALSE.
       END IF
     END DO
     DO I=1,NVAL
       STRINGVAL(I+1) = TRIM(FRAG(I))
     END DO
     RETURN
   END IF

!
!  CHECK IF IT IS A FLOAT
!

   DOTLOC = INDEX(VARVAL,".")
   IF(DOTLOC /= 0) THEN
     VARTYPE = "float"
   ELSE
     VARTYPE = "integer"
   END IF
!
!-----------------FRAGMENT INTO STRINGS FOR MULTIPLE VALUES---------------------!
!
   NP = 1
   ONFRAG = .TRUE.
   DO I=1,LVARVAL
     IF(VARVAL(I:I) /= " ")THEN 
       FRAG(NP) = TRIM(FRAG(NP))//VARVAL(I:I)
       ONFRAG = .TRUE.
     ELSE
       IF(ONFRAG) NP = NP + 1
       ONFRAG = .FALSE.
     END IF
   END DO
!
!-----------------EXTRACT NUMBER(S) FROM CHARACTER STRINGS----------------------!
!
   
   NVAL = NP
   DO I=1,NP
     TEMP = TRIM(FRAG(I))
     IF(VARTYPE == "float") THEN 
       READ(TEMP,*)REALVAL(I)
     ELSE
       READ(TEMP,*)INTVAL(I)
     END IF
   END DO

END SUBROUTINE GET_VALUE


!==============================================================================|

    FUNCTION SCAN_FILE(UNIT,VNAME,ISCAL,FSCAL,IVEC,FVEC,CVEC,NSZE,CVAL,LVAL)           

!==============================================================================|
!   Scan an Input File for a Variable                                          |
!   RETURN VALUE:                                                              |
!        0 = FILE FOUND, VARIABLE VALUE FOUND                                  |
!       -1 = FILE DOES NOT EXIST OR PERMISSIONS ARE INCORRECT                  |
!       -2 = VARIABLE NOT FOUND OR IMPROPERLY SET                              |
!       -3 = VARIABLE IS OF DIFFERENT TYPE, CHECK INPUT FILE                   |
!       -4 = VECTOR PROVIDED BUT DATA IS SCALAR TYPE                           |
!       -5 = NO DATATYPE DESIRED, EXITING                                      |
!							                       |
!   REQUIRED INPUT:		        				       |
!        UNIT = File UNIT					               |
!        FSIZE = Length of Filename					       |
!                                                                              | 
!   OPTIONAL (MUST PROVIDE ONE)        					       | 
!        ISCAL = INTEGER SCALAR					               |
!        FSCAL = FLOAT SCALAR  						       | 
!        CVAL = CHARACTER VARIABLE                                             |
!        LVAL = LOGICAL VARIABLE                                               |
!        IVEC = INTEGER VECTOR **                                              |
!        FVEC = FLOAT VECTOR **                                                |
!        CVEC = STRING VECTOR ** (STRINGS OF LENGTH 80)                        |
!      **NSZE = ARRAY SIZE (MUST BE PROVIDED WITH IVEC/FVEC)                   |
!                                                                              | 
!==============================================================================|

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: UNIT
   CHARACTER(LEN=*) :: VNAME
   INTEGER, INTENT(INOUT), OPTIONAL :: ISCAL,IVEC(*)
   REAL(SP),INTENT(INOUT), OPTIONAL :: FSCAL,FVEC(*)
   CHARACTER(LEN=80), OPTIONAL      :: CVAL,CVEC(*)
   LOGICAL, INTENT(INOUT), OPTIONAL :: LVAL
   INTEGER, INTENT(INOUT), OPTIONAL :: NSZE 
   
!------------------------------------------------------------------------------|
   INTEGER :: SCAN_FILE
   REAL(DP) REALVAL(150)
   INTEGER  INTVAL(150)
   CHARACTER(LEN=40 ) :: VARNAME
   CHARACTER(LEN=80 ) :: STRINGVAL(150),TITLE
   CHARACTER(LEN=400 ) :: INPLINE
   CHARACTER(LEN=800) :: TLINE
   CHARACTER(LEN=7  ) :: VARTYPE, ENDLINE
   CHARACTER(LEN=20 ), DIMENSION(200)  :: SET
   INTEGER I,NVAL,J,NSET,NLINE,NREP,bgn,nd, LEL
   LOGICAL SETYES,ALLSET,CHECK,LOGVAL


   SCAN_FILE = 0

!==============================================================================|
!            SCAN THE FILE FOR THE VARIABLE NAME                               |
!==============================================================================|

# if defined (COMPAQ) || defined (IFORT) || defined (IRIX)
   ENDLINE='\\'
# else
   ENDLINE= '\\\\'
# endif
   LEL = LEN_TRIM(ENDLINE)-1

   REWIND(UNIT)

   NSET = 0
   NLINE = 0
   DO WHILE(.TRUE.)

      if(NLINE >200) then
         CALL Warning("Read 200 lines of header with out finding parameters! ")
         SCAN_FILE=-2
         return
      end if

     TLINE(1:LEN(TLINE)) = ' ' 
     NREP  = 0
     NLINE = NLINE + 1
     READ(UNIT,'(a)',END=20) INPLINE
     TLINE = TRIM(INPLINE)

!----PROCESS LINE CONTINUATIONS------------------------------------------------!
     DO
!        write(ipt,*) '"'//TRIM(INPLINE)//'"'
        I = LEN_TRIM(INPLINE)
        IF(I > 1)THEN

           IF( INPLINE(I-LEL:I) == trim(endline))THEN
              
              NREP = NREP + 1
              READ(UNIT,'(a)',END=20) INPLINE
              NLINE = NLINE + 1
              bgn = LEN_TRIM(TLINE)+1
              nd = bgn +len_trim(INPLINE)
              
              TLINE(bgn:nd) = TRIM(INPLINE)
           ELSE
              EXIT
              
           END IF
        ELSE
           EXIT
        END IF
     END DO

!----REMOVE LINE CONTINUATION CHARACTER \\-------------------------------------!
     IF(NREP > 0)THEN
       DO I=2,LEN_TRIM(TLINE)
         IF( TLINE(I-LEL:I) == ENDLINE) TLINE(I-LEL:I) = '  '

       END DO
     END IF
       
!----PROCESS THE LINE----------------------------------------------------------!
     CALL GET_VALUE(NLINE,LEN_TRIM(TLINE),ADJUSTL(TLINE),VARNAME,VARTYPE,LOGVAL,&
                 STRINGVAL,REALVAL,INTVAL,NVAL)

!----IF VARNAME MATCHES, PROCESS VARIABLE AND ERROR-CHECK----------------------!

     IF(TRIM(VARNAME) == TRIM(VNAME))THEN

       IF(PRESENT(ISCAL))THEN
          IF(VARTYPE == 'integer')THEN
           ISCAL = INTVAL(1)
           RETURN
         ELSE
           SCAN_FILE = -3
           return
         END IF
       ELSE IF(PRESENT(FSCAL))THEN
         IF(VARTYPE == 'float')THEN
           FSCAL = REALVAL(1)
           RETURN
         ELSE
           SCAN_FILE = -3
           return
         END IF
       ELSE IF(PRESENT(CVAL))THEN
         IF(VARTYPE == 'string')THEN
           CVAL = STRINGVAL(1) 
           RETURN
         ELSE
           SCAN_FILE = -3
           return
         END IF
       ELSE IF(PRESENT(LVAL))THEN
         IF(VARTYPE == 'logical')THEN
           LVAL = LOGVAL 
           RETURN
         ELSE
           SCAN_FILE = -3
           return
         END IF
       ELSE IF(PRESENT(IVEC))THEN
         IF(NVAL > 1)THEN
           IF(VARTYPE == 'integer')THEN
             IVEC(1:NVAL) = INTVAL(1:NVAL) 
             NSZE = NVAL         
             RETURN
           ELSE
             SCAN_FILE = -3
             return
           END IF
           ELSE
           SCAN_FILE = -4
           return
         END IF
       ELSE IF(PRESENT(FVEC))THEN
         IF(NVAL > 1)THEN
           IF(VARTYPE == 'float')THEN
             FVEC(1:NVAL) = REALVAL(1:NVAL) 
             NSZE = NVAL           
             RETURN
           ELSE
             SCAN_FILE = -3
             return
           END IF
         ELSE
           SCAN_FILE = -4
           return
         END IF
       ELSE IF(PRESENT(CVEC))THEN
         IF(NVAL > 0)THEN
           IF(VARTYPE == 'string')THEN
             CVEC(1:NVAL) = STRINGVAL(2:NVAL+1)
             NSZE = NVAL 
             RETURN
           ELSE
             SCAN_FILE = -3
             return
           END IF
         ELSE
           SCAN_FILE = -4
           return
         END IF
       ELSE
         SCAN_FILE = -5
         return
       END IF
     END IF  !!VARIABLE IS CORRECT
            
   END DO !!LOOP OVER INPUT FILE

20   SCAN_FILE = -2

   RETURN 
   END FUNCTION SCAN_FILE


   SUBROUTINE SPLIT_STRING(instring,delim, outstrings)
     IMPLICIT NONE
     character(len=*), intent(in) :: instring
     character, intent(in) :: delim
     character(len=*), intent(OUT), ALLOCATABLE :: outstrings(:)
     integer :: nlen, i, cnt, prev,next, idx, outlen, lgn

!     character(len=len(outstrings)), ALLOCATABLE :: out_temp(:)
     character(len=len(instring)), ALLOCATABLE :: out_temp(:)
    

     ! Get the length of the string
     lgn = len_trim(instring)
     outlen = len(outstrings)
 
     ! CHECK FOR DEGERNERATE CASE (EMPTY STRING!)
     IF(LGN==0) THEN
        ALLOCATE(outSTRINGS(1))
        outStrings=""
        RETURN
     END IF


     ! Count the number of seperations
     cnt = 0
     do i = 1,lgn
        if(instring(I:I) == delim) cnt=cnt+1
     end do
     ! If the string is not terminated, count the last entry too...
     if(instring(lgn:lgn) /= delim) cnt=cnt+1
     
     ! Allocate space
     ALLOCATE(outSTRINGS(CNT))
     
     ! Split the string
     prev=1
     next=0
     DO I = 1,CNT
!        write(*,*) "*** '"//trim(instring(prev:lgn))//"'"
        
        ! Find the first seperation
        idx = index(instring(prev:lgn),delim) 
        if(idx==0) then
           ! IF none found, use end of string
           idx = lgn+1
        else
           ! Get the index into the real string length
           idx =idx + prev-1
        end if
        
        ! Set that last value to take
        next = idx-1
!        write(*,*) I, prev, idx, next

        if(outlen .le. next-prev) Call WARNING&
             ("Insufficent room to split string!")

        ! Copy it into the string array
        outstrings(I) = TRIM(adjustl(instring(prev:next)))

        if(outlen .le. next-prev) Call WARNING&
             ("Insufficent room to split string!","'"//trim(outstrings(I))//"'")

        
        ! Set the first character of the next string
        prev=idx+1
        
!        write(*,*) "! '"//trim(strings(I))//"'"
     END DO
     
     ! REMOVE DEGENERATE CASE FOR 'space' delimiter

     IF(delim == ' ') THEN
        
        idx = 0
        ALLOCATE(OUT_TEMP(CNT))
        DO I = 1,CNT
           
           IF(len_trim(outstrings(I)) > 0) THEN
              idx = idx + 1
              OUT_TEMP(idx) = outstrings(I)
           END IF
           
        END DO

        DEALLOCATE(outstrings)
        allocate(outstrings(idx))
        outstrings=out_temp


     END IF

   END SUBROUTINE SPLIT_STRING

   ! UTILITY TO SPLIT A PATH AND FILE NAME INTO THE DIRECTORY, THE
   ! FILE NAME AND THE FILE EXTENSION
   SUBROUTINE PATH_SPLIT(STRING,PATH,FILE,EXTENSION)
     IMPLICIT NONE
     
     CHARACTER(LEN=*), INTENT(IN) :: STRING
     CHARACTER(LEN=*), INTENT(OUT):: PATH,FILE,EXTENSION

     INTEGER :: IDX, LGN, I

     LGN = LEN_TRIM(STRING)

     PATH = ''
     FILE = ''
     EXTENSION = ''

     IDX = 1
     DO I=1,LGN
        if(string(I:I) == '/') idx = i
     END DO
     

     IF(IDX>1) THEN
        PATH = string(1:idx)
        idx = idx+1
     ELSE ! HANDLE THE CASE OF NO PATH
        PATH = './'
     END IF

     
     DO I = idx,LGN
        IF(string(I:I) == '.') THEN
           FILE = STRING(idx:(I-1))
           EXTENSION =STRING(I:LGN)
           EXIT
        END IF
     END DO

     ! HANDLE THE CASE OF NO EXTENSION
     IF (LEN_TRIM(FILE) == 0) FILE = STRING(idx:LGN)


   END SUBROUTINE PATH_SPLIT


   SUBROUTINE TEST_SPLIT_STRINGS
     implicit none
     CHARACTER(len=200) TESTIN
     CHARACTER(len=50),allocatable::testout(:)
     integer I
     
!     TESTIN = "Hello world"
!     TESTIN = "Hello world, Hello world2"
!     TESTIN = "Hello worldmore that fiftymore that fiftymore that fiftymore that fiftymore that fiftymore that fiftymore that fiftymore that fifty,, Hello world2,"
     TESTIN = "Hello world, Hello world2, Hello world2, Hello world2, Hello world2"
     TESTIN = "Hello world"//achar(10)//" Hello world2"//achar(10)//" Hello world2, Hello world2, Hello world2"
     
     call split_string(testin,achar(10),testout)
     
     write(ipt,*) "! "
     write(ipt,*) "! TESTING SPLIT STRINGS"
     write(ipt,*) "! "
     
     do i=1,size(testout)
        write(ipt,*) "! ",I,"'"//trim(testout(I))//"'"
     end do
     
   END SUBROUTINE TEST_SPLIT_STRINGS

!==========================================================================
! Calculate LED Limiter L(u,v)
!==========================================================================
  Real(sp) Function LimLED(a,b,q) Result(lim)

  IMPLICIT NONE

  real(sp) a,b
  real(sp) q,R
  real(sp) eps
  eps = epsilon(eps)

 ! exponent
 ! q = 0. !1st order
 ! q = 1. !minmod
 ! q = 2. !van leer

  R = abs(   (a-b)/(abs(a)+abs(b)+eps) )**q
  lim = .5*(1-R)*(a+b)

  End Function LimLED

  Real(sp) Function LimLED1(a,b,alpha) Result(lim)
  IMPLICIT NONE

  real(sp) a,b,alpha

  lim = 0.5_sp*(sign(1.,a)+sign(1.,b))*max(min(alpha*abs(a),abs(b)),min(abs(a),alpha*abs(b)))

  End Function LimLED1

  Real(sp) Function LimLED2(a,b,alpha) Result(lim)
  IMPLICIT NONE

  real(sp) a,b,alpha

  lim = 0.5_sp*(sign(1.,a)+sign(1.,b))*min(0.5_sp*abs(a+b),alpha*abs(a),alpha*abs(b))

  End Function LimLED2

  REAL(DP) FUNCTION READ_FLOAT(ITEM,IERR)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: ITEM
    INTEGER, INTENT(OUT) :: IERR
    
    LOGICAL :: ISFLOAT
    INTEGER :: I

    IERR = -1
    ISFLOAT = .FALSE.
    READ_FLOAT = -99999.9
    DO I = 1,Len_trim(ITEM)
       IF (ITEM(I:I) == ".") THEN
          ISFLOAT = .true.
          CYCLE
       END IF
       
       IF( LGT("0",ITEM(I:I)) .or. LLT("9",ITEM(I:I)) ) THEN
          ISFLOAT = .FALSE.
          EXIT
       END IF
       
    END DO

    IF(ISFLOAT) THEN
       READ(ITEM,*,IOSTAT=ierr) READ_FLOAT
    END IF

  END FUNCTION READ_FLOAT

  REAL(SP) FUNCTION READ_INT(ITEM,IERR)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: ITEM
    INTEGER, INTENT(OUT) :: IERR
    
    INTEGER :: I

    IERR = -1
    READ_INT = -99999.9
    DO I = 1,Len_trim(ITEM)
       
       IF( LGT("0",ITEM(I:I)) .or. LLT("9",ITEM(I:I)) ) THEN
          RETURN
       END IF
       
    END DO

    READ(ITEM,*,IOSTAT=ierr) READ_INT
    
  END FUNCTION READ_INT

  FUNCTION SCAN_FILE2(FNAME,VNAME,ISCAL,FSCAL,IVEC,FVEC,CVEC,NSZE,CVAL,LVAL)           

!==============================================================================|
!   Scan an Input File for a Variable                                          |
!   RETURN VALUE:                                                              |
!        0 = FILE FOUND, VARIABLE VALUE FOUND                                  |
!       -1 = FILE DOES NOT EXIST OR PERMISSIONS ARE INCORRECT                  |
!       -2 = VARIABLE NOT FOUND OR IMPROPERLY SET                              |
!       -3 = VARIABLE IS OF DIFFERENT TYPE, CHECK INPUT FILE                   |
!       -4 = VECTOR PROVIDED BUT DATA IS SCALAR TYPE                           |
!       -5 = NO DATATYPE DESIRED, EXITING                                      |
!							                       |
!   REQUIRED INPUT:		        				       |
!        FNAME = File Name					               |
!        FSIZE = Length of Filename					       |
!                                                                              | 
!   OPTIONAL (MUST PROVIDE ONE)        					       | 
!        ISCAL = INTEGER SCALAR					               |
!        FSCAL = FLOAT SCALAR  						       | 
!        CVAL = CHARACTER VARIABLE                                             |
!        LVAL = LOGICAL VARIABLE                                               |
!        IVEC = INTEGER VECTOR **                                              |
!        FVEC = FLOAT VECTOR **                                                |
!        CVEC = STRING VECTOR ** (STRINGS OF LENGTH 80)                        |
!      **NSZE = ARRAY SIZE (MUST BE PROVIDED WITH IVEC/FVEC)                   |
!                                                                              | 
!==============================================================================|

   USE MOD_PREC
   USE OCPCOMM4
   IMPLICIT NONE
   CHARACTER(LEN=*) :: FNAME,VNAME
   INTEGER, INTENT(INOUT), OPTIONAL :: ISCAL,IVEC(*)
   REAL(SP),INTENT(INOUT), OPTIONAL :: FSCAL,FVEC(*)
   CHARACTER(LEN=80), OPTIONAL      :: CVAL,CVEC(*)
   LOGICAL, INTENT(INOUT), OPTIONAL :: LVAL
   INTEGER, INTENT(INOUT), OPTIONAL :: NSZE 
   
!------------------------------------------------------------------------------|

   INTEGER :: SCAN_FILE2
   REAL(DP) REALVAL(150)
   INTEGER  INTVAL(150)
   CHARACTER(LEN=40 ) :: VARNAME
   CHARACTER(LEN=80 ) :: STRINGVAL(150),TITLE
   CHARACTER(LEN=80 ) :: INPLINE
   CHARACTER(LEN=400) :: TLINE
   CHARACTER(LEN=7  ) :: VARTYPE
   CHARACTER(LEN=20 ), DIMENSION(200)  :: SET
   INTEGER I,NVAL,J,NSET,NLINE,NREP
   LOGICAL SETYES,ALLSET,CHECK,LOGVAL


   SCAN_FILE2 = 0
!==============================================================================|
!            OPEN THE INPUT FILE                                               |
!==============================================================================|
   INQUIRE(FILE=TRIM(FNAME),EXIST=CHECK)
   IF(.NOT.CHECK)THEN
     SCAN_FILE2 = -1
     RETURN
   END IF

   OPEN(INPUTF,FILE=TRIM(FNAME)) ; REWIND(INPUTF) 

!==============================================================================|
!            SCAN THE FILE FOR THE VARIABLE NAME                               |
!==============================================================================|

   NSET = 0
   NLINE = 0
   DO WHILE(.TRUE.)
     TLINE(1:LEN(TLINE)) = ' ' 
     NREP  = 0
     NLINE = NLINE + 1
     READ(INPUTF,'(a)',END=20) INPLINE
     TLINE(1:80) = INPLINE(1:80)

!----PROCESS LINE CONTINUATIONS------------------------------------------------!
 110 CONTINUE
     I = LEN_TRIM(INPLINE)
     IF(I /= 0)THEN
#    if defined (COMPAQ) || defined (INTEL) || defined (IRIX)
     IF( INPLINE(I-1:I) == '\\')THEN
#    else
     IF( INPLINE(I-1:I) == '\\\\')THEN
#    endif
       NREP = NREP + 1
       READ(INPUTF,'(a)',END=20) INPLINE
       NLINE = NLINE + 1
       TLINE( NREP*80 + 1 : NREP*80 +80) = INPLINE(1:80)
       GOTO 110
     END IF
     END IF
     IF(NREP > 4)CALL Fatal_Error("CANNOT HAVE > 4 LINE CONTINUATIONS")

!----REMOVE LINE CONTINUATION CHARACTER \\-------------------------------------!
     IF(NREP > 0)THEN
       DO I=2,LEN_TRIM(TLINE)
#    if defined (COMPAQ) || defined (INTEL) || defined (IRIX)
         IF( TLINE(I-1:I) == '\\') TLINE(I-1:I) = '  '
#        else
         IF( TLINE(I-1:I) == '\\\\') TLINE(I-1:I) = '  '
#        endif
       END DO
     END IF
       
!----PROCESS THE LINE----------------------------------------------------------!
     CALL GET_VALUE(NLINE,LEN_TRIM(TLINE),ADJUSTL(TLINE),VARNAME,VARTYPE,LOGVAL,&
                 STRINGVAL,REALVAL,INTVAL,NVAL)

!----IF VARNAME MATCHES, PROCESS VARIABLE AND ERROR-CHECK----------------------!

     IF(TRIM(VARNAME) == TRIM(VNAME))THEN

       IF(PRESENT(ISCAL))THEN
         IF(VARTYPE == 'integer')THEN
           ISCAL = INTVAL(1)
           RETURN
         ELSE
           SCAN_FILE2 = -3
         END IF
       ELSE IF(PRESENT(FSCAL))THEN
         IF(VARTYPE == 'float')THEN
           FSCAL = REALVAL(1)
           RETURN
         ELSE
           SCAN_FILE2 = -3
         END IF
       ELSE IF(PRESENT(CVAL))THEN
         IF(VARTYPE == 'string')THEN
           CVAL = STRINGVAL(1) 
           RETURN
         ELSE
           SCAN_FILE2 = -3
         END IF
       ELSE IF(PRESENT(LVAL))THEN
         IF(VARTYPE == 'logical')THEN
           LVAL = LOGVAL 
           RETURN
         ELSE
           SCAN_FILE2 = -3
         END IF
       ELSE IF(PRESENT(IVEC))THEN
         IF(NVAL > 1)THEN
           IF(VARTYPE == 'integer')THEN
             IVEC(1:NVAL) = INTVAL(1:NVAL) 
             NSZE = NVAL         
             RETURN
           ELSE
             SCAN_FILE2 = -3
           END IF
           ELSE
           SCAN_FILE2 = -4 
         END IF
       ELSE IF(PRESENT(FVEC))THEN
         IF(NVAL > 1)THEN
           IF(VARTYPE == 'float')THEN
             FVEC(1:NVAL) = REALVAL(1:NVAL) 
             NSZE = NVAL           
             RETURN
           ELSE
             SCAN_FILE2 = -3
           END IF
         ELSE
           SCAN_FILE2 = -4 
         END IF
       ELSE IF(PRESENT(CVEC))THEN
         IF(NVAL > 0)THEN
           IF(VARTYPE == 'string')THEN
             CVEC(1:NVAL) = STRINGVAL(2:NVAL+1)
             NSZE = NVAL 
             RETURN
           ELSE
             SCAN_FILE2 = -3
           END IF
         ELSE
           SCAN_FILE2 = -4
         END IF
       ELSE
         SCAN_FILE2 = -5
       END IF
     END IF  !!VARIABLE IS CORRECT
            
   END DO !!LOOP OVER INPUT FILE
 20 CLOSE(INPUTF) 
   SCAN_FILE2 = -2

   RETURN 
   END FUNCTION SCAN_FILE2


    FUNCTION SCAN_FILE3(FNAME,VNAME,ISCAL,FSCAL,IVEC,FVEC,CVEC,NSZE,CVAL,LVAL)           

!==============================================================================|
!   Scan an Input File for a Variable                                          |
!   RETURN VALUE:                                                              |
!        0 = FILE FOUND, VARIABLE VALUE FOUND                                  |
!       -1 = FILE DOES NOT EXIST OR PERMISSIONS ARE INCORRECT                  |
!       -2 = VARIABLE NOT FOUND OR IMPROPERLY SET                              |
!       -3 = VARIABLE IS OF DIFFERENT TYPE, CHECK INPUT FILE                   |
!       -4 = VECTOR PROVIDED BUT DATA IS SCALAR TYPE                           |
!       -5 = NO DATATYPE DESIRED, EXITING                                      |
!							                       |
!   REQUIRED INPUT:		        				       |
!        FNAME = File Name					               |
!        FSIZE = Length of Filename					       |
!                                                                              | 
!   OPTIONAL (MUST PROVIDE ONE)        					       | 
!        ISCAL = INTEGER SCALAR					               |
!        FSCAL = FLOAT SCALAR  						       | 
!        CVAL = CHARACTER VARIABLE                                             |
!        LVAL = LOGICAL VARIABLE                                               |
!        IVEC = INTEGER VECTOR **                                              |
!        FVEC = FLOAT VECTOR **                                                |
!        CVEC = STRING VECTOR ** (STRINGS OF LENGTH 80)                        |
!      **NSZE = ARRAY SIZE (MUST BE PROVIDED WITH IVEC/FVEC)                   |
!                                                                              | 
!==============================================================================|

   USE MOD_PREC
   USE OCPCOMM4
   IMPLICIT NONE
   CHARACTER(LEN=*) :: FNAME,VNAME
   INTEGER, INTENT(INOUT), OPTIONAL :: ISCAL,IVEC(*)
   REAL(SP),INTENT(INOUT), OPTIONAL :: FSCAL,FVEC(*)
   CHARACTER(LEN=80), OPTIONAL      :: CVAL,CVEC(*)
   LOGICAL, INTENT(INOUT), OPTIONAL :: LVAL
   INTEGER, INTENT(INOUT), OPTIONAL :: NSZE 
   
!------------------------------------------------------------------------------|

   INTEGER :: SCAN_FILE3
   REAL(DP) REALVAL(150)
   INTEGER  INTVAL(150)
   CHARACTER(LEN=40 ) :: VARNAME
   CHARACTER(LEN=80 ) :: STRINGVAL(150),TITLE
   CHARACTER(LEN=80 ) :: INPLINE
   CHARACTER(LEN=400) :: TLINE
   CHARACTER(LEN=7  ) :: VARTYPE
   CHARACTER(LEN=20 ), DIMENSION(200)  :: SET
   INTEGER I,NVAL,J,NSET,NLINE,NREP
   LOGICAL SETYES,ALLSET,CHECK,LOGVAL


   SCAN_FILE3 = 0
!==============================================================================|
!            OPEN THE INPUT FILE                                               |
!==============================================================================|
   INQUIRE(FILE=TRIM(FNAME),EXIST=CHECK)
   IF(.NOT.CHECK)THEN
     SCAN_FILE3 = -1
     RETURN
   END IF

   OPEN(INPUTF,FILE=TRIM(FNAME)) ; REWIND(INPUTF) 

!==============================================================================|
!            SCAN THE FILE FOR THE VARIABLE NAME                               |
!==============================================================================|

   NSET = 0
   NLINE = 0
   DO WHILE(.TRUE.)
     TLINE(1:LEN(TLINE)) = ' ' 
     NREP  = 0
     NLINE = NLINE + 1
     READ(INPUTF,'(a)',END=20) INPLINE
     TLINE(1:80) = INPLINE(1:80)

!----PROCESS LINE CONTINUATIONS------------------------------------------------!
 110 CONTINUE
     I = LEN_TRIM(INPLINE)
     IF(I > 2)THEN
#    if defined (COMPAQ) || defined (INTEL) || defined (IRIX)
     IF( INPLINE(I-1:I) == '\\')THEN
#    else
     IF( INPLINE(I-1:I) == '\\')THEN
#    endif
       NREP = NREP + 1
       READ(INPUTF,'(a)',END=20) INPLINE
       NLINE = NLINE + 1
       TLINE( NREP*80 + 1 : NREP*80 +80) = INPLINE(1:80)
       GOTO 110
     END IF
     END IF
     IF(NREP > 4)CALL Fatal_Error("CANNOT HAVE > 4 LINE CONTINUATIONS")

!----REMOVE LINE CONTINUATION CHARACTER \\-------------------------------------!
     IF(NREP > 0)THEN
       DO I=2,LEN_TRIM(TLINE)
#    if defined (COMPAQ) || defined (INTEL) || defined (IRIX)
         IF( TLINE(I-1:I) == '\\') TLINE(I-1:I) = '  '
#        else
         IF( TLINE(I-1:I) == '\\') TLINE(I-1:I) = '  '
#        endif
       END DO
     END IF

      
!----PROCESS THE LINE----------------------------------------------------------!
     CALL GET_VALUE(NLINE,LEN_TRIM(TLINE),ADJUSTL(TLINE),VARNAME,VARTYPE,LOGVAL,&
                 STRINGVAL,REALVAL,INTVAL,NVAL)

!----IF VARNAME MATCHES, PROCESS VARIABLE AND ERROR-CHECK----------------------!

     IF(TRIM(VARNAME) == TRIM(VNAME))THEN

       IF(PRESENT(ISCAL))THEN
         IF(VARTYPE == 'integer')THEN
           ISCAL = INTVAL(1)
           CLOSE(INPUTF)
           RETURN
         ELSE
           SCAN_FILE3 = -3
         END IF
       ELSE IF(PRESENT(FSCAL))THEN
         IF(VARTYPE == 'float')THEN
           FSCAL = REALVAL(1)
           CLOSE(INPUTF)
           RETURN
         ELSE
           SCAN_FILE3 = -3
         END IF
       ELSE IF(PRESENT(CVAL))THEN
         IF(VARTYPE == 'string')THEN
           CVAL = STRINGVAL(1) 
           CLOSE(INPUTF)
           RETURN
         ELSE
           SCAN_FILE3 = -3
         END IF
       ELSE IF(PRESENT(LVAL))THEN
         IF(VARTYPE == 'logical')THEN
           LVAL = LOGVAL 
           CLOSE(INPUTF)
           RETURN
         ELSE
           SCAN_FILE3 = -3
         END IF
       ELSE IF(PRESENT(IVEC))THEN
         IF(NVAL > 1)THEN
           IF(VARTYPE == 'integer')THEN
             IVEC(1:NVAL) = INTVAL(1:NVAL) 
             NSZE = NVAL
             CLOSE(INPUTF)         
             RETURN
           ELSE
             SCAN_FILE3 = -3
           END IF
           ELSE
           SCAN_FILE3 = -4 
         END IF
       ELSE IF(PRESENT(FVEC))THEN
         IF(NVAL > 1)THEN
           IF(VARTYPE == 'float')THEN
             FVEC(1:NVAL) = REALVAL(1:NVAL) 
             NSZE = NVAL 
             CLOSE(INPUTF)          
             RETURN
           ELSE
             SCAN_FILE3 = -3
           END IF
         ELSE
           SCAN_FILE3 = -4 
         END IF
       ELSE IF(PRESENT(CVEC))THEN
         IF(NVAL > 0)THEN
           IF(VARTYPE == 'string')THEN
             CVEC(1:NVAL) = STRINGVAL(2:NVAL+1)
             NSZE = NVAL 
             CLOSE(INPUTF)
             RETURN
           ELSE
             SCAN_FILE3 = -3
           END IF
         ELSE
           SCAN_FILE3 = -4
         END IF
       ELSE
         SCAN_FILE3 = -5
       END IF
     END IF  !!VARIABLE IS CORRECT
            
   END DO !!LOOP OVER INPUT FILE
 20 CLOSE(INPUTF) 
   SCAN_FILE3 = -2

   RETURN 
   END FUNCTION SCAN_FILE3

!==============================================================================|


END module mod_utils

