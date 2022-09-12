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

Module Mod_Esmf_Nesting
# if defined (MULTIPROCESSOR)
  use all_vars
  use mod_utils
  use mod_par
  USE SWCOMM3, ONLY : MSC,MDC
  implicit none

  SAVE

  ! For data IO
  TYPE NEST_DATA
     REAL(SP), POINTER :: AC2_BLK(:,:,:)
  END TYPE NEST_DATA

  ! GRID AND DATA FOR RUNNING AS SUBDOMAIN
  TYPE(GRID), POINTER :: NESTING_GRID
  TYPE(NEST_DATA), POINTER :: NESTING_DATA

  !========================================================
  ! Used in NESTING CODE
  ! This memory holds the data for each variable at the current time.
  ! This is where the data that is updated is stored in arrays which
  ! are the size of the nesting boundary, and indexed to the FVCOM domain. 
  REAL(SP), ALLOCATABLE :: EL_NEST1(:)
  REAL(SP), ALLOCATABLE :: UA_NEST1(:)
  REAL(SP), ALLOCATABLE :: VA_NEST1(:)
  REAL(SP), ALLOCATABLE :: U_NEST1(:,:)
  REAL(SP), ALLOCATABLE :: V_NEST1(:,:)
  REAL(SP), ALLOCATABLE :: T_NEST1(:,:)
  REAL(SP), ALLOCATABLE :: S_NEST1(:,:)
!lwang
  REAL(SP), ALLOCATABLE :: ZZ_NEST1(:,:)
  REAL(SP), ALLOCATABLE :: ZZ1_NEST1(:,:)
!lwang
  REAL(SP), ALLOCATABLE :: ELT_NEST1(:)
  REAL(SP), ALLOCATABLE :: UAT_NEST1(:)
  REAL(SP), ALLOCATABLE :: VAT_NEST1(:)
  REAL(SP), ALLOCATABLE :: UT_NEST1(:,:)
  REAL(SP), ALLOCATABLE :: VT_NEST1(:,:)
  REAL(SP), ALLOCATABLE :: TT_NEST1(:,:)
  REAL(SP), ALLOCATABLE :: ST_NEST1(:,:)

  REAL(SP), ALLOCATABLE :: ETT(:)
  REAL(SP), ALLOCATABLE :: UAT(:)
  REAL(SP), ALLOCATABLE :: VAT(:)
  REAL(SP), ALLOCATABLE :: UT(:,:)
  REAL(SP), ALLOCATABLE :: VT(:,:)
  REAL(SP), ALLOCATABLE :: TT1(:,:)
  REAL(SP), ALLOCATABLE :: ST1(:,:)

  REAL(SP), ALLOCATABLE :: AC2_NEST1(:,:,:)
  REAL(SP), ALLOCATABLE :: AC2T_NEST1(:,:,:)
  !========================================================

  ! GRID AND DATA FOR OUTPUT
!JQI  INTEGER :: NCNEST_NUM
  INTEGER :: NEST_NUM, NEST_NUM_CELL
  CHARACTER(LEN=80), ALLOCATABLE :: NCNEST_FNAMES(:)
  TYPE(GRID), POINTER :: NEST_GRIDS(:)
  TYPE(NEST_DATA), POINTER :: NCNEST_DATA(:)
  
  LOGICAL, PRIVATE :: NEED_INIT_NEST = .TRUE.

  !--Parameters in NameList NML_NESTING
  LOGICAL ESMF_NESTING_ON
  LOGICAL TWO_WAY_ON
  INTEGER IN_FEEDBACK_DOMAIN
  CHARACTER(LEN=160) NESTING_NODE_FILES
  CHARACTER(LEN=160) NESTING_CELL_FILES

  NAMELIST /NML_ESMF_NESTING/        &
       & ESMF_NESTING_ON,            &
       & TWO_WAY_ON,                 &
       & IN_FEEDBACK_DOMAIN,         &
       & NESTING_NODE_FILES,         &
       & NESTING_CELL_FILES

   TYPE(MAP), POINTER, DIMENSION(:) :: E_NEST_MAP,N_NEST_MAP

   INTEGER, POINTER :: NESTID(:),ELMS_GL(:),NESTCID(:,:)

   TYPE(MAP), POINTER, DIMENSION(:) :: N_NEST_MAP1, N_NEST_MAP2 
   TYPE(MAP), POINTER, DIMENSION(:) :: N_NEST_MAP3, N_NEST_MAP4 

   TYPE(MAP), POINTER, DIMENSION(:) :: E_NEST_MAP1, E_NEST_MAP2
   TYPE(MAP), POINTER, DIMENSION(:) :: E_NEST_MAP3, E_NEST_MAP4

   INTEGER, POINTER :: NESTID1(:), NESTID2(:), NESTID3(:),NESTID4(:)
   INTEGER, POINTER :: NESTCID1(:,:), NESTCID2(:,:), NESTCID3(:,:),NESTCID4(:,:)

   INTEGER, POINTER :: NESTCELL1(:), NESTCELL2(:), NESTCELL3(:), NESTCELL4(:)

   INTEGER :: nnode, ncell

CONTAINS
  !==============================================================================!
  !
  !==============================================================================!
  SUBROUTINE NAME_LIST_INITIALIZE_NEST
    USE CONTROL

    IMPLICIT NONE

    !--Parameters in NameList NML_NESTING
    ESMF_NESTING_ON    = .FALSE.
    TWO_WAY_ON         = .FALSE.
    IN_FEEDBACK_DOMAIN = 0
    NESTING_NODE_FILES = "none"
    NESTING_CELL_FILES = "none"

    RETURN
  END SUBROUTINE NAME_LIST_INITIALIZE_NEST
  !==============================================================================!
  !  
  !==============================================================================!  
  SUBROUTINE NAME_LIST_PRINT_NEST
    USE CONTROL

    IMPLICIT NONE

    WRITE(UNIT=IPT,NML=NML_ESMF_NESTING)

    RETURN
  END SUBROUTINE NAME_LIST_PRINT_NEST
  !==============================================================================!
  !
  !==============================================================================!
  SUBROUTINE NAME_LIST_READ_NEST
    USE MOD_UTILS
    USE CONTROL

    IMPLICIT NONE

    INTEGER :: ios,I
    CHARACTER(LEN=120) :: FNAME

    IF(DBG_SET(dbg_sbr)) write(IPT,*) "Subroutine Begins: name_list_read_nest;"

    ios = 0
    FNAME = "./"//trim(casename)//"_run.nml"
    IF(DBG_SET(dbg_io)) write(IPT,*) "Get_nestpar: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NESTING FLAG
    REWIND(NMLUNIT)

    !READ NESTING FLAG
    READ(UNIT=NMLUNIT, NML=NML_ESMF_NESTING,IOSTAT=ios)  
    IF(ios /= 0)THEN
       IF(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_ESMF_NESTING)
       CALL FATAL_ERROR("Can Not Read NameList NML_ESMF_NESTING from file: "//trim(FNAME))
    END IF

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:NML_ESMF_NESTING"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_ESMF_NESTING)

    CLOSE(NMLUNIT)

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Ends: name_list_read_nest;"    

  END SUBROUTINE NAME_LIST_READ_NEST
  !==============================================================================!
  SUBROUTINE ALLOCATE_ESMF
  
  IMPLICIT NONE
  
  TYPE(GRID), POINTER :: G1
  INTEGER :: IOB, IOB_NODE
  
  IF(ESMF_NESTING_ON)THEN
 
    G1=>NEST_GRIDS(1)
  
    ALLOCATE(EL_NEST1(G1%MGL))
    ALLOCATE(UA_NEST1(G1%NGL))
    ALLOCATE(VA_NEST1(G1%NGL))
    ALLOCATE(U_NEST1(G1%NGL,KB))
    ALLOCATE(V_NEST1(G1%NGL,KB))
    ALLOCATE(T_NEST1(G1%MGL,KB))
    ALLOCATE(S_NEST1(G1%MGL,KB))

    ALLOCATE(ELT_NEST1(G1%MGL))
    ALLOCATE(UAT_NEST1(G1%NGL))
    ALLOCATE(VAT_NEST1(G1%NGL))
    ALLOCATE(UT_NEST1(G1%NGL,KB))
    ALLOCATE(VT_NEST1(G1%NGL,KB))
    ALLOCATE(TT_NEST1(G1%MGL,KB))
    ALLOCATE(ST_NEST1(G1%MGL,KB))

    ALLOCATE(AC2_NEST1(G1%MGL,MSC,MDC))
    ALLOCATE(AC2T_NEST1(G1%MGL,MSC,MDC))

!    ALLOCATE(ETT(0:MT))
!    ALLOCATE(UT(0:NT,KB))
!    ALLOCATE(VT(0:NT,KB))
!    ALLOCATE(TT1(0:MT,KB))
!    ALLOCATE(ST1(0:MT,KB))

!    DO IOB = 1,G1%MGL/2
!      if(G1%M /= 0)then
!        iob_node = NLID(NESTID1(iob))
!        if(iob_node /= 0)then
!          ISONB(iob_node) = 3
!        end if
!      end if
!    END DO

  ELSE

    ALLOCATE(EL_NEST1(0))
    ALLOCATE(UA_NEST1(0))
    ALLOCATE(VA_NEST1(0))
    ALLOCATE(U_NEST1(0,0))
    ALLOCATE(V_NEST1(0,0))
    ALLOCATE(T_NEST1(0,0))
    ALLOCATE(S_NEST1(0,0))

    ALLOCATE(ELT_NEST1(0))
    ALLOCATE(UAT_NEST1(0))
    ALLOCATE(VAT_NEST1(0))
    ALLOCATE(UT_NEST1(0,0))
    ALLOCATE(VT_NEST1(0,0))
    ALLOCATE(TT_NEST1(0,0))
    ALLOCATE(ST_NEST1(0,0))

    ALLOCATE(AC2_NEST1(0,0,0))
    ALLOCATE(AC2T_NEST1(0,0,0))
  END IF
  
  RETURN
  END SUBROUTINE ALLOCATE_ESMF
  !==============================================================================!
  SUBROUTINE SET_NEST_DOMAIN_ESMF
    
    IMPLICIT NONE
!-----JQI

CHARACTER(LEN=80), ALLOCATABLE :: NEST_FNAMES(:),NEST_FNAMES_CELL(:)
!JQIINTEGER :: NEST_NUM, NEST_NUM_CELL
!JQITYPE(GRID), POINTER :: NEST_GRIDS(:)

!-----JQI
    
    LOGICAL :: FOUND
    ! NESTING VARIABLES
    TYPE(GRID), POINTER :: GLOBAL_NG
    character(len=160) :: pathnfile
    
    ! NCNEST VARIABLES
    INTEGER :: STATUS,I
!JQI    INTEGER, POINTER :: NID(:)
    
    if(DBG_SET(dbg_sbr))  write(IPT,*) "START Setup_nest_domain"
    
    if(DBG_SET(dbg_log))  write(IPT,*) "! "
    
    if (dbg_set(dbg_LOG)) write(IPT,*) "!           SETTING UP NESTING IO"

    CALL SPLIT_STRING(NESTING_NODE_FILES,",",NEST_FNAMES)
    CALL SPLIT_STRING(NESTING_CELL_FILES,",",NEST_FNAMES_CELL)
    NEST_NUM = SIZE(NEST_FNAMES)
    NEST_NUM_CELL = SIZE(NEST_FNAMES_CELL)
    
    IF(NEST_NUM /= NEST_NUM_CELL) &
      & CALL FATAL_ERROR("THE NUMBER OF NESTING NODE FILES MUST EQUAL TO THE NUMBER OF NESTING CELL FILES!")
       
    ALLOCATE(NEST_GRIDS(NEST_NUM),STAT=status)
    IF(STATUS /=0) CALL FATAL_ERROR("COULD NOT ALLOCATE NEST_GRIDS!")
       
    if(DBG_SET(dbg_log))  write(IPT,*) "! READING NEST FILES:",NEST_NUM
    
    DO I = 1, NEST_NUM
          
      if(DBG_SET(dbg_log))  &
        & write(IPT,*) "! READING NEST FILE:"//TRIM(NEST_FNAMES(I))//" AND "//TRIM(NEST_FNAMES_CELL(I))
          
      ! NESTID IS ALLOCATED AND DEALLOCATED INTERNALLY!
!JQI      CALL LOAD_NESTING_NODE_FILE(NEST_FNAMES(I),NID)
      CALL LOAD_NESTING_NODE_FILE(NEST_FNAMES(I))
      CALL LOAD_NESTING_CELL_FILE(NEST_FNAMES_CELL(I))
          
!JQI      CALL GENMAP_NCNEST(NID,NCNEST_GRIDS(I))
      CALL GENMAP_NEST(NEST_GRIDS(I))
      
      if(i == 1)then
        N_NEST_MAP1 => N_NEST_MAP
        E_NEST_MAP1 => E_NEST_MAP
	NESTID1     => NESTID
	NESTCID1    => NESTCID
	NESTCELL1   => ELMS_GL
      end if	

      if(i == 2)then
        N_NEST_MAP2 => N_NEST_MAP
        E_NEST_MAP2 => E_NEST_MAP
	NESTID2     => NESTID
	NESTCID2    => NESTCID
	NESTCELL2   => ELMS_GL
      end if	

      if(i == 3)then
        N_NEST_MAP3 => N_NEST_MAP
        E_NEST_MAP3 => E_NEST_MAP
	NESTID3     => NESTID
	NESTCID3    => NESTCID
	NESTCELL3   => ELMS_GL
      end if	

      if(i == 4)then
        N_NEST_MAP4 => N_NEST_MAP
        E_NEST_MAP4 => E_NEST_MAP
	NESTID4     => NESTID
	NESTCID4    => NESTCID
	NESTCELL4   => ELMS_GL
      end if	

      write(*,*) "! SET NEST DOMAIN:"
      write(*,*) "! DIMENSIONS: MGL =",NEST_GRIDS(I)%MGL
      write(*,*) "! DIMENSIONS: NGL =",NEST_GRIDS(I)%NGL
      write(*,*) "! DIMENSIONS: M   =",NEST_GRIDS(I)%M
      write(*,*) "! DIMENSIONS: N   =",NEST_GRIDS(I)%N
      
!JQI      ALLOCATE(AC2_NEST(MDC,MSC,NEST_GRIDS(I)%MGL)); AC2_NEST = 0.0_SP
      

      if(DBG_SET(dbg_log))  write(IPT,*) "! SET NEST DOMAIN:"
      if(DBG_SET(dbg_log))  write(IPT,*) "! DIMENSIONS: MGL =",NEST_GRIDS(I)%MGL
      if(DBG_SET(dbg_log))  write(IPT,*) "! DIMENSIONS: NGL =",NEST_GRIDS(I)%NGL

    END DO
       
    if(DBG_SET(dbg_sbr))  write(IPT,*) "END Setup_nest_domain"
  END SUBROUTINE SET_NEST_DOMAIN_ESMF
  !==============================================================================!
  !
  !==============================================================================!
  SUBROUTINE LOAD_NESTING_NODE_FILE(FNAME)
    USE CONTROL
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: FNAME
!JQI    INTEGER, INTENT(OUT), POINTER :: nodes(:)

!JQI    INTEGER :: nnode
    INTEGER, POINTER :: ns(:)

    INTEGER :: I,J,n1,n2,n3
    CHARACTER(LEN=80) :: PATHNFILE, TEMP,TEMP2
    INTEGER :: ISCAN,IOS,sender,IERR

    IF(MSR)THEN

       !============================================
       !
       !  READ THE NESTING NODES
       !
       !============================================

       PATHNFILE = TRIM(INPUT_DIR)//TRIM(FNAME)
       CALL FOPEN(NESTUNIT,trim(pathnfile),'cfr')

       ISCAN = SCAN_FILE(NESTUNIT,"Node_Nest Number", ISCAL = nnode)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of NESTING NODE FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'The header must contain: "Node_Nest Number = "', &
               & 'Followed by an integer number of Nesting Nodes')
       END IF


       REWIND NESTUNIT

       DO WHILE(.TRUE.)
          READ(NESTUNIT,*,IOSTAT=IOS,END=199)N1,N2,N3
          if (IOS == 0) then
             BackSpace NESTUNIT
             exit
          end if

          CYCLE

199       Call FATAL_ERROR('Improper formatting of NESTING NODE FILE:',&
               &'Reached end of file with out finding NESTING data?',&
               &'FORMAT: NEST# SUBDOMAIN# LARGEDOMAIN# (ALL INTEGERS)')

       END DO
       
       ALLOCATE(NESTID(nnode))

       I = 0
       DO WHILE(.TRUE.)

          READ(NESTUNIT,*,IOSTAT=IOS) N1,N2,N3
          IF(IOS < 0) exit

          I = I + 1
          IF(I > NNODE) CALL FATAL_ERROR('Number of rows of data in the NESTING NODE file &
               &exceeds the stated number of nodes in the header ?')

          NESTID(I) = N2

       END DO

       CLOSE(NESTUNIT)

       ! CHECK TO MAKE SURE VALUES ARE REASONABLE!
       IF( 1 > MINVAL(NESTID)) THEN
          write(temp,'(I8)') MINLOC(NESTID)
          write(temp2,'(I8)') MGL
          CALL FATAL_ERROR('NESTING NODE NUMBER'//trim(temp)//&
               & ' IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE NESTING NODES ARE 1 <= '//trim(temp2))
       END IF

       IF(  MAXVAL(NESTID) > MGL) THEN
          write(temp,'(I8)') MAXLOC(NESTID)
          write(temp2,'(I8)') MGL
          CALL FATAL_ERROR('NESTING NODE NUMBER'//trim(temp)//&
               & ' IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE NESTING NODES ARE 1 <= '//trim(temp2))
       END IF

    END IF

    ! BROADCAST TO OTHER PROCS
    
    SENDER = MSRID -1 ! SEND FROM MASTER
    CALL MPI_BCAST(NNODE,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
    
    IF( .NOT. MSR) THEN
       ALLOCATE(NESTID(NNODE))
    END IF
    
    CALL MPI_BCAST(NESTID,NNODE,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)

     RETURN
  END SUBROUTINE LOAD_NESTING_NODE_FILE
  !==============================================================================!
  !==============================================================================!
  SUBROUTINE LOAD_NESTING_CELL_FILE(FNAME)
    USE CONTROL
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: FNAME

    INTEGER, POINTER :: ns(:)

    INTEGER :: I,J,c1,c2,c3
    CHARACTER(LEN=80) :: PATHNFILE, TEMP,TEMP2
    INTEGER :: ISCAN,IOS,sender,IERR
    INTEGER, PARAMETER :: NESTUNIT1=32

    IF(MSR)THEN

       !============================================
       !
       !  READ THE NESTING CELLS
       !
       !============================================

       PATHNFILE = TRIM(INPUT_DIR)//TRIM(FNAME)
       CALL FOPEN(NESTUNIT1,trim(pathnfile),'cfr')

       ISCAN = SCAN_FILE(NESTUNIT1,"Cell_Nest Number", ISCAL = ncell)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of NESTING CELL FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'The header must contain: "Cell_Nest Number = "', &
               & 'Followed by an integer number of Nesting Cells')
       END IF


       REWIND NESTUNIT1

       DO WHILE(.TRUE.)
          READ(NESTUNIT1,*,IOSTAT=IOS,END=199)C1,C2,C3
          if (IOS == 0) then
             BackSpace NESTUNIT1
             exit
          end if

          CYCLE

199       Call FATAL_ERROR('Improper formatting of NESTING CELL FILE:',&
               &'Reached end of file with out finding NESTING data?',&
               &'FORMAT: NEST# SUBDOMAIN# LARGEDOMAIN# (ALL INTEGERS)')

       END DO
       
       ALLOCATE(NESTCID(ncell,2))

       I = 0
       DO WHILE(.TRUE.)

          READ(NESTUNIT1,*,IOSTAT=IOS) C1,C2,C3
          IF(IOS < 0) exit

          I = I + 1
          IF(I > NCELL) CALL FATAL_ERROR('Number of rows of data in the NESTING CELL file &
               &exceeds the stated number of cells in the header ?')

          NESTCID(I,1) = C2
	  NESTCID(I,2) = C3

       END DO

       CLOSE(NESTUNIT1)

       ! CHECK TO MAKE SURE VALUES ARE REASONABLE!
       IF( 1 > MINVAL(NESTCID(:,1))) THEN
          write(temp,'(I8)') MINLOC(NESTCID(:,1))
          write(temp2,'(I8)') NGL
          CALL FATAL_ERROR('NESTING CELL NUMBER'//trim(temp)//&
               & ' IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE NESTING CELLS ARE 1 <= '//trim(temp2))
       END IF

       IF(  MAXVAL(NESTCID(:,1)) > NGL) THEN
          write(temp,'(I8)') MAXLOC(NESTCID(:,1))
          write(temp2,'(I8)') NGL
          CALL FATAL_ERROR('NESTING CELL NUMBER'//trim(temp)//&
               & ' IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE NESTING CELLS ARE 1 <= '//trim(temp2))
       END IF

    END IF

    ! BROADCAST TO OTHER PROCS
    
    SENDER = MSRID -1 ! SEND FROM MASTER
    CALL MPI_BCAST(NCELL,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
    
    IF( .NOT. MSR) THEN
       ALLOCATE(NESTCID(NCELL,2))
    END IF
    
    CALL MPI_BCAST(NESTCID,NCELL*2,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)

     RETURN
  END SUBROUTINE LOAD_NESTING_CELL_FILE
  !==============================================================================!
  !
  !==============================================================================!
SUBROUTINE GENMAP_NEST(NG)
!==============================================================================|
!
! CREATE A GLOBAL TO LOCAL MAP FOR DOMAIN NESTING OUTPUT
! USES DATA READ INTO: NID - THE NESTING BOUNDARY NODE STRING
!
!
!                     
! Creates:             MAP LINK LIST ENTRY FOR IO
!                      THE GRID FOR THE NESTING OUTPUT
!                          The transform from the local nesting list to the local domain
!
!==============================================================================|
  USE MOD_PAR
  USE MOD_TYPES
  USE LIMS
  USE CONTROL
  USE ALL_VARS
  IMPLICIT NONE

  TYPE(GRID), INTENT(OUT) :: NG !The Local Grid, The Local Nesting Grid
!JQI  INTEGER,INTENT(INOUT), POINTER :: NID(:)

  integer :: SENDER,RECVER, ierr, I, J,NCNT, NSZE, I1, status,CNT,CNT_L,K,lb,ub

  INTEGER, POINTER :: TEMP1(:),TEMP2(:)

  INTEGER, POINTER :: ELMS(:)     !, ELMS_GL(:)

!JQI  TYPE(MAP), POINTER, DIMENSION(:) :: E_MAP,N_MAP

  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "START: GENMAP_NCNEST"

  IF(.not.ASSOCIATED(NESTID)) CALL FATAL_ERROR&
       &('Called GENMAP_NCNEST, but NESTID is not associated!')

  ! FIND THE NESTING ELEMENTS
  IF(.not. IOPROC) THEN
     ALLOCATE(ELMS(0:NT)); ELMS=0 ! TEMPORARY STORAGE
     CNT=0
     DO I = 1,NT
        IF(ANY(NESTID==NGID_X(NV(I,1))))THEN ! FIND ANY CELL THAT HAS
          ! A LISTED NODE
           IF(ANY(NESTID==NGID_X(NV(I,2))) .and. ANY(NESTID==NGID_X(NV(I,3))))THEN
              ! IF ALL THE NODES IN THAT CELL ARE LISTED...
              CNT = CNT+1
              ELMS(CNT) = I
           END IF
        END IF
     END DO
     
     
     ! MUST COUNT ELEMENTS
     CNT_L=0
     DO I=1,CNT
        IF(ELMS(I)<=N) CNT_L = CNT_L+1
     END DO
  ELSE
     CNT_L = 0
     CNT = 0
  END IF
  
  ! SET DIMENSIONS...
  NG%MGL = ubound(NESTID,1) ! nodes are easy
   
  CALL MPI_ALLREDUCE(CNT_L,NG%NGL,1,MPI_INTEGER,MPI_SUM,MPI_FVCOM_GROUP,IERR)
  
  IF(NG%NGL < 1) CALL FATAL_ERROR&
       &("GENMAP_NCNEST: THERE IS A PROBLEM WITH YOUR NESTING FILE",&
       & "THERE ARE NO ELEMENTS SELCTED BASED ON THE NODES YOU SPECIFIED!")

!  WRITE(IPT,*) "MYID,LCL,LXL,GL",MYID,CNT_L,CNT,NG%NGL,NT

  IF(.not. IOPROC) THEN
     
     ! BUILD ONE GLOBAL LIST OF NESTING ELEMENTS
     ALLOCATE(ELMS_GL(NG%NGL))
     K = 1
     DO I = 1, NPROCS
        
        ! SEND THE NUMBER OF ELEMENTS YOU HAVE
        SENDER = I -1
        IF(MYID == I)THEN
           I1 = CNT_L
        END IF
        CALL MPI_BCAST(I1,1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
        
        ! COLLECT THE GLOBAL IDS
        IF(I1 > 0)THEN
           ALLOCATE(TEMP1(I1))
           IF(MYID == I)THEN
              CNT_L =0
              DO J =1,CNT
                 IF(ELMS(J)<=N) THEN
                    CNT_L = CNT_L+1
                    TEMP1(CNT_L)=EGID_X(ELMS(J))
                 END IF
              END DO
!              WRITE(IPT,*) "ELMS_LCL",TEMP1
              
           END IF
           
           
           CALL MPI_BCAST(TEMP1,I1,MPI_INTEGER,SENDER,MPI_FVCOM_GROUP,IERR)
           
           ELMS_GL(K:(K+I1-1))=TEMP1
           DEALLOCATE(TEMP1)
           
           K = K+I1
        END IF
        
     END DO
     
!     WRITE(IPT,*) "ELMS_GL",ELMS_GL
     
     DEALLOCATE(ELMS)
!============================================
! Make a list of the local Nesting nodes
!============================================
     
     !!SET UP LOCAL NESTING NODES
     ALLOCATE(TEMP1(0:NG%MGL));      TEMP1=0
     ALLOCATE(TEMP2(0:NG%MGL));      TEMP2=0
     
     NCNT = 0
     DO I=1,NG%MGL
        I1 = NLID( NESTID(I) )
        IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF BOUNDARY NODES
     NG%M = NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(NG%NGID(0:NG%M),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NCNEST: can not allocate:NG%NGID")

     ALLOCATE(NG%NLID(0:NG%M),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NCNEST: can not allocate:NG%NLID")

     NG%NLID = TEMP1(0:NCNT)
     NG%NGID = TEMP2(0:NCNT)

     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)

     !!SET UP LOCAL+HALO NESTING NODES
     ALLOCATE(TEMP1(0:NG%MGL));      TEMP1=0
     ALLOCATE(TEMP2(0:NG%MGL));      TEMP2=0
     
     NCNT = 0
     DO I=1,NG%MGL
        I1 = NLID_X( NESTID(I) )
        IF(I1  > M)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF BOUNDARY NODES
     NG%MT = NG%M + NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(NG%NGID_X(0:NG%MT),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NCNEST: can not allocate:NG%NGID")

     ALLOCATE(NG%NLID_X(0:NG%MT),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NCNEST: can not allocate:NG%NLID")

     NG%NLID_X(0:NG%M) = NG%NLID

     NG%NGID_X(0:NG%M) = NG%NGID

     lb = NG%M+1
     ub = NG%MT
     NG%NLID_X(lb:ub) = TEMP1(1:NCNT)

     NG%NGID_X(lb:ub) = TEMP2(1:NCNT)

     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)
!     DEALLOCATE(NID)

!============================================
! Make a list of the local Nesting elements
!============================================

     !!SET UP LOCAL NESTING ELEMENTS
     ALLOCATE(TEMP1(0:NG%NGL));      TEMP1=0
     ALLOCATE(TEMP2(0:NG%NGL));      TEMP2=0

     NCNT = 0
     DO I=1,NG%NGL
        I1 =  ELID(ELMS_GL(I)) 
        IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF BOUNDARY NODES
     NG%N = NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN

     ALLOCATE(NG%EGID(0:NG%N),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NCNEST: can not allocate:NG%EGID")

     ALLOCATE(NG%ELID(0:NG%N),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NCNEST: can not allocate:NG%ELID")

     NG%ELID = TEMP1(0:NCNT)
     NG%EGID = TEMP2(0:NCNT)

     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)



     !!SET UP LOCAL+HALO NESTING ELEMENTS
     ALLOCATE(TEMP1(0:NG%NGL));      TEMP1=0
     ALLOCATE(TEMP2(0:NG%NGL));      TEMP2=0

     NCNT = 0
     DO I=1,NG%NGL
        I1 =  ELID_X(ELMS_GL(I)) 
        IF(I1 > N)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF BOUNDARY NODES
     NG%NT = NG%N + NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN

     ALLOCATE(NG%EGID_X(0:NG%NT),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NCNEST: can not allocate:NG%EGID")

     ALLOCATE(NG%ELID_X(0:NG%NT),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NCNEST: can not allocate:NG%ELID")

     NG%ELID_X(0:NG%N) = NG%ELID
     NG%EGID_X(0:NG%N) = NG%EGID

     lb = NG%N+1
     ub = NG%NT
     NG%ELID_X(lb:ub) = TEMP1(1:NCNT)

     NG%EGID_X(lb:ub) = TEMP2(1:NCNT)
     
     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)
!JQI     DEALLOCATE(ELMS_GL)

  END IF

  !==============================================================================|
  !   SET UP ELEMENT MAPPING FOR GLOBAL 2 LOCAL TRANSFER OF BC'S                 | 
  !   BOUNDARY MAP :: BCMAP(NPROCS)                                              |
  !     BCMAP(1-->NPROCS)%NSIZE  :: NUMBER OF BOUNDARY NODES IN EACH DOM         |
  !     BCMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM |
  !==============================================================================|

  ! SET UP TRANSFER FROM GLOBAL NESTING BOUNDARY TO LOCAL NESTING
  ! BOUNDARY - NOT DIRECTLY TO THE LOCAL MESH

  ! ELEMENTS:

  ! ADD THE MAP TO THE LOCAL DOMAIN AND THE ONE2ONE MAP
  E_NEST_MAP => MAKE_MAP(MYID,NPROCS,NG%NGL,NT,NG%EGID_X,NG%ELID_X)
!  CALL PRINT_MAP(E_MAP,"GLOBAL 2 GRID")
  CALL ADD_MAP2LIST(INTERNAL_MAPS,E_NEST_MAP)
  NULLIFY(E_NEST_MAP)

  E_NEST_MAP => MAKE_MAP(MYID,NPROCS,NG%NGL,NG%NT,NG%EGID_X)
!  CALL PRINT_MAP(E_MAP,"GLOBAL 2 DATA")
  CALL ADD_MAP2LIST(INTERNAL_MAPS,E_NEST_MAP)
  NULLIFY(E_NEST_MAP)

  E_NEST_MAP => MAKE_MAP(MYID,NPROCS,NG%NGL,N,NG%EGID,NG%ELID)
!  CALL PRINT_MAP(E_MAP,"GLOBAL 2 GRID")
  CALL ADD_MAP2LIST(INTERNAL_MAPS,E_NEST_MAP)
  NULLIFY(E_NEST_MAP)

  E_NEST_MAP => MAKE_MAP(MYID,NPROCS,NG%NGL,NG%N,NG%EGID)
!  CALL PRINT_MAP(E_MAP,"GLOBAL 2 DATA")
  CALL ADD_MAP2LIST(INTERNAL_MAPS,E_NEST_MAP)
!  NULLIFY(E_MAP)



  ! NODES
  ! ADD THE MAP TO THE LOCAL DOMAIN AND THE ONE2ONE MAP
  N_NEST_MAP => MAKE_MAP(MYID,NPROCS,NG%MGL,MT,NG%NGID_X,NG%NLID_X)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_NEST_MAP)
!  NULLIFY(N_NEST_MAP)

  N_NEST_MAP => MAKE_MAP(MYID,NPROCS,NG%MGL,NG%MT,NG%NGID_X)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_NEST_MAP)
!  NULLIFY(N_NEST_MAP)

  N_NEST_MAP => MAKE_MAP(MYID,NPROCS,NG%MGL,M,NG%NGID,NG%NLID)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_NEST_MAP)
!  NULLIFY(N_NEST_MAP)

  N_NEST_MAP => MAKE_MAP(MYID,NPROCS,NG%MGL,NG%M,NG%NGID)
  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_NEST_MAP)
!  NULLIFY(N_NEST_MAP)

!  CALL PRINT_MAP(N_MAP,"TIME MAP")

!Q  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_MAP)
!Q  NULLIFY(N_MAP)
!Q  DEALLOCATE(TEMP1)

  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "END: GENMAP_NCNEST"   

  RETURN
END SUBROUTINE GENMAP_NEST
# endif
END Module Mod_Esmf_Nesting
