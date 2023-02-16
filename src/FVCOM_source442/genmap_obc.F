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

!==============================================================================|
SUBROUTINE GENMAP_OBC
!==============================================================================|
!      OBC Node Number
!      TS OBC Type
!      OBC Node list and types
!==============================================================================|
  USE MOD_OBCS
  USE BCS
  USE MOD_PAR
  USE LIMS
  USE CONTROL
  IMPLICIT NONE
  integer :: SENDER,RECVER, ierr, I, NCNT, NSZE, I1

  INTEGER, POINTER :: TEMP1(:),TEMP2(:)
  INTEGER, POINTER :: TEMP3(:),TEMP4(:)

  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "START: SETUP_OBCMAP"

  IF(SERIAL) THEN
     IOBCN    = IOBCN_GL

     ALLOCATE(I_OBC_N(IOBCN))
     I_OBC_N = I_OBC_GL
     ALLOCATE(TYPE_OBC(IOBCN))
     TYPE_OBC = TYPE_OBC_GL

     ALLOCATE(I_OBC_N_OUTPUT(IOBCN))
     I_OBC_N_OUTPUT = I_OBC_GL

     if (dbg_set(dbg_sbr)) &
          & write(IPT,*) "END: GENMAP_OBC - serial"
     return
  END IF

# if defined (MULTIPROCESSOR)

  ALLOCATE(TEMP1(IOBCN_GL))
  ALLOCATE(TEMP2(IOBCN_GL))
  ALLOCATE(TEMP3(IOBCN_GL))
  
  IOBCN = 0
  NCNT = 0
  IF(.not. IOPROC) THEN
     !!SET UP LOCAL OPEN BOUNDARY NODES
     DO I=1,IOBCN_GL
        I1 = NLID( I_OBC_GL(I) )
        IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = TYPE_OBC_GL(I)
           TEMP3(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF BOUNDARY NODES
     IOBCN = NCNT
     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(I_OBC_N(NCNT),TYPE_OBC(NCNT),I_OBC_N_OUTPUT(NCNT))
     I_OBC_N  = TEMP1(1:NCNT)
     TYPE_OBC = TEMP2(1:NCNT)
     I_OBC_N_OUTPUT = NGID(I_OBC_N)
  END IF


  !==============================================================================|
  !   SET UP ELEMENT MAPPING FOR GLOBAL 2 LOCAL TRANSFER OF BC'S                 | 
  !   BOUNDARY MAP :: BCMAP(NPROCS)                                              |
  !     BCMAP(1-->NPROCS)%NSIZE  :: NUMBER OF BOUNDARY NODES IN EACH DOM         |
  !     BCMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM |
  !==============================================================================|

  TEMP4 => TEMP3(1:NCNT)
  BCMAP => MAKE_MAP(MYID,NPROCS_TOTAL,IOBCN_GL,NCNT,TEMP4)
  NULLIFY(TEMP4)

!!$  ALLOCATE(BCMAP(NPROCS)); BCMAP(:)%NSIZE=0
!!$
!!$  !--Determine Number of Elements for Each Processor
!!$  DO I=1,NPROCS    
!!$     IF(MYID == I) BCMAP(I)%NSIZE = NCNT
!!$     SENDER = I - 1
!!$     CALL MPI_BCAST(BCMAP(I)%NSIZE,1,MPI_INTEGER,SENDER,MPI_COMM_FVCOM,IERR)
!!$  END DO
!!$
!!$  !--Allocate Mapping Array for Each Processor
!!$  DO I=1,NPROCS
!!$     ALLOCATE(BCMAP(I)%LOC_2_GL(0:BCMAP(I)%NSIZE))
!!$     BCMAP(I)%LOC_2_GL(0) = 0
!!$  END DO
!!$
!!$  !--Construct Mapping Array for Each Processor 
!!$  DO I=1,NPROCS
!!$     NSZE = BCMAP(I)%NSIZE
!!$     if(myid == I) BCMAP(I)%LOC_2_GL(1:NSZE) =  TEMP3(1:NSZE)
!!$
!!$     CALL MPI_BCAST(BCMAP(I)%LOC_2_GL(1:NSZE),NSZE,MPI_INTEGER,I-1,MPI_COMM_FVCOM,IERR)
!!$
!!$  END DO

  DEALLOCATE(TEMP1,TEMP2,TEMP3)

  ALLOCATE(TEMP1(IOBCN_GL))
  TEMP1 = I_OBC_GL


  ALLOCATE(TEMP3(NCNT))
  TEMP3 = 0 ! THIS ASSIGNMENT MAY CAUSE A PROBLEM ON THE IOPROC WHERE NCNT == 0
  
  SENDER =MSRID
  if (USE_MPI_IO_MODE) SENDER = IOPROCID ! TEST DEAL FROM IOPROC
!!$  CALL PDEAL(MYID,SENDER,NPROCS,BCMAP,TEMP1,TEMP3)
  CALL PDEAL_IO(MYID,SENDER,NPROCS_TOTAL,BCMAP,TEMP1,TEMP3)

  if (.not. IOPROC) then
     DO I =1, NCNT
        IF (I_OBC_N(I) /= NLID(TEMP3(I)))&
             & CALL FATAL_ERROR("VEC_INT_DEAL BC TEST  : FAILED")
     END DO
  end if

  if(dbg_set(dbg_log)) then
     WRITE(IPT,*)  '!  BC DEAL TEST          :    PASSED    '
     !      WRITE(IPT,*)  '!  BC COLLECT TEST       :     '
  end if

!!$  TEMP1 = 0
!!$  DO I = 1, NPROCS_total
!!$     RECVER= I
!!$     if( .not. IOPROC .or. (RECVER .EQ. MYID)) &
!!$          CALL PCOLLECT(MYID,RECVER,NPROCS,BCMAP,TEMP3,TEMP1)
!!$  END DO

  TEMP1 = huge(TEMP1)
  IF(USE_MPI_IO_MODE)THEN
    RECVER = IOPROCID
    CALL PCOLLECT_IO(MYID,RECVER,NPROCS_TOTAL,BCMAP,TEMP3,TEMP1)
    IF(.not. IOPROC)THEN
      DO I = 1, NPROCS_FVCOM
        RECVER = I
        CALL PCOLLECT(MYID,RECVER,NPROCS_TOTAL,BCMAP,TEMP3,TEMP1)
      END DO
    ELSE
      IF(NPROCS_FVCOM > 1)THEN
        ALLOCATE(TEMP4(IOBCN_GL))
        CALL MPI_ALLREDUCE(TEMP1,TEMP4,IOBCN_GL,MPI_INTEGER,MPI_MIN,MPI_FVCOM_GROUP,IERR)
        TEMP1=TEMP4
        DEALLOCATE(TEMP4)
      END IF
    END IF
  ELSE
    DO I = 1, NPROCS_IO
      RECVER = I
      CALL PCOLLECT_IO(MYID,RECVER,NPROCS_TOTAL,BCMAP,TEMP3,TEMP1)
    END DO
  END IF
          	  	
  !      write(IPT,*)"temp1= ",temp1
  !      write(IPT,*)"I_OBC_GL= ",I_OBC_GL

  DO I = 1,IOBCN_GL
     IF( I_OBC_GL(I) .NE. TEMP1(I)) &
          & CALL FATAL_ERROR("VEC_INT_COLLECT BC TEST: FAILED")
  END DO

  if(dbg_set(dbg_log)) &
       & WRITE(IPT,*)  '!  BC COLLECT TEST       :    PASSED  '


  DEALLOCATE(TEMP1,TEMP3)

  ! Add the BOUNDARY CONDITION map to the Halo list
  CALL ADD_MAP2LIST(HALO_MAPS,BCMAP)



# endif

  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "END: SETUP_OBCMAP - parallel"   
  RETURN
END SUBROUTINE GENMAP_OBC

