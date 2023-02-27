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
SUBROUTINE GENMAP_LSF
!==============================================================================|
!      OBC Node Number
!      TS OBC Type
!      OBC Node list and types
!==============================================================================|
  USE MOD_PAR
  USE LIMS
  USE CONTROL
  USE ALL_VARS
  IMPLICIT NONE
  integer :: SENDER,RECVER, ierr, I, NCNT, NSZE, I1

  INTEGER, POINTER :: TEMP1(:),TEMP2(:)
  INTEGER, POINTER :: TEMP3(:),TEMP4(:)

  REAL(SP), POINTER :: FTMP1(:),FTMP2(:)
  REAL(SP), POINTER :: FTMP3(:),FTMP4(:)


  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "START: SETUP_LSFMAP"

  IF(SERIAL) THEN
     NOBCLSF    = NOBCLSF_GL

     ALLOCATE(IBCLSF(NOBCLSF))
     IBCLSF = IBCLSF_GL
     DEALLOCATE(IBCLSF_GL)

     ALLOCATE(RBC_GEO(NOBCLSF))
     RBC_GEO = RBC_GEO_GL
     DEALLOCATE(RBC_GEO_GL)

     ALLOCATE(RBC_WDF(NOBCLSF))
     RBC_WDF = RBC_WDF_GL
     DEALLOCATE(RBC_WDF_GL)

     ALLOCATE(IBCLSF_OUTPUT(NOBCLSF))
     IBCLSF_OUTPUT = IBCLSF

     if (dbg_set(dbg_sbr)) &
          & write(IPT,*) "END: GENMAP_LSF - serial"
     return
  END IF

# if defined (MULTIPROCESSOR)

  ALLOCATE(TEMP1(NOBCLSF_GL))
  ALLOCATE(FTMP2(NOBCLSF_GL))
  ALLOCATE(FTMP3(NOBCLSF_GL))
  ALLOCATE(TEMP4(NOBCLSF_GL))
  
  NOBCLSF = 0
  NCNT = 0
  IF(.not. IOPROC) THEN
     FTMP2=0.0_SP
     FTMP3=0.0_SP

     !!SET UP LOCAL OPEN BOUNDARY NODES
     DO I=1,NOBCLSF_GL

        I1 = NLID( IBCLSF_GL(I) ) ! MUST NOT INCLUD HALO
        IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           FTMP2(NCNT) = RBC_GEO_GL(I)
           FTMP3(NCNT) = RBC_WDF_GL(I)
           TEMP4(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF BOUNDARY NODES
     NOBCLSF = NCNT
     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(IBCLSF(NCNT),RBC_GEO(NCNT),RBC_WDF(NCNT),IBCLSF_OUTPUT(NCNT))
     IBCLSF  = TEMP1(1:NCNT)
     RBC_GEO = FTMP2(1:NCNT)
     RBC_WDF = FTMP3(1:NCNT)
     IBCLSF_OUTPUT = NGID(IBCLSF)
  END IF

  !==============================================================================|
  !   SET UP ELEMENT MAPPING FOR GLOBAL 2 LOCAL TRANSFER OF BC'S                 | 
  !   BOUNDARY MAP :: BCMAP(NPROCS)                                              |
  !     BCMAP(1-->NPROCS)%NSIZE  :: NUMBER OF BOUNDARY NODES IN EACH DOM         |
  !     BCMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM |
  !==============================================================================|

  TEMP2 => TEMP4(1:NCNT)
  LSFMAP => MAKE_MAP(MYID,NPROCS_TOTAL,NOBCLSF_GL,NCNT,TEMP2)
  NULLIFY(TEMP2)

!!$  ALLOCATE(LSFMAP(NPROCS)); LSFMAP(:)%NSIZE=0
!!$
!!$  !--Determine Number of Elements for Each Processor
!!$  DO I=1,NPROCS    
!!$     IF(MYID == I) LSFMAP(I)%NSIZE = NCNT
!!$     SENDER = I - 1
!!$     CALL MPI_BCAST(LSFMAP(I)%NSIZE,1,MPI_INTEGER,SENDER,MPI_COMM_FVCOM,IERR)
!!$  END DO
!!$
!!$  !--Allocate Mapping Array for Each Processor
!!$  DO I=1,NPROCS
!!$     ALLOCATE(LSFMAP(I)%LOC_2_GL(0:LSFMAP(I)%NSIZE))
!!$     LSFMAP(I)%LOC_2_GL(0) = 0
!!$  END DO
!!$
!!$  !--Construct Mapping Array for Each Processor 
!!$  DO I=1,NPROCS
!!$     NSZE = LSFMAP(I)%NSIZE
!!$     if(myid == I) LSFMAP(I)%LOC_2_GL(1:NSZE) =  TEMP4(1:NSZE)
!!$     SENDER = I -1
!!$     CALL MPI_BCAST(LSFMAP(I)%LOC_2_GL(1:NSZE),NSZE,MPI_INTEGER,Sender,MPI_COMM_FVCOM,IERR)
!!$
!!$  END DO

  DEALLOCATE(TEMP1,FTMP2,FTMP3,TEMP4)

  ! TEST THE MAP

  ALLOCATE(TEMP1(NOBCLSF_GL))
  TEMP1 = IBCLSF_GL


  ALLOCATE(TEMP3(NCNT))
  TEMP3 = 0 ! THIS ASSIGNMENT MAY CAUSE A PROBLEM ON THE IOPROC WHERE NCNT == 0


  SENDER =MSRID
  if (USE_MPI_IO_MODE) SENDER = IOPROCID ! TEST DEAL FROM IOPROC
!!$  CALL PDEAL(MYID,SENDER,NPROCS,LSFMAP,TEMP1,TEMP3)
  CALL PDEAL_IO(MYID,SENDER,NPROCS_TOTAL,LSFMAP,TEMP1,TEMP3)

  if (.not. IOPROC) then
     DO I =1, NCNT
        IF (IBCLSF(I) .NE. NLID(TEMP3(I)) ) CALL FATAL_ERROR&
             & ("VEC_INT_DEAL LSF TEST : FAILED")
     END DO
  end if

  if(dbg_set(dbg_log)) then
     WRITE(IPT,*)  '!  LSF DEAL TEST         :    PASSED    '
     !      WRITE(IPT,*)  '!  BC COLLECT TEST       :     '
  end if


!!$  TEMP1 = 0
!!$  DO I = 1, NPROCS_total
!!$     RECVER= I
!!$     if( .not. IOPROC .or. (RECVER .EQ. MYID)) &
!!$          CALL PCOLLECT(MYID,RECVER,NPROCS,LSFMAP,TEMP3,TEMP1)
!!$  END DO
   
   TEMP1 = huge(TEMP1)
   IF(USE_MPI_IO_MODE)THEN
     RECVER = IOPROCID
     CALL PCOLLECT_IO(MYID,RECVER,NPROCS_TOTAL,LSFMAP,TEMP3,TEMP1)
     IF(.not. IOPROC)THEN
       DO I = 1, NPROCS_FVCOM
         RECVER = I
         CALL PCOLLECT(MYID,RECVER,NPROCS_TOTAL,LSFMAP,TEMP3,TEMP1)
       END DO
     ELSE
       IF(NPROCS_FVCOM > 1)THEN
         ALLOCATE(TEMP4(NOBCLSF_GL))
         CALL MPI_ALLREDUCE(TEMP1,TEMP4,NOBCLSF_GL,MPI_INTEGER,MPI_MIN,MPI_FVCOM_GROUP,IERR)
         TEMP1 = TEMP4
         DEALLOCATE(TEMP4)
       END IF
     END IF
   ELSE
     DO I = 1, NPROCS_IO
       RECVER = I
       CALL PCOLLECT_IO(MYID,RECVER,NPROCS_TOTAL,LSFMAP,TEMP3,TEMP1)
     END DO
   END IF        	  	 

   DO I = 1, NPROCS_IO
     RECVER = I
     CALL PCOLLECT_IO(MYID,RECVER,NPROCS_TOTAL,LSFMAP,TEMP3,TEMP1)
   END DO

  !      write(IPT,*)"temp1= ",temp1
  !      write(IPT,*)"I_OBC_GL= ",I_OBC_GL

  DO I = 1,NOBCLSF_GL
     IF( IBCLSF_GL(I) .NE. TEMP1(I)) &
          & CALL FATAL_ERROR("VEC_INT_COLLECT LSF TEST: FAILED")
  END DO

  if(dbg_set(dbg_log)) &
       & WRITE(IPT,*)  '!  LSF COLLECT TEST      :    PASSED  '


  DEALLOCATE(TEMP1,TEMP3)


  ! Add it to the Halo list only
  CALL ADD_MAP2LIST(HALO_MAPS,LSFMAP)


# endif

  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "END: SETUP_LSFMAP - parallel"   
  RETURN
END SUBROUTINE GENMAP_LSF

