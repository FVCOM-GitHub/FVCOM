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
   SUBROUTINE PRINT_VALS          

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_OBCS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (DATA_ASSIM)
   USE MOD_ASSIM
#  endif
   IMPLICIT NONE
   INTEGER :: I,K,IOUTTMP,ierr
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: UATMP,VATMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: T1TMP,S1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP,DTMP

!==============================================================================|
   

    IF(SERIAL)THEN

# if defined (TIDE_OUTPUT)
      if(iint >= OBC_TIDEOUT_INITIAL .and. mod(iint, OBC_TIDEOUT_INTERVAL)==0) then
         write(911,'(300f8.4)')    (EL(I_OBC_GL(I)),I=1,IOBCN_GL)
         write(912,'(I7,300f8.4)') iint,(EL(I_OBCNODE_GL(I)),I=1,IOBCNODE_GL)
         write(913,'(I7,300f8.4)') iint,(UA(I_OBCELL_GL(I)),I=1,IOBCELL_GL)
         write(913,'(I7,300f8.4)') iint,(VA(I_OBCELL_GL(I)),I=1,IOBCELL_GL)
         do K = 1,kbm1
            write(913,'(I7,300f8.4)') iint,(U(I_OBCELL_GL(I),k),I=1,IOBCELL_GL)
            write(913,'(I7,300f8.4)') iint,(V(I_OBCELL_GL(I),k),I=1,IOBCELL_GL)
         end do
      end if

      if(iint == OBC_TIDEOUT_INITIAL) then
         do I = 1,MGL
            write(914,*) EL(I)
         end do
         do I = 1,NGL
            write(915,'(300f8.4)') (U(I,K),V(I,K),K=1,kbm1)
         end do
         write(916,*) IOBCNODE_GL
         do I = 1,IOBCNODE_GL
            write(916,*) I_OBCNODE_GL(I)
         enddo
         write(917,*) IOBCELL_GL
         do I = 1,IOBCELL_GL
            write(917,*) I_OBCELL_GL(I)
         enddo
      end if
# endif

    END IF

#  if defined (MULTIPROCESSOR) && (TIDE_OUTPUT)
    IF(PAR)THEN
     ALLOCATE(UTMP(0:NGL,KB))
     ALLOCATE(VTMP(0:NGL,KB))
     ALLOCATE(T1TMP(0:MGL,KB))
     ALLOCATE(S1TMP(0:MGL,KB))
     ALLOCATE(ELTMP(0:MGL))
     ALLOCATE(DTMP(0:MGL))
     CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
     ALLOCATE(UATMP(NGL),VATMP(NGL))
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,U,  UTMP)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,V,  VTMP)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,UA, UATMP)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,VA, VATMP)
     CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,EL, ELTMP)
!     CALL GATHER(LBOUND(U,1),  UBOUND(U,1),  N,NGL,KB,MYID,NPROCS,EMAP,U,  UTMP)
!     CALL GATHER(LBOUND(V,1),  UBOUND(V,1),  N,NGL,KB,MYID,NPROCS,EMAP,V,  VTMP)
!     CALL GATHER(LBOUND(UA,1), UBOUND(UA,1), N,NGL, 1,MYID,NPROCS,EMAP,UA, UATMP)
!     CALL GATHER(LBOUND(VA,1), UBOUND(VA,1), N,NGL, 1,MYID,NPROCS,EMAP,VA, VATMP)
!!     CALL GATHER(LBOUND(T1,1),  UBOUND(T1,1),  M,MGL,KB,MYID,NPROCS,NMAP,T1,  T1TMP)
!!     CALL GATHER(LBOUND(S1,1),  UBOUND(S1,1),  M,MGL,KB,MYID,NPROCS,NMAP,S1,  S1TMP)
!     CALL GATHER(LBOUND(EL,1), UBOUND(EL,1), M,MGL, 1,MYID,NPROCS,NMAP,EL, ELTMP)

!      if(msr) then
!         write(900,'(i8,6f10.5)')iint,eltmp(5),eltmp(90),eltmp(1280)
!         write(901,'(i8,6f10.5)')iint,uatmp(5),vatmp(5),uatmp(90),vatmp(90),uatmp(2300),vatmp(2300)
!      end if

!!# if defined (TIDE_OUTPUT)
      if(msr .and. iint >= OBC_TIDEOUT_INITIAL .and. mod(iint, OBC_TIDEOUT_INTERVAL)==0) then
         write(911,'(300f8.4)')    (ELTMP(I_OBC_GL(I)),I=1,IOBCN_GL)
         write(912,'(I7,300f8.4)') iint,(ELTMP(I_OBCNODE_GL(I)),I=1,IOBCNODE_GL)
         write(913,'(I7,300f8.4)') iint,(UATMP(I_OBCELL_GL(I)),I=1,IOBCELL_GL)
         write(913,'(I7,300f8.4)') iint,(VATMP(I_OBCELL_GL(I)),I=1,IOBCELL_GL)
         do K = 1,kbm1
            write(913,'(I7,300f8.4)') iint,(UTMP(I_OBCELL_GL(I),k),I=1,IOBCELL_GL)
            write(913,'(I7,300f8.4)') iint,(VTMP(I_OBCELL_GL(I),k),I=1,IOBCELL_GL)
         end do
      end if

      if(msr .and. iint == OBC_TIDEOUT_INITIAL) then
         do I = 1,MGL
            write(914,*) ELTMP(I)
         end do
         do I = 1,NGL
            write(915,'(300f8.4)') (UTMP(I,K),VTMP(I,K),K=1,kbm1)
         end do
         write(916,*) IOBCNODE_GL
         do I = 1,IOBCNODE_GL
            write(916,*) I_OBCNODE_GL(I)
         enddo
         write(917,*) IOBCELL_GL
         do I = 1,IOBCELL_GL
            write(917,*) I_OBCELL_GL(I)
         enddo
      end if
!!# endif
   
   DEALLOCATE(T1TMP,S1TMP,ELTMP,DTMP)

   END IF
#  endif
 
    END
