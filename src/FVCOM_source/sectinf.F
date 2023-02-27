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
!  this subroutine is used to specify the output files used for the            !
!  graphics developed by the Ocean Ecosystem Dynamics Laboratory at            !
!  SMAST/UMASSD.                                                               !
!==============================================================================|

   SUBROUTINE SECTINF             

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   INTEGER :: I,K,J,N1,N2
   CHARACTER(LEN=80)   :: TEMPDIR
   REAL(SP), ALLOCATABLE,DIMENSION(:,:) :: FTEMP1,FTEMP2,FTEMP3,FTEMP4
   INTEGER , ALLOCATABLE,DIMENSION(:,:) :: NTEMP1,NTEMP2
   INTEGER , ALLOCATABLE,DIMENSION(:,:) :: NVTEMP,NBETEMP

!==============================================================================|
   
   TEMPDIR = TRIM(OUTDIR)//"/medm/"
!
!-----------------OUTPUT SURFACE NODE COORDINATES------------------------------|
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'xy_node.dat',STATUS='unknown') ; REWIND(1)

   IF(SERIAL)THEN
     DO I=1,M
       WRITE(1,*) VX(I)+VXMIN,VY(I)+VYMIN
     END DO
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     ALLOCATE(FTEMP1(MGL,1))
     ALLOCATE(FTEMP2(MGL,1))
     CALL GATHER(LBOUND(VX,1),UBOUND(VX,1),M,MGL,1,MYID,NPROCS,NMAP,VX,FTEMP1)
     CALL GATHER(LBOUND(VY,1),UBOUND(VY,1),M,MGL,1,MYID,NPROCS,NMAP,VY,FTEMP2)
     IF(MSR)THEN
       DO I=1,MGL
         WRITE(1,*) FTEMP1(I,1)+VXMIN,FTEMP2(I,1)+VYMIN
       END DO
     END IF
     DEALLOCATE(FTEMP1,FTEMP2)
   END IF
#  endif

   IF(MSR)CLOSE(1)

!
!------------------OUTPUT SURFACE ELEMENT COORDINATES--------------------------|
!

   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'xy_cell.dat',STATUS='unknown') ; REWIND(1)

   IF(SERIAL)THEN
     DO I=1,N
       WRITE(1,*) XC(I)+VXMIN,YC(I)+VYMIN
     END DO
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     ALLOCATE(FTEMP1(NGL,1))
     ALLOCATE(FTEMP2(NGL,1))
     CALL GATHER(LBOUND(XC,1),UBOUND(XC,1),N,NGL,1,MYID,NPROCS,EMAP,XC,FTEMP1)
     CALL GATHER(LBOUND(XC,1),UBOUND(XC,1),N,NGL,1,MYID,NPROCS,EMAP,YC,FTEMP2)
     IF(MSR)THEN
       DO I=1,NGL
         WRITE(1,*) FTEMP1(I,1)+VXMIN,FTEMP2(I,1)+VYMIN
       END DO
     END IF
     DEALLOCATE(FTEMP1,FTEMP2)
   END IF
#  endif

   IF(MSR)CLOSE(1)

!
!------------------OUTPUT EDGES AND VERTICES FOR EACH ELEMENT------------------|
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'mesh.inf',STATUS='unknown')  ; REWIND(1)
   IF(MSR)WRITE(1,*) 'nbe(i,j),j=1,3; nv(i,j),j=1,3'

   IF(SERIAL)THEN
     DO I=1,N
       WRITE(1,100) (NBE(I,J),J=1,3),(NV(I,J),J=1,3)
     END DO
   END IF
  
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     ALLOCATE(NTEMP1(NGL,3))
     ALLOCATE(NTEMP2(NGL,3))
     ALLOCATE(NBETEMP(LBOUND(NBE,1):UBOUND(NBE,1),3))
     ALLOCATE( NVTEMP(LBOUND(NV ,1):UBOUND(NV ,1),3))

!    TRANSFORM NBE AND NV ARRAYS TO GLOBAL INDEXING
     DO J=1,3
     DO I=1,N
       IF(NBE(I,J) > N) THEN
         NBETEMP(I,J) = HE_LST(NBE(I,J)-N)
       ELSE
         NBETEMP(I,J) = EGID(NBE(I,J))
       END IF
       NVTEMP(I,J)  = NGID(NV(I,J))
     END DO
     END DO

     CALL IGATHER(LBOUND(NBE,1),UBOUND(NBE,1),N,NGL,3,MYID,NPROCS,EMAP,NBETEMP,NTEMP1)
     CALL IGATHER(LBOUND(NV ,1),UBOUND(NV ,1),N,NGL,3,MYID,NPROCS,EMAP,NVTEMP ,NTEMP2)
     IF(MSR)THEN
       DO I=1,NGL
       WRITE(1,100) (NTEMP1(I,J),J=1,3),(NTEMP2(I,J),J=1,3)
       END DO
     END IF
     DEALLOCATE(NTEMP1,NTEMP2,NBETEMP,NVTEMP)
   END IF
#  endif

   IF(MSR)CLOSE(1)

!
!------------------SHAPE FACTORS-----------------------------------------------|
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'shape.inf',STATUS='unknown') ; REWIND(1)
   IF(MSR)WRITE(1,*) 'au(i,j),j=1,4; av(i,j),j=1,4'

   IF(SERIAL)THEN
     DO I=1,N
       WRITE(1,200) (A1U(I,J),J=1,4),(A2U(I,J),J=1,4)
     END DO
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     ALLOCATE(FTEMP1(NGL,4))
     ALLOCATE(FTEMP2(NGL,4))
     CALL GATHER(LBOUND(A1U,1),UBOUND(A1U,1),N,NGL,4,MYID,NPROCS,EMAP,A1U,FTEMP1)
     CALL GATHER(LBOUND(A2U,1),UBOUND(A2U,1),N,NGL,4,MYID,NPROCS,EMAP,A2U,FTEMP2)
     IF(MSR)THEN
       DO I=1,NGL
       WRITE(1,200) (FTEMP1(I,J),J=1,4),(FTEMP2(I,J),J=1,4)
       END DO
     END IF
     DEALLOCATE(FTEMP1,FTEMP2)
   END IF
#  endif


   IF(MSR)CLOSE(1)

!
!-------------------LINEAR FUNCTIONS-------------------------------------------|
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'awxcof.inf',STATUS='unknown') ; REWIND(1)
   IF(MSR)WRITE(1,*) 'aw0(i,j),awx(i,j),awy(i,j),j=1,3'

   IF(SERIAL)THEN
     DO I=1,N
       WRITE(1,300) (AW0(I,J),AWX(I,J),AWY(I,J),J=1,3)
     END DO
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     ALLOCATE(FTEMP1(NGL,3))
     ALLOCATE(FTEMP2(NGL,3))
     ALLOCATE(FTEMP3(NGL,3))
     CALL GATHER(LBOUND(AW0,1),UBOUND(AW0,1),N,NGL,3,MYID,NPROCS,EMAP,AW0,FTEMP1)
     CALL GATHER(LBOUND(AWX,1),UBOUND(AWX,1),N,NGL,3,MYID,NPROCS,EMAP,AWX,FTEMP2)
     CALL GATHER(LBOUND(AWY,1),UBOUND(AWY,1),N,NGL,3,MYID,NPROCS,EMAP,AWY,FTEMP3)
     IF(MSR)THEN
       DO I=1,NGL
       WRITE(1,300) (FTEMP1(I,J),FTEMP2(I,J),FTEMP3(I,J),J=1,3)
       END DO
     END IF
     DEALLOCATE(FTEMP1,FTEMP2,FTEMP3)
   END IF
#  endif


   IF(MSR)CLOSE(1)

!
!------------------OUTPUT DEPTH AT NODE POINTS---------------------------------!
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'depth.xy',STATUS='unknown') ; REWIND(1)
   IF(MSR)WRITE(1,*) 'scat2d'
   IF(MSR)WRITE(1,500) MGL,1

   IF(SERIAL)THEN
     DO I=1,M
       WRITE(1,'(3E20.10)') VX(I)+VXMIN,VY(I)+VYMIN,H(I)
     END DO
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     ALLOCATE(FTEMP1(MGL,1))
     ALLOCATE(FTEMP2(MGL,1))
     ALLOCATE(FTEMP3(MGL,1))
     CALL GATHER(LBOUND(VX,1),UBOUND(VX,1),M,MGL,1,MYID,NPROCS,NMAP,VX,FTEMP1)
     CALL GATHER(LBOUND(VY,1),UBOUND(VY,1),M,MGL,1,MYID,NPROCS,NMAP,VY,FTEMP2)
     CALL GATHER(LBOUND(H ,1),UBOUND(H ,1),M,MGL,1,MYID,NPROCS,NMAP,H ,FTEMP3)
     IF(MSR)THEN
       DO I=1,MGL
       WRITE(1,'(3E20.10)') FTEMP1(I,1)+VXMIN,FTEMP2(I,1)+VYMIN,FTEMP3(I,1)
       END DO
     END IF
     DEALLOCATE(FTEMP1,FTEMP2,FTEMP3)
   END IF
#  endif

   IF(MSR)CLOSE(1)

!
!------------------OUTPUT SIGMA DISTRIBUTION-----------------------------------!
!
   IF(MSR)OPEN(1,FILE=TRIM(TEMPDIR)//'sigma.dat',STATUS='unknown') ; REWIND(1)
  
     IF(SERIAL)THEN
       DO I=1,M
         WRITE(1,'(I10,100E20.10)') I,(Z(I,K),K=1,KB)
       END DO
     END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     ALLOCATE(FTEMP4(MGL,KB))
     CALL GATHER(LBOUND(Z ,1),UBOUND(Z ,1),M,MGL,KB,MYID,NPROCS,NMAP,Z ,FTEMP4)
     IF(MSR)THEN
       DO I=1,MGL
       WRITE(1,'(I10,100E20.10)') I,(FTEMP4(I,K),K=1,KB)
       END DO
     END IF
     DEALLOCATE(FTEMP4)
   END IF
#  endif

     IF(MSR)CLOSE(1)
!
!-------------------------FORMATTING-------------------------------------------!
!

100 FORMAT(6I10)
200 FORMAT(8E18.8)
300 FORMAT(9E18.8)
500 FORMAT('xyd ',i10,' depth ',i3,' h')

   RETURN
   END SUBROUTINE SECTINF
!==============================================================================|
