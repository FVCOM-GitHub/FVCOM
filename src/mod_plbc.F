
MODULE MOD_PERIODIC_LBC
# if defined (PLBC)
   USE ALL_VARS
   USE MOD_PREC
   USE MOD_PAR
   USE VARS_WAVE

   IMPLICIT NONE
   SAVE
   
   INTEGER              :: nplbccell_GL, nplbccell
   INTEGER              :: nplbcnode_GL, nplbcnode
   INTEGER, ALLOCATABLE :: I_PLBCNODE_GL(:,:),I_PLBCNODE_N(:,:)
   INTEGER, ALLOCATABLE :: I_PLBCCELL_GL(:,:),I_PLBCCELL_N(:,:)

   CONTAINS

!========================================================================
     SUBROUTINE FIND_NODE_CELL
     USE ALL_VARS
     IMPLICIT NONE
     INTEGER :: I,J,I1,NCNT
     INTEGER, ALLOCATABLE :: TEMP(:)

!----Read in Node Data----------------------------------
!
     REWIND(INTNODE1)
     READ(INTNODE1,*) nplbcnode_GL

     nplbcnode = 0
  IF (nplbcnode_GL > 0) THEN
     ALLOCATE(I_PLBCNODE_GL(nplbcnode_GL,5))
     DO I=1,nplbcnode_GL
        READ(INTNODE1,*) (I_PLBCNODE_GL(I,J),J=1,5)
     ENDDO
  END IF
     CLOSE(INTNODE1)
!
!---Map to Local Domain----------------------------------------
!
     IF(SERIAL) THEN
        nplbcnode    = nplbcnode_GL
        ALLOCATE(I_PLBCNODE_N(nplbcnode,5))
        I_PLBCNODE_N(:,:) = I_PLBCNODE_GL(:,:)
     ENDIF

#    if defined (MULTIPROCESSOR)
!     IF(PAR)THEN
!       ALLOCATE(TEMP(nplbcnode_GL))
!       NCNT = 0
!       DO I=1,nplbcnode_GL
!         I1=NLID(I_PLBCNODE_GL1(I))
!         IF(I1 /= 0)THEN
!           NCNT = NCNT + 1
!           TEMP(NCNT) = I1
!         END IF
!       END DO
!       nplbcnode = NCNT
!       IF(nplbcnode > 0)THEN       
!         ALLOCATE(I_PLBCNODE_N1(nplbcnode))
!         I_PLBCNODE_N1(1:nplbcnode) = TEMP(1:nplbcnode)
!       END IF
!!
!       NCNT=0
!       DO I=1,nplbcnode_GL
!         I1=NLID(I_PLBCNODE_GL2(I))
!         IF(I1 /= 0)THEN
!           NCNT = NCNT + 1
!           TEMP(NCNT) = I1
!         END IF
!       END DO
!       nplbcnode = NCNT
!       IF(nplbcnode > 0)THEN
!         ALLOCATE(I_PLBCNODE_N2(nplbcnode))
!         I_PLBCNODE_N2(1:nplbcnode) = TEMP(1:nplbcnode)
!       END IF
!       DEALLOCATE(TEMP)
!     END IF
#    endif

!----Read in Cell Data----------------------------------
     REWIND(INTCELL1)
     READ(INTCELL1,*) nplbccell_GL

     nplbccell = 0
  IF (nplbccell_GL > 0) THEN

     ALLOCATE(I_PLBCCELL_GL(nplbccell_GL,4))
     DO I=1,nplbccell_GL
        READ(INTCELL1,*)(I_PLBCCELL_GL(I,J),J=1,4)
     ENDDO
  END IF
     CLOSE(INTCELL1)

!
!---Map to Local Domain----------------------------------------
!
     IF(SERIAL) THEN
        nplbccell = nplbccell_GL
        ALLOCATE(I_PLBCCELL_N(nplbccell,4))
        I_PLBCCELL_N(:,:) = I_PLBCCELL_GL(:,:)
     ENDIF

#    if defined (MULTIPROCESSOR)
!     IF(PAR)THEN
!       ALLOCATE(TEMP(nplbccell_GL))
!       NCNT = 0
!       DO I=1,nplbccell_GL
!         I1=ELID(I_PLBCCELL_GL1(I))
!         IF(I1 /= 0)THEN
!           NCNT = NCNT + 1
!           TEMP(NCNT) = I1
!         END IF
!       END DO
!       nplbccell = NCNT
!       IF(nplbccell > 0)THEN                        
!         ALLOCATE(I_PLBCCELL_N1(nplbccell))
!         I_PLBCCELL_N1(1:nplbccell) = TEMP(1:nplbccell)
!       END IF
!
!       NCNT = 0
!       DO I=1,nplbccell_GL
!         I1=ELID(I_PLBCCELL_GL2(I))
!         IF(I1 /= 0)THEN
!           NCNT = NCNT + 1
!           TEMP(NCNT) = I1
!         END IF
!       END DO
!       nplbccell = NCNT
!       IF(nplbccell > 0)THEN
!         ALLOCATE(I_PLBCCELL_N2(nplbccell))
!         I_PLBCCELL_N2(1:nplbccell) = TEMP(1:nplbccell)
!       END IF
!
!       DEALLOCATE(TEMP)
!     ENDIF
#    endif
END SUBROUTINE FIND_NODE_CELL


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine replace_vel_2D
  implicit none
  integer :: i,j1,j2,j3,j4
  do i=1,nplbccell
     j1=I_PLBCCELL_N(i,1)
     j2=I_PLBCCELL_N(i,2)
     j3=I_PLBCCELL_N(i,3)
     j4=I_PLBCCELL_N(i,4)
     UAF(j1)=UAF(j2)
     UAF(j3)=UAF(j2)
     UAF(j4)=UAF(j2)
     VAF(j1)=VAF(j2)
     VAF(j3)=VAF(j2)
     VAF(j4)=VAF(j2)
  end do
  !VAF=0.0_SP
  !UAF=0.0_SP
  end subroutine replace_vel_2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine replace_vel_3D
  implicit none
  integer :: i,j1,j2,j3,j4
  do i=1,nplbccell
     j1=I_PLBCCELL_N(i,1)
     j2=I_PLBCCELL_N(i,2)
     j3=I_PLBCCELL_N(i,3)
     j4=I_PLBCCELL_N(i,4)

     UF(j1,1:KBM1)=UF(j2,1:KBM1)
     UF(j3,1:KBM1)=UF(j2,1:KBM1)
     UF(j4,1:KBM1)=UF(j2,1:KBM1)
     VF(j1,1:KBM1)=VF(j2,1:KBM1)
     VF(j3,1:KBM1)=VF(j2,1:KBM1)
     VF(j4,1:KBM1)=VF(j2,1:KBM1)
  end do
   !VF=0.0_SP
   !UF=0.0_SP
  end subroutine replace_vel_3D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine replace_ele
  implicit none
  integer :: i,j1,j2,j3,j4,j5

  do i=1,nplbcnode
     j1=I_PLBCNODE_N(i,1)
     j2=I_PLBCNODE_N(i,2)
     j3=I_PLBCNODE_N(i,3)
     j4=I_PLBCNODE_N(i,4)
     j5=I_PLBCNODE_N(i,5)
     ELF(j1)=ELF(j3)
     ELF(j5)=ELF(j3)
     ELF(j4)=ELF(j2)
  end do
     ELF(198)=ELF(149)
     ELF(99)=ELF(149)
     ELF(51)=ELF(100)
     ELF(150)=ELF(100)
     !ELF=0.0
  end subroutine replace_ele
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine replace_vel_w
  implicit none
 integer :: i,j1,j2,j3,j4,j5

  do i=1,nplbcnode
     j1=I_PLBCNODE_N(i,1)
     j2=I_PLBCNODE_N(i,2)
     j3=I_PLBCNODE_N(i,3)
     j4=I_PLBCNODE_N(i,4)
     j5=I_PLBCNODE_N(i,5)
     wts(j1,2:KB)=wts(j3,2:KB)
     wts(j5,2:KB)=wts(j3,2:KB)
     wts(j2,2:KB)=wts(j4,2:KB)
  end do
     wts(198,2:KB)=wts(149,2:KB)
     wts(99,2:KB)=wts(149,2:KB)
     wts(51,2:KB)=wts(100,2:KB)
     wts(150,2:KB)=wts(100,2:KB)
  end subroutine replace_vel_w
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine replace_q2
  implicit none
  integer :: i,j1,j2,j3,j4,j5

  do i=1,nplbcnode
     j1=I_PLBCNODE_N(i,1)
     j2=I_PLBCNODE_N(i,2)
     j3=I_PLBCNODE_N(i,3)
     j4=I_PLBCNODE_N(i,4)
     j5=I_PLBCNODE_N(i,5)
     q2f(j1,:)=q2f(j3,:)
     q2f(j5,:)=q2f(j3,:)
     q2f(j2,:)=q2f(j4,:)
  end do
  q2f(198,:)=q2f(149,:)
  q2f(99,:)=q2f(149,:)
  q2f(51,:)=q2f(100,:)
  q2f(150,:)=q2f(100,:)
  end subroutine replace_q2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine replace_q2l
  implicit none
  integer :: i,j1,j2,j3,j4,j5
  
  do i=1,nplbcnode
     j1=I_PLBCNODE_N(i,1)
     j2=I_PLBCNODE_N(i,2)
     j3=I_PLBCNODE_N(i,3)
     j4=I_PLBCNODE_N(i,4)
     j5=I_PLBCNODE_N(i,5)
     q2lf(j1,:)=q2lf(j3,:)
     q2lf(j5,:)=q2lf(j3,:)
     q2lf(j2,:)=q2lf(j4,:)
  end do
  q2lf(198,:)=q2lf(149,:)
  q2lf(99,:)=q2lf(149,:)
  q2lf(51,:)=q2lf(100,:)
  q2lf(150,:)=q2lf(100,:)
  end subroutine replace_q2l


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine replace_ac2(ID,ISS)
  implicit none
  integer :: ID,ISS
  integer :: i,j1,j2,j3,j4,j5
  do i=1,nplbcnode
     j1=I_PLBCNODE_N(i,1)
     j2=I_PLBCNODE_N(i,2)
     j3=I_PLBCNODE_N(i,3)
     j4=I_PLBCNODE_N(i,4)
     j5=I_PLBCNODE_N(i,5)
     AC2(ID,ISS,j1)=AC2(ID,ISS,j3)
     AC2(ID,ISS,j5)=AC2(ID,ISS,j3)
     AC2(ID,ISS,j2)=AC2(ID,ISS,j4)
  end do
     AC2(ID,ISS,198)=AC2(ID,ISS,149)
     AC2(ID,ISS,99)=AC2(ID,ISS,149)
  end subroutine replace_ac2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine replace_N32(N32,ID,ISS)
   USE ALL_VARS, ONLY : MT
   USE SWCOMM3, ONLY : MDC,MSC
  IMPLICIT NONE
  real,dimension(MDC,MSC,0:MT) :: N32
  integer :: i,j1,j2,j3,j4,j5,ID,ISS
  do i=1,nplbcnode
     j1=I_PLBCNODE_N(i,1)
     j2=I_PLBCNODE_N(i,2)
     j3=I_PLBCNODE_N(i,3)
     j4=I_PLBCNODE_N(i,4)
     j5=I_PLBCNODE_N(i,5)
     N32(ID,ISS,j1)=N32(ID,ISS,j3)
     N32(ID,ISS,j5)=N32(ID,ISS,j3)
     N32(ID,ISS,j2)=N32(ID,ISS,j4)
  end do
     N32(ID,ISS,198)=N32(ID,ISS,149)
     N32(ID,ISS,99)=N32(ID,ISS,149)

  end subroutine replace_N32
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine replace_node(x)
   USE ALL_VARS, ONLY : MT
  IMPLICIT NONE
  real,dimension(0:MT) :: x
  integer :: i,j1,j2,j3,j4,j5
  do i=1,nplbcnode
     j1=I_PLBCNODE_N(i,1)
     j2=I_PLBCNODE_N(i,2)
     j3=I_PLBCNODE_N(i,3)
     j4=I_PLBCNODE_N(i,4)
     j5=I_PLBCNODE_N(i,5)
     x(j1)=x(j3)
     x(j5)=x(j3)
     x(j2)=x(j4)
  end do
     x(198)=x(149)
     x(99)=x(149)
     x(51)=x(100)
     x(150)=x(100)

  end subroutine replace_node

# endif
END MODULE MOD_PERIODIC_LBC

