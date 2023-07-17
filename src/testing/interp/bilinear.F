PROGRAM interp_test
  USE MOD_INTERP
  ! SEE PARAMETER SETTINGS IN MOD_SST2GRD.F90
  IMPLICIT NONE
 
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  INTEGER :: STATUS, TEST_CODE, I, J, n,nsq,nr
  CHARACTER(len=80) :: FNAME

  REAL(SP), allocatable :: X(:,:),Y(:,:),XZ(:),YZ(:),Z(:,:),ZZ(:)
  INTEGER,ALLOCATABLE :: Land_Mask(:,:)

  real(SP) :: SY,SX

  TYPE(INTERP_WEIGHTS),POINTER :: MYWEIGHTS
  
  CALL INITIALIZE_CONTROL("INTERP TEST")

  ALLOCATE(Land_Mask(2,2))
!  ALLOCATE(X(3,3))
!  ALLOCATE(Y(3,3))
  ALLOCATE(X(2,2))
  ALLOCATE(Y(2,2))
  ALLOCATE(Z(2,2))
  
  
!!$	x(:,1)=(/ 0, 5,10/)
!!$	x(:,2)=(/ 0, 5,10/)
!!$	x(:,3)=(/ 0, 5,10/)
!!$
!!$	y(:,1)=(/ 0, 0, 0/)
!!$	y(:,2)=(/ 5, 5, 5/)
!!$	y(:,3)=(/10,10,10/)
  
  

  Land_Mask(:,1)=(/ 1, 0/)
  Land_Mask(:,2)=(/ 0, 0/)


  x(:,1)=(/ 150.0, 170.0/)
  x(:,2)=(/ 150.0, 190.0/)
 
  y(:,1)=(/   15.0,  15.0/)
  y(:,2)=(/ -10.0, -10.0/)

! ROTATED
!!$  x(:,1)=(/ .0, 0./)
!!$  x(:,2)=(/ 1.0, 1./)
!!$  
!!$  y(:,1)=(/  .0, 1./)
!!$  y(:,2)=(/ .0,1./)
  
  
  Z(:,1)=(/ 1., .5/)
  Z(:,2)=(/ 0.,  .75/)
  
  PRINT*,X
  PRINT*,Y
  PRINT*,Z
  PRINT*,LAND_MASK
  
  
! TEST POINT
!!$  Allocate(XZ(1))
!!$  ALLOCATE(YZ(1))
!!$  do while(.true.)
!!$     
!!$     ALLOCATE(MYWEIGHTS)
!!$
!!$     print *,&
!!$          & 'input position for testing coefficients(-9999 -9999 to continue)'
!!$     read(*,*)xz(1),yz(1)
!!$     
!!$     if(xz(1).eq.-9999) exit
!!$     
!!$     print *,'input',xz(1),yz(1)
!!$!     CALL SETUP_INTERP_BILINEAR_A(X,Y,XZ,YZ,MYWEIGHTS,Land_Mask)
!!$     CALL SETUP_INTERP_BILINEAR_A(X,Y,XZ,YZ,MYWEIGHTS)
!!$     
!!$     CALL PRINT_WEIGHTS(MYWEIGHTS,"TADA")
!!$
!!$     CALL KILL_WEIGHTS(MYWEIGHTS)
!!$     
!!$     PRINT*,"done"
!!$     
!!$  enddo
!!$  DEAllocate(XZ)
!!$  DEALLOCATE(YZ)
  

!  CALL PSHUTDOWN

  n = 61
  nr = 121
  nsq= nr**2

! TEST GRID
  Allocate(XZ(nsq)); XZ= 0.0_SP
  ALLOCATE(YZ(nsq)); YZ= 0.0_SP
  ALLOCATE(ZZ(nsq)); ZZ= 0.0_SP
  
  
  DO J = 1,nr
     DO I = 1,nr
        XZ(I +nr*(J-1)) = REAL(I-n)
        YZ(I +nr*(J-1)) = REAL(J-n)
     END DO
  END DO
  
  XZ = XZ/real(n-1)
  YZ = YZ/real(n-1)

  Print*,"min/max Xz",Minval(xz),maxval(Xz)
  Print*,"min/max YZ",Minval(Yz),maxval(Yz)

  SX = abs(sum(X(1,:)) - sum(X(2,:)))/3.5_SP
  SY =   abs(sum(Y(:,1)) - sum(Y(:,2)))/3.5_SP

  PRINT*,"SX,SY",SX,SY

  XZ = XZ * SX

  YZ = YZ * SY

  Print*,"min/max Xz",Minval(xz),maxval(Xz)
  Print*,"min/max YZ",Minval(Yz),maxval(Yz)

  XZ = sum(sum(X,1))/4.0_SP + XZ

  YZ = sum(sum(Y,1))/4.0_SP + YZ

  Print*,"min/max Xz",Minval(xz),maxval(Xz)
  Print*,"min/max YZ",Minval(Yz),maxval(Yz)
  
!!$  WHERE(YZ>90.0_SP)
!!$     YZ = 90.0_SP
!!$  ELSEWHERE(YZ<-90.0_SP)
!!$     YZ = -90.0_SP
!!$  END WHERE


  ALLOCATE(MYWEIGHTS)
  CALL SETUP_INTERP_BILINEAR_A(X,Y,XZ,YZ,MYWEIGHTS,LAND_MASK)
!  CALL SETUP_INTERP_BILINEAR_A(X,Y,XZ,YZ,MYWEIGHTS)
  
  
  CALL interp_bilinear_A(Z,MYWEIGHTS,ZZ)
  
  DO I =1,nsq
     write(11,*) XZ(I)
     write(12,*) YZ(I)
     write(13,*) ZZ(I)
  END DO

  write(14,*) X
  write(15,*) Y

  CALL KILL_WEIGHTS(MYWEIGHTS)

END PROGRAM interp_test

