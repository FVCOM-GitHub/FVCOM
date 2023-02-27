PROGRAM interp_test
  USE MOD_INTERP
  IMPLICIT NONE
 
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  INTEGER :: STATUS, TEST_CODE, I, J, n,nsq
  CHARACTER(len=80) :: FNAME

  REAL(SP), allocatable :: X(:),Y(:),Z(:),XZ(:),YZ(:),ZZ(:)

  TYPE(INTERP_WEIGHTS) :: MYWEIGHTS
  
  CALL INITIALIZE_CONTROL("INTERP TEST QUAD")

  ALLOCATE(X(10))
  ALLOCATE(Y(10))
  ALLOCATE(Z(10))
  
  ! X and Y are generated randomly from matlab
  x=(/ 0.814723686393179, 0.905791937075619, 0.126986816293506, 0.913375856139019,&
       & 0.632359246225410, 0.097540404999410, 0.278498218867048, 0.546881519204984,&
       & 0.957506835434298, 0.964888535199277/)
  
  y=(/0.157613081677548, 0.970592781760616, 0.957166948242946,&
       & 0.485375648722841, 0.800280468888800, 0.141886338627215, 0.421761282626275,&
       & 0.915735525189067, 0.792207329559554, 0.959492426392903/)

  ! z = x^2 - y^2
  Z=(/0.638932801654198,  -0.121591314734608,  -0.900042915296353,&
       & 0.598665934204567,  -0.240570612598109,  -0.010617602481588,  -0.100321321610442,&
       & -0.539492156045345,   0.289226886895523,   0.010384169053666 /)
  
  PRINT*,X
  PRINT*,Y
  PRINT*,Z
  
  
  n = 120
  nsq= n**2

! TEST GRID
  Allocate(XZ(0:nsq)); XZ= 0.0_SP
  ALLOCATE(YZ(0:nsq)); YZ= 0.0_SP
  ALLOCATE(ZZ(0:nsq)); ZZ= 0.0_SP
  
  
  DO J = 1,n
     DO I = 1,n
        XZ(I +n*(J-1)) = REAL(I-1)
        YZ(I +n*(J-1)) = REAL(J-1)
     END DO
  END DO
  
  XZ = XZ/100.0_SP
  YZ = YZ/100.0_SP
  
  CALL SETUP_INTERP_QUAD_A(X,Y,XZ,YZ,MYWEIGHTS,0.1_sp)
  
  
  CALL interp_quad_A(Z,MYWEIGHTS,ZZ)
  
  DO I =1,nsq
     write(11,*) XZ(I)
     write(12,*) YZ(I)
     write(13,*) ZZ(I)
  END DO


END PROGRAM interp_test

