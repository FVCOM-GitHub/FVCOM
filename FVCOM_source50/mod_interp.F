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

MODULE MOD_INTERP
  USE MOD_UTILS
  USE MOD_PREC
  USE MOD_CLOCK
  USE MOD_SPHERICAL
  IMPLICIT NONE  



  !/===========================================================================/
  !
  !  LAND MASK IN BILINEAR INTERP!
  !
  !  MASK == 0  OCEAN
  !  MASK == 1  LAND
  !
  !/===========================================================================/

  ! SEARCH RADIUS
  REAL(SP) :: search=80000.0_SP

  ! SLOPE TOLERANCE FOR CURVALINEAR INTERP PROGRAM
  REAL(SP), PARAMETER :: TB_TOL=1000.0_sp
  REAL(SP), PARAMETER :: LR_TOL=1000.0_sp

  REAL(SP), PARAMETER :: small= 1E-6

  ! NUMBER OF POINTS PER QUADRENT
  !  INTEGER, parameter :: nppq=2
  !  INTEGER, parameter :: quad=4

  TYPE R2PTS
     REAL(SP) :: WEIGHTS(4,2)
     INTEGER :: QUAD_NUM(4,2)
  END TYPE R2PTS

  TYPE INTERP_WEIGHTS
     INTEGER :: Nin
     INTEGER :: Nout
     INTEGER, ALLOCATABLE :: INDEX(:)
     TYPE(R2PTS), ALLOCATABLE:: PTW(:)

     ! HANDLE LAND MASK TO FIND NEAREST NEIGHBOR IN MESH WITH DATA
     ! THIS IS A HACK FOR GETTING VALUES ON AN INLAND RIVER...
     INTEGER :: N_TYPE
     INTEGER, POINTER :: N_INDEX(:,:)
     INTEGER, POINTER :: N_CNT(:)
     INTEGER, POINTER :: N_ORDER(:)
     INTEGER, POINTER :: N_FOUND(:)


  END TYPE INTERP_WEIGHTS

  INTEGER, POINTER :: N_FOUND(:)
  
  INTEGER, PARAMETER :: TYPE_NODE = 1
  INTEGER, PARAMETER :: TYPE_ELEMENT = 2



  TYPE(WATCH) :: MIN_W

  TYPE(WATCH) :: BOX_W

  TYPE(WATCH) :: COND_W

  TYPE(WATCH) :: CASE_W

  TYPE(WATCH) :: TOT_W


CONTAINS

  SUBROUTINE KILL_WEIGHTS(MYW)
    IMPLICIT NONE
    TYPE(INTERP_WEIGHTS),POINTER :: MYW


    IF(.not.ASSOCIATED(MYW)) RETURN

    IF(ALLOCATED(MYW%INDEX)) DEALLOCATE(MYW%INDEX)

    IF(ALLOCATED(MYW%PTW)) DEALLOCATE(MYW%PTW)
    
    IF(ASSOCIATED(MYW%N_INDEX)) DEALLOCATE(MYW%N_INDEX)
    IF(ASSOCIATED(MYW%N_CNT)) DEALLOCATE(MYW%N_CNT)
    IF(ASSOCIATED(MYW%N_ORDER)) DEALLOCATE(MYW%N_ORDER)
    IF(ASSOCIATED(MYW%N_FOUND)) DEALLOCATE(MYW%N_FOUND)

    DEALLOCATE(MYW)

  END SUBROUTINE KILL_WEIGHTS


  SUBROUTINE PRINT_WEIGHTS(MYW,STR)
    IMPLICIT NONE
    TYPE(INTERP_WEIGHTS) :: MYW
    CHARACTER(LEN=*) :: STR
    INTEGER :: I

    WRITE(IPT,*) "===== PRINTING WEIGHTS:"//TRIM(STR)//"========"
    WRITE(IPT,*)
    WRITE(IPT,*) "Nin=",MYW%Nin
    WRITE(IPT,*) "Nout=",MYW%Nout
    IF(ALLOCATED(MYW%index)) THEN
       WRITE(IPT,*) "Size(index)=",size(MYW%index)
    ELSE
       WRITE(IPT,*) "Index not allocated!"
    END IF

    IF(ALLOCATED(MYW%PTW)) THEN
       DO I = 1,MYW%Nout
          CALL PRINT_PTW(MYW%PTW(I),I)
       END DO
    ELSE
       WRITE(IPT,*) "PTW not ALLOCATED!"       
    END IF

    WRITE(IPT,*) "NType=",MYW%N_TYPE

    IF(Associated(MYW%N_INDEX)) THEN
       WRITE(IPT,*) "Size(N_INDEX)=",size(MYW%N_INDEX)
    ELSE
       WRITE(IPT,*) "N_INDEX not allocated!"
    END IF

    IF(Associated(MYW%N_CNT)) THEN
       WRITE(IPT,*) "Size(N_CNT)=",size(MYW%N_CNT)
    ELSE
       WRITE(IPT,*) "N_CNT not allocated!"
    END IF

    IF(Associated(MYW%N_ORDER)) THEN
       WRITE(IPT,*) "Size(N_ORDER)=",size(MYW%N_ORDER)
    ELSE
       WRITE(IPT,*) "N_ORDER not allocated!"
    END IF

    IF(Associated(MYW%N_FOUND)) THEN
       WRITE(IPT,*) "Size(N_FOUND)=",size(MYW%N_FOUND)
    ELSE
       WRITE(IPT,*) "N_FOUND not allocated!"
    END IF
    WRITE(IPT,*) ""
    WRITE(IPT,*)"========== END WEIGHTS ========"


  END SUBROUTINE PRINT_WEIGHTS

  SUBROUTINE PRINT_PTW(PTW,CNT)
    IMPLICIT NONE
    TYPE(R2PTS) :: PTW
    INTEGER :: CNT
    WRITE(IPT,*) "===== PRINTING PTW#",CNT
    WRITE(IPT,*) "WEIGHTS="
    WRITE(IPT,'(4F12.6)')PTW%WEIGHTS(:,1)
    WRITE(IPT,'(4F12.6)')PTW%WEIGHTS(:,2)
    WRITE(IPT,*) "QUADNUM="
    WRITE(IPT,'(4I8)')PTW%QUAD_NUM(:,1)
    WRITE(IPT,'(4I8)')PTW%QUAD_NUM(:,2)


  END SUBROUTINE PRINT_PTW

  ! THE BILINEAR METHOD USED HERE IS DESIGNED TO FUNCTION ON
  ! CURVALINEAR GRIDDED DATA. ANY 



  ! GET THE LOCATION OF THE NEAREST POINT AND RETURN AN INDEX TO IT
  ! ALLONG WITH ITS NEIGHBORS IN THE MESH FOR THE CELL THAT CONTAINS
  ! IT. THE CONDITION INDICATES HOW TO INTERPOLATE IN THIS CELL TO
  ! HANDLE A LAND MASK OR A VALUE OUTSIDE THE DOMAIN BOUNDARY
  SUBROUTINE GET_LOC(X_bnd,Y_bnd,msze,nsze,RSQ,MASK,X,Y,BOX,CONDITION)
    IMPLICIT NONE 
    REAL(SP), POINTER, INTENT(IN) :: X_bnd(:,:)
    REAL(SP), POINTER, INTENT(IN) :: Y_bnd(:,:)
    REAL(SP), INTENT(IN) :: X,Y
    INTEGER, INTENT(IN) :: MSZE,NSZE
    REAL(SP), POINTER :: Rsq(:,:)
    INTEGER, POINTER,INTENT(IN) :: MASK(:,:)


    INTEGER, INTENT(OUT) :: BOX(4,2)
    INTEGER, INTENT(OUT) :: CONDITION

    INTEGER :: MISSING(4)

    INTEGER, DIMENSION(2) :: NN, NNT
    INTEGER :: err,I,j,CNT
    REAL(SP) :: Xn,Yn

    REAL(SP) :: PTEST

    REAL(SP),POINTER::XIN(:,:),YIN(:,:)

    XIN => X_BND(1:msze,1:nsze)
    YIN => Y_BND(1:msze,1:nsze)

    CALL WATCH_START_LAP(MIN_W)
    ! Find Nearest point

!!$!---------------------------------------------------------------
!!$# if defined(SPHERICAL)
!!$!---------------------------------------------------------------
!!$
!!$
!!$    IF(ASSOCIATED(MASK))THEN;
!!$       
!!$       DO I=1,msze
!!$          DO J=1,nsze
!!$             IF(MASK(I,J)==0) CALL ARC(Xin(i,j),Yin(i,j),X,Y,RSQ(i,j))
!!$          END DO
!!$       END DO
!!$       
!!$       NN = minloc(Rsq,MASK==0)
!!$
!!$       ! AT THE NORTH POLE, ONLY THE LONGITUDE MATTERS...
!!$       IF( YIN(NN(1),NN(2)) == 90.0_SP .or. Y == 90.0_SP) THEN
!!$          
!!$          NN = minloc(abs(Xin - X),Yin == 90.0_SP)
!!$          
!!$       ELSEIF ( YIN(NN(1),NN(2)) == -90.0_SP .or. Y == -90.0_SP) THEN
!!$          NN = minloc(abs(Xin - X),Yin == -90.0_SP)
!!$          
!!$       END IF
!!$       
!!$    ELSE
!!$       DO I=1,msze
!!$          DO J=1,nsze
!!$             CALL ARC(Xin(i,j),Yin(i,j),X,Y,RSQ(i,j))
!!$          END DO
!!$       END DO
!!$       
!!$       NN = minloc(Rsq)
!!$ 
!!$       ! AT THE NORTH POLE, ONLY THE LONGITUDE MATTERS...
!!$       IF( YIN(NN(1),NN(2)) == 90.0_SP .or. Y == 90.0_SP) THEN
!!$          
!!$          NN = minloc(abs(Xin - X),Yin == 90.0_SP)
!!$          
!!$       ELSEIF ( YIN(NN(1),NN(2)) == -90.0_SP .or. Y == -90.0_SP) THEN
!!$          NN = minloc(abs(Xin - X),Yin == -90.0_SP)
!!$          
!!$       END IF
!!$ 
!!$
!!$   END IF
!!$# else
    IF(ASSOCIATED(MASK))THEN;
       
       WHERE(MASK==0)
          Rsq = (Xin -X)**2 + (Yin -Y)**2
       END WHERE
       
       NN = minloc(Rsq,RSQ>-1.0_SP)
       
    ELSE
       
       Rsq = (Xin -X)**2 + (Yin -Y)**2
       
       NN = minloc(Rsq,RSQ>-1.0_SP)
       
    END IF
!!$!---------------------------------------------------------------
!!$# endif
!!$!---------------------------------------------------------------
    
 
    CALL WATCH_STOP_LAP(MIN_W)

    !    PRINT*, "NN=",NN


    CALL WATCH_START_LAP(BOX_W)
    BOX = Get_BOX(X,Y,NN,X_bnd,Y_bnd,err)


    ! IN SOME STRANGE GEOMETRY, THE CLOSEST POINT IS NOT IN THE SAME CELL!
    NNT = NN
    CNT = 0
    DO WHILE(ERR/=0)
       PTEST = Rsq(NNT(1),NNT(2))
       NNT = minloc(Rsq,RSQ>PTEST)

       IF(NNT(1) == 0) EXIT
!       IF(NNT(1) == 0) THEN
!          WRITE(IPT,*) "SEARCHED WHOLE GRID - NO VALID CELLS CONTAIN THE POINT"
!          exit
!       END IF
       BOX = Get_BOX(X,Y,NNT,X_bnd,Y_bnd,err)
       
       IF(CNT > 12) EXIT
       CNT = CNT +1

    END DO
    
    ! SET THE FINALL VALUE FOR NN OR HANDLE ERRROR
    IF( ERR == 0) THEN
       NN = NNT
    ELSE
       
       IF(ASSOCIATED(MASK)) THEN
          ! IF THE MASK BLANKS OUT ALL THE NODES IN THE CELL
          ! CONTAINING THE POINT, JSUT USE THE VALUE AT THE NEAREST POINT


          ! THIS IS NOT A GOOD SOLUTION - NEED TO HAVE A BACK UP
          ! USING NEAREST NEIGHBOR WITH DATA...
          CONDITION = -7
          BOX(1,:) = NN
          return
       ELSE
          ! GET THE LOCATION
          Xn = X_bnd(nn(1),nn(2))
          Yn = Y_bnd(nn(1),nn(2))

          ! THIS SHOULD NOT HAPPEN IF THRE IS NO MASK
          write(ipt,*) "Error searching for Bilinear Interpolation to Point:"
          write(ipt,*) "X=",X,"; Y=",Y
          write(ipt,*) "Found Nearest Grid index:",NN
          write(ipt,*) "Bounds(",MSZE,NSZE,")"
          write(ipt,*) "Grid location:",Xn,Yn
          
!          WRITE(IPT,*) RSQ(:,1)
!          WRITE(IPT,*) RSQ(:,2)

          CALL FATAL_ERROR("Could not find the triangle of the four nearest nodes",&
               & "containing the point x,y ?")

       END IF
    END IF


    CALL WATCH_STOP_LAP(BOX_W)


    ! ASSESS CONDITION - IS THE POINT OUTSIDE THE DOMAIN?

    ! ASSUME CONDITION 0
    CONDITION = 0

    ! ANY VALUES IN THE BOUNDARY SET TO ZERO
    BOX(:,1) = MOD(BOX(:,1),MSZE+1)

    BOX(:,2) = MOD(BOX(:,2),NSZE+1)


    IF(ASSOCIATED(MASK)) THEN
       MISSING = 0
       DO I = 1,4

          IF(BOX(I,1)*BOX(I,2)==0) THEN
             BOX(I,:)=0
             MISSING(I) = 1
             CYCLE
          END IF


!         WRITE(IPT,*) "I",I,"; MASK(BOX(I)):",MASK(BOX(I,1),BOX(I,2))

          IF(MASK(BOX(I,1),BOX(I,2))==1) THEN
             BOX(I,:)=0
             MISSING(I) = 1
          END IF
!          WRITE(IPT,*) "I",I,"; missing;",missing

       END DO

       select case(SUM(MISSING))
       CASE(1)
          IF(ALL(MISSING ==(/0,0,0,1/)) ) THEN
             CONDITION=-1
          ELSEIF(ALL(MISSING ==(/0,0,1,0/)) ) THEN
             CONDITION=-2
          ELSEIF(ALL(MISSING ==(/0,1,0,0/)) ) THEN
             CONDITION=-3
          ELSEIF(ALL(MISSING ==(/1,0,0,0/)) ) THEN
             CONDITION=-4
          END IF
          RETURN
       CASE(2)

          IF(ALL(MISSING == (/1,0,1,0/)) ) THEN
             CONDITION=-5
          ELSEIF(ALL(MISSING == (/0,1,0,1/)) ) THEN
             CONDITION=-6
          ELSEIF(ALL(MISSING ==(/1,1,0,0/)) ) THEN
             CONDITION=6
          ELSEIF(ALL(MISSING ==(/0,1,1,0/)) ) THEN
             CONDITION=4
          ELSEIF(ALL(MISSING ==(/0,0,1,1/)) ) THEN
             CONDITION=2
          ELSEIF(ALL(MISSING ==(/1,0,0,1/)) ) THEN
             CONDITION=8
          END IF
          RETURN
       CASE(3)

          IF(ALL(MISSING ==(/1,1,1,0/)) ) THEN
             CONDITION=5
          ELSEIF(ALL(MISSING ==(/1,1,0,1/)) ) THEN
             CONDITION=7
          ELSEIF(ALL(MISSING ==(/1,0,1,1/)) ) THEN
             CONDITION=1
          ELSEIF(ALL(MISSING ==(/0,1,1,1/)) ) THEN
             CONDITION=3
          END IF
          RETURN
       CASE(4)
          CALL FATAL_ERROR&
               ("THE LAND MASK HAS LEFT NO VALID POINTS - THIS SHOULD BE IMPOSSIBLE!")
       END select

    END IF

    CALL WATCH_START_LAP(COND_W)

    ! IF ANY LESS THAN ONE, VALUES ARE OUTSIDE!
    IF(ANY(BOX==0)) THEN
       ! OUTSIDE BOUNDARY OR ON THE LAND MASK!!!


       ! TEST FOR 8,6,4,2 while one value will correctly determine value
       IF(COUNT(BOX==0)==2) THEN

          ! TEST FOR ZERO: TOP OR BOTTOM, LEFT OR RIGHT
          IF(BOX(1,1) == 0) CONDITION = 8
          IF(BOX(1,2) == 0) CONDITION = 6
          IF(BOX(3,1) == 0) CONDITION = 4
          IF(BOX(3,2) == 0) CONDITION = 2
       END IF


       ! SET BOTH INDICIES TO ZERO IF EITHER IS
       DO I = 1,4
          IF(PRODUCT(BOX(I,:))==0) BOX(I,:)=0
       END DO


       IF(COUNT(BOX==0)==6)THEN
          ! TEST THE NON ZERO CORNERS
          IF(BOX(1,1) /= 0) CONDITION = 3
          IF(BOX(2,1) /= 0) CONDITION = 1
          IF(BOX(3,1) /= 0) CONDITION = 7
          IF(BOX(4,1) /= 0) CONDITION = 5


       END IF

    END IF

    CALL WATCH_STOP_LAP(COND_W)

  END SUBROUTINE GET_LOC

  ! GET THE CELL THAT CONTAINS THE POINT, IT COULD BE ONE OF FOUR
  ! THAT SHARE A COMMON NODE....
  FUNCTION GET_BOX(X,Y,NN,X_bnd,Y_bnd,err) RESULT(BOX)
    IMPLICIT NONE
    REAL(SP), POINTER, INTENT(IN) :: X_BND(:,:)
    REAL(SP), POINTER, INTENT(IN) :: Y_BND(:,:)
    REAL(SP), INTENT(IN) :: X,Y
    INTEGER, INTENT(IN):: NN(2)
    INTEGER, INTENT(OUT) :: ERR

    INTEGER:: BOX(4,2)

    REAL(SP) :: Xn,Yn,Xt(3),Yt(3)
    INTEGER :: I,J,TRI


    INTEGER:: TN(4,2)
    INTEGER:: TNO(4,2) ! TRIANGLE OPOSITE
    INTEGER:: BN(4,2)
    INTEGER:: DN(4,2)


    ! INDEX OF POINTS TO TRINAGLE
    Tn(1,:) = (/ 1 , 0 /)
    Tn(2,:) = (/ 0 , 1 /)
    Tn(3,:) = (/-1 , 0 /)
    Tn(4,:) = (/ 0 ,-1 /)

    ! INDEX OF POINTS TO OPOSITE CORNER
    Tno(1,:) = (/ 1 , 1 /)
    Tno(2,:) = (/-1 , 1 /)
    Tno(3,:) = (/-1 ,-1 /)
    Tno(4,:) = (/ 1 ,-1 /)

    ! INDEX OF POINTS TO BOX
    Bn(1,:) = (/ 0 , 1 /)
    Bn(2,:) = (/ 1 , 1 /)
    Bn(3,:) = (/ 1 , 0 /)
    Bn(4,:) = (/ 0 , 0 /)

    ! DELTA BOX DEPENDING ON WHICH TRIANGLE
    Dn(1,:) = (/ 0 , 0 /)
    Dn(2,:) = (/-1 , 0 /)
    Dn(3,:) = (/-1 ,-1 /)
    Dn(4,:) = (/ 0 ,-1 /)


    Xt(1) = X_bnd(nn(1),nn(2))

    Yt(1) = Y_bnd(nn(1),nn(2))

    err = -1
    BOX=-1


    ! TEST EACH OF THE FOUR TRIANGLE IN CELLS SURROUNDING THE NEAREST NODE
    DO I=1,4

       J = mod(I,4)+1

       Xt(2) = X_bnd(NN(1)+TN(I,1),NN(2)+TN(I,2))
       Yt(2) = Y_bnd(NN(1)+TN(I,1),NN(2)+TN(I,2))

       Xt(3) = X_bnd(NN(1)+TN(J,1),NN(2)+TN(J,2))
       Yt(3) = Y_bnd(NN(1)+TN(J,1),NN(2)+TN(J,2))

!       WRITE(IPT,*) "======FOUND TRI==========="
!       write(ipt,*) "XT",XT
!       write(ipt,*) "YT",YT
!       write(ipt,*) "X,Y",X,Y


       IF(ISINTRI(X,Y,XT,YT)) THEN
          !          print*, "TRI#",I


          ERR = 0
          BN(:,1) = BN(:,1) +  Dn(I,1)
          BN(:,2) = BN(:,2) +  Dn(I,2)

          BOX(:,1) = BN(:,1) + NN(1) 
          BOX(:,2) = BN(:,2) + NN(2) 


          RETURN
       END IF

    END DO

    ! TEST EACH OF THE FOUR TRIANGLES IN THE OPOSITE CORNER OF CELL
    ! WITH THE NEAREST NODE
    DO I=1,4

       J = mod(I,4)+1

       !OPOSITE CORNER
       Xt(1) = X_bnd(nn(1)+TNO(I,1),nn(2)+TNO(I,2))
       Yt(1) = Y_bnd(nn(1)+TNO(I,1),nn(2)+TNO(I,2))


       ! TWO OTHERS
       Xt(2) = X_bnd(NN(1)+TN(I,1),NN(2)+TN(I,2))
       Yt(2) = Y_bnd(NN(1)+TN(I,1),NN(2)+TN(I,2))

       Xt(3) = X_bnd(NN(1)+TN(J,1),NN(2)+TN(J,2))
       Yt(3) = Y_bnd(NN(1)+TN(J,1),NN(2)+TN(J,2))

!          WRITE(IPT,*) "======FOUND TRI==========="
!          write(ipt,*) "XT",XT
!          write(ipt,*) "YT",YT
!          write(ipt,*) "X,Y",X,Y


       IF(ISINTRI(X,Y,XT,YT)) THEN
          !          print*, "TRI OPOSITE! #",I


          ERR = 0
          BN(:,1) = BN(:,1) +  Dn(I,1)
          BN(:,2) = BN(:,2) +  Dn(I,2)

          BOX(:,1) = BN(:,1) + NN(1) 
          BOX(:,2) = BN(:,2) + NN(2) 


          RETURN
       END IF

    END DO




  END FUNCTION GET_BOX

  SUBROUTINE SETUP_INTERP_BILINEAR_A(Xin,Yin,Xout,Yout,WEIGHTS, land_mask)
    IMPLICIT NONE 
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Xin(:,:)
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Yin(:,:)
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Xout(:)
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Yout(:)

    TYPE(INTERP_WEIGHTS), INTENT(OUT) :: WEIGHTS

    INTEGER, OPTIONAL, ALLOCATABLE, TARGET, INTENT(IN) :: LAND_MASK(:,:)

    REAL(SP), POINTER :: XinP(:,:)
    REAL(SP), POINTER :: YinP(:,:)
    REAL(SP), POINTER :: XoutP(:)
    REAL(SP), POINTER :: YoutP(:)

    INTEGER, POINTER :: LAND_MASKP(:,:)

    NULLIFY(XinP,YinP,XoutP,YoutP,Land_MaskP)


    IF(ALLOCATED(Xin)) XinP => Xin

    IF(ALLOCATED(Yin)) YinP => Yin

    IF(ALLOCATED(Xout)) XoutP => Xout

    IF(ALLOCATED(Yout)) YoutP => Yout

    IF(PRESENT(LAND_MASK)) THEN

       IF(ALLOCATED(LAND_MASK)) THEN
          LAND_MASKP => land_mask
       ELSE
          CALL FATAL_ERROR("PASSED LAND_MASK BUT IT IS NOT ALLOCATED") 
       END IF

       CALL SETUP_INTERP_BILINEAR_P(XinP,YinP,XoutP,YoutP,WEIGHTS,LAND_MASKP)
    ELSE
       CALL SETUP_INTERP_BILINEAR_P(XinP,YinP,XoutP,YoutP,WEIGHTS)
    END IF
  END SUBROUTINE SETUP_INTERP_BILINEAR_A


  SUBROUTINE SETUP_INTERP_BILINEAR_P(Xin,Yin,Xout,Yout,WEIGHTS,land_mask)
    USE ALL_VARS, only :NBSN,M,N,MT,NT,PAR

    IMPLICIT NONE 
    REAL(SP), POINTER, INTENT(IN) :: Xin(:,:)
    REAL(SP), POINTER, INTENT(IN) :: Yin(:,:)
    REAL(SP), POINTER, INTENT(IN) :: Xout(:)
    REAL(SP), POINTER, INTENT(IN) :: Yout(:)

    TYPE(INTERP_WEIGHTS), INTENT(OUT) :: WEIGHTS
    INTEGER, OPTIONAL, POINTER, INTENT(IN) :: LAND_MASK(:,:)

    INTEGER, POINTER :: MASK(:,:)

    INTEGER, POINTER :: N_INDEX(:,:)
    INTEGER, POINTER :: N_FOUND(:)
    INTEGER, POINTER :: N_CNT(:)
    INTEGER, POINTER :: N_ORDER(:)
    INTEGER :: mk_cnt


    INTEGER :: N_TYPE

    INTEGER :: I, J, status, lb, ub, msze,nsze

    REAL(SP) :: Xbox(4)
    REAL(SP) :: Ybox(4)
    REAL(SP) :: wghts(4)


    REAL(SP) :: DRA,DRB,DRT

    ! DELTA X/Y for each edge
    REAL(SP) :: DXT, DXB
    REAL(SP) :: DYT, DYB

    REAL(SP) :: ARCX1, ARCX2,ARCX3

    REAL(SP),POINTER::X_bnd(:,:),Y_bnd(:,:),RSQ(:,:)

    LOGICAL :: BOUNDS_FLAG = .false.

    REAL(SP) :: X,Y
    INTEGER :: BOX(4,2)
    INTEGER :: CONDITION

    ! CHECK INPUT ALLOCATION STATUS
    if(.not.associated(Xin)) CALL FATAL_ERROR("SETUP_INTERP: INPUT ARGU&
         &MENTS MUST BE ALLOCATED!")
    if(.not.associated(Yin)) CALL FATAL_ERROR("SETUP_INTERP: INPUT ARGU&
         &MENTS MUST BE ALLOCATED!")
    if(.not.associated(Xout)) CALL FATAL_ERROR("SETUP_INTERP: INPUT ARGU&
         &MENTS MUST BE ALLOCATED!")
    if(.not.associated(Yout)) CALL FATAL_ERROR("SETUP_INTERP: INPUT ARGU&
         &MENTS MUST BE ALLOCATED!")

    BOUNDS_FLAG = .FALSE.


    nullify(N_INDEX)
    nullify(N_CNT)
    nullify(N_ORDER)
    nullify(N_FOUND)

    ! GET THE DIMENSIONS
    WEIGHTS%Nout= ubound(Xout,1)
    msze=size(Xin,1)
    nsze=size(Xin,2)

    WEIGHTS%Nin  = msze*nsze


    ! ALLOCATE THE SPACE FOR THE WEIGHTS OUTPUT AND INDEX
    ALLOCATE(WEIGHTS%PTW(WEIGHTS%Nout), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR("SETUP_INTERP: COULD NOT ALLOCATE SPACE")


    ALLOCATE(RSQ(MSZE,NSZE)); RSQ = 0.0_SP
    ALLOCATE(X_bnd(0:MSZE+1,0:NSZE+1)); X_bnd=0.0_sp
    ALLOCATE(Y_bnd(0:MSZE+1,0:NSZE+1)); Y_bnd=0.0_sp


    ! WE NEED A LOCATION MATRIX WHICH CAN HANDLE VALUES OUTSIDE THE
    ! DOMAIN SO USE THE LIMITS AT EACH EDGE...

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
    IF(any(Xin < 0.0_sp) .or. any(Xin > 360.0)) CALL FATAL_ERROR&
         &("SETUP_INTERP_BILINEAR: Not expecting longitude less than 0.0",& 
         & "or greater than 360.0 in gridded data")

    IF(any(Yin < -90.0_sp) .or. any(Yin > 90.0)) CALL FATAL_ERROR&
         &("SETUP_INTERP_BILINEAR: Not expecting latitude greater than 90.0",& 
         & "or less than -90.0 in gridded data")
    
!    IF(any(Xin < 0.0_sp) .or. any(Xin > 360.0)) CALL FATAL_ERROR&
!         &("SETUP_INTERP_BILINEAR: Not expecting longitude less than 0.0",& 
!         & "or greater than 360.0")

!!$    IF(any(Yout < -90.0_sp) .or. any(Yout > 90.0)) CALL FATAL_ERROR&
!!$         &("SETUP_INTERP_BILINEAR: Not expecting latitude greater than 90.0",& 
!!$         & "or less than -90.0 in interpolated data")

    
    ! Set center values
    X_bnd(1:MSZE,1:NSZE)=Xin

    ! Set edge values

    WHERE(X_bnd(1,:) - X_bnd(MSZE,:) < 0.0_SP)
       ! WEST TO EAST
       X_bnd(0,:) = X_bnd(1,:) - 90.0_SP ! LEFT
       X_bnd(MSZE+1,:) = X_bnd(MSZE,:) + 90.0_SP ! RIGHT
    ELSEWHERE
       ! EAST TO WEST
       X_bnd(0,:) = X_bnd(1,:) + 90.0_SP ! LEFT
       X_bnd(MSZE+1,:) = X_bnd(MSZE,:) - 90.0_SP ! RIGHT
    END WHERE

    X_bnd(:,0) = X_bnd(:,1) !TOP
    X_bnd(:,NSZE+1) = X_bnd(:,NSZE) !BOTTOM

    ! Set center values
    Y_bnd(1:MSZE,1:NSZE)=Yin

    WHERE(Y_bnd(:,1) - Y_bnd(:,NSZE) < 0.0_SP)
       ! NORTH TO SOUTH
       Y_bnd(:,0) = -95.0_SP !TOP
       Y_bnd(:,NSZE+1) = 95.0_SP !BOTTOM
    ELSEWHERE
       ! SOUTH TO NORTH
       Y_bnd(:,0) = 95.0_SP !TOP
       Y_bnd(:,NSZE+1) = -95.0_SP !BOTTOM
    END WHERE

    ! Set edge values
    Y_bnd(0,:) = Y_bnd(1,:) ! LEFT
    Y_bnd(MSZE+1,:) = Y_bnd(MSZE,:) ! RIGHT

            
# else
    ! Set center values
    X_bnd(1:MSZE,1:NSZE)=Xin

    ! Set edge values
    X_bnd(0,:) = 2*X_bnd(1,:) - X_bnd(MSZE,:) ! LEFT
    X_bnd(MSZE+1,:) = 2*X_bnd(MSZE,:) - X_bnd(1,:) ! RIGHT

    X_bnd(:,0) = X_bnd(:,1) !TOP
    X_bnd(:,NSZE+1) = X_bnd(:,NSZE) !BOTTOM

    ! Set center values
    Y_bnd(1:MSZE,1:NSZE)=Yin

    ! Set edge values
    Y_bnd(0,:) = Y_bnd(1,:) ! LEFT
    Y_bnd(MSZE+1,:) = Y_bnd(MSZE,:) ! RIGHT

    Y_bnd(:,0) = 2*Y_bnd(:,1) - Y_bnd(:,NSZE) !TOP
    Y_bnd(:,NSZE+1) = 2*Y_bnd(:,NSZE) - Y_bnd(:,1) !BOTTOM



    ! CORNERS
    Y_bnd(0,0) = Y_bnd(1,0)
    Y_bnd(MSZE+1,0) = Y_bnd(MSZE,0)

    Y_bnd(0,NSZE+1) = Y_bnd(1,NSZE+1)
    Y_bnd(MSZE+1,NSZE+1) = Y_bnd(MSZE,NSZE+1)

    X_bnd(0,0) = X_bnd(0,1)
    X_bnd(MSZE+1,0) = X_bnd(MSZE+1,1)

    X_bnd(0,NSZE+1) = X_bnd(0,NSZE)
    X_bnd(MSZE+1,NSZE+1) = X_bnd(MSZE+1,NSZE)

!---------------------------------------------------------------
# endif
!---------------------------------------------------------------
# if defined(SPHERICAL)
    IF(X_bnd(1,1) - X_bnd(2,1) > 0.0_SP) THEN
       DO I = 0,MSZE
          IF(X_bnd(I,1) - X_bnd(I+1,1) < 0.0_SP) CALL FATAL_ERROR&
               &("Bilinear Interp: X_BND INPUT IS NOT MONOTONIC!")
          
          IF(X_bnd(I,NSZE) - X_bnd(I+1,NSZE) < 0.0_SP) CALL FATAL_ERROR&
               &("Bilinear Interp: X_BND INPUT IS NOT MONOTONIC!")
       END DO
    ELSE
       DO I = 0,MSZE
          IF(X_bnd(I,1) - X_bnd(I+1,1) > 0.0_SP) CALL FATAL_ERROR&
               &("Bilinear Interp: X_BND INPUT IS NOT MONOTONIC!")
          
          IF(X_bnd(I,NSZE) - X_bnd(I+1,NSZE) > 0.0_SP) CALL FATAL_ERROR&
               &("Bilinear Interp: X_BND INPUT IS NOT MONOTONIC!")
       END DO
    END IF

    IF(ANY(ABS(X_bnd(0:MSZE,1)-X_bnd(1:MSZE+1,1))>180.0_SP)) CALL FATAL_ERROR&
         &("Bilinear Interp: X Grid data crosses 0 Longitude!")

    IF(ANY(ABS(X_bnd(0:MSZE,NSZE)-X_bnd(1:MSZE+1,NSZE))>180.0_SP)) CALL FATAL_ERROR&
         &("Bilinear Interp: X Grid data crosses 0 Longitude!")


    IF(Y_bnd(1,1) - Y_bnd(1,2) > 0.0_SP) THEN
       DO I = 0,NSZE
          IF(Y_bnd(1,I) - Y_bnd(1,I+1) < 0.0_SP) CALL FATAL_ERROR&
               &("Bilinear Interp: Y_BND INPUT IS NOT MONOTONIC!")
          
          IF(Y_bnd(MSZE,I) - Y_bnd(MSZE,I+1) < 0.0_SP) CALL FATAL_ERROR&
               &("Bilinear Interp: Y_BND INPUT IS NOT MONOTONIC!")
       END DO
    ELSE
       DO I = 0,NSZE
          IF(Y_bnd(1,I) - Y_bnd(1,I+1) > 0.0_SP) CALL FATAL_ERROR&
               &("Bilinear Interp: Y_BND INPUT IS NOT MONOTONIC!")
          
          IF(Y_bnd(MSZE,I) - Y_bnd(MSZE,I+1) > 0.0_SP) CALL FATAL_ERROR&
               &("Bilinear Interp: Y_BND INPUT IS NOT MONOTONIC!")
       END DO
    END IF
# endif


! TO PRINT THE X_BND AND Y_BND ARRAY
! Uncommented by Jadon Ge
!    DO I=0,UBOUND(X_bnd,2)
!       write(ipt,*) "X",X_bnd(:,I)
!    END DO

!    DO I=0,UBOUND(X_bnd,2)
!       write(ipt,*) "Y",Y_bnd(:,I)
!    END DO

    !===================================================
    ! TIMING STUFF TO TEST OPTIMIZATION
    !===================================================

    CALL WATCH_INIT(MIN_W)
    CALL WATCH_INIT(BOX_W)
    CALL WATCH_INIT(COND_W)
    CALL WATCH_INIT(CASE_W)

    CALL WATCH_INIT(TOT_W)

    !===================================================

    N_TYPE = 0
    NULLIFY(MASK)
    IF(PRESENT(LAND_MASK))THEN
       MASK => LAND_MASK

       IF(ALL(MASK==1)) CALL FATAL_ERROR&
            &("BILINEAR INTERP CAN NOT WORK WITH EVERY POINT MASKED")

       IF(ANY(MASK<0 .OR. MASK>1)) CALL FATAL_ERROR&
            &("BILINEAR INTERP MASK VALUES ARE INCORRECT:",&
            & "SET 1 FOR MASKED (LAND)",&
            & "SET 0 FOR UNMASKED (VALID/OCEAN)")

       IF(ALL(MASK==0)) NULLIFY(MASK)


       IF(WEIGHTS%Nout == MT) THEN

          IF (.not. allocated(NBSN)) CALL FATAL_ERROR&
               &("CAN NOD USE MASK FOR A GRID VARIABLE WITHOUT RUNNING TGE FIRST")

          N_TYPE = TYPE_NODE
          ! allocate :: M, max number of nodes surrounding a node
          ALLOCATE(N_INDEX(ubound(nbsn,1),ubound(nbsn,2))); N_INDEX=0
          ALLOCATE(N_FOUND(0:MT)); N_FOUND=0
          ALLOCATE(N_CNT(M)); N_CNT=0

          IF(PAR) call warning("Possible interpolation using land mask",&
               &"This procedure works in parallel but may not report"&
               &,"Domain boundary errors which stop the code",&
               &"If your code stops with no message try running single processor mode")


       ELSEIF(WEIGHTS%Nout == NT) THEN

          IF (.not. allocated(NBSN)) CALL FATAL_ERROR&
               &("CAN NOD USE MASK FOR A GRID VARIABLE WITHOUT RUNNING TGE FIRST")

          N_TYPE = TYPE_ELEMENT
          ALLOCATE(N_INDEX(N,3)); N_INDEX=0
          ALLOCATE(N_FOUND(0:NT)); N_FOUND=0
          ALLOCATE(N_CNT(N)); N_CNT=0

          IF(PAR) call warning("Possible interpolation using land mask",&
               &"This procedure works in parallel but may not report"&
               &,"Domain boundary errors which stop the code",&
               &"If your code stops with no message try running single processor mode")



       END IF


    END IF




    do i=1,WEIGHTS%Nout

       X = Xout(I) ! THE LOCATION WE ARE SOLVING FOR!
       Y = Yout(I) 

       RSQ = -1.0_SP
       CALL GET_LOC(X_BND,Y_BND,MSZE,NSZE,RSQ,MASK,X,Y,BOX,CONDITION)


       ! GET_LOC RETURNS THE INDICIES TO THE GRID CELL CONTAINING X,Y
       ! IN A CLOCKWISE ORDER STARTING FROM UPPER LEFT


       !               |                  |
       !      VII      |       VI         |           V
       !               |                  |
       !---------------------------------------------------------
       !               |                  |
       !               |                  |
       !      VIII     |       box        |           IV        
       !               |                  |
       !---------------------------------------------------------
       !               |                  |
       !      I        |       II         |           III


       ! NEGATIVE CONDITION VALUES OCCUR WHEN A LAND MASK IS APPLIED
       !
       ! o - ocean
       ! x - land
       ! 
       ! BILINEAR INTERP IN TRIANGLE, WEIGHTED AVERAGE IN MISSING QUADRANT
       !           o o
       !   I   -   x o 
       ! 
       !           o o
       !  II   -   o x
       ! 
       !           o x
       ! III   -   o o
       ! 
       !           x o
       !  IV   -   o o
       !
       ! WEIGHTED AVERAGE BETEEN POINTS
       !           o x 
       !   V   -   x o
       ! 
       !           x o 
       !  VI   -   o x
       ! OTHER CASES OF TWO MISSING CORRISPOND TO CASES 2/4/6/8
       !
       ! ALL CASES OF THREE MISSING CORRISPOND TO CASES 1/3/5/7
       !
       ! VII  - IF THE MASK BLANKS OUT ALL NODES IN THE CELL
       ! CONTAINING THE POINT JUST USE THE NEAREST VALUE



!      WRITE(IPT,*) "CONDITION=",CONDITION

       WEIGHTS%PTW(i)%QUAD_NUM= box
       WEIGHTS%PTW(i)%WEIGHTS = 0.0_SP

       CALL WATCH_START_LAP(CASE_W)

       IF(CONDITION/=0) BOUNDS_FLAG = .true.

       !   Normal situation
       SELECT CASE(CONDITION)

       CASE(0) ! INSIDE BOUNDARIES

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
          ! DELTA X
!          DXT = (XIN(BOX(2,1),BOX(2,2)) -  XIN(BOX(1,1),BOX(1,2)))
          CALL ARCX(XIN(BOX(1,1),BOX(1,2)),YIN(BOX(1,1),BOX(1,2)),&
               & XIN(BOX(2,1),BOX(2,2)),YIN(BOX(2,1),BOX(2,2)),DXT )

!          DXB = (XIN(BOX(4,1),BOX(4,2)) -  XIN(BOX(3,1),BOX(3,2)))
          CALL ARCX(XIN(BOX(3,1),BOX(3,2)),YIN(BOX(3,1),BOX(3,2)),&
               & XIN(BOX(4,1),BOX(4,2)),YIN(BOX(4,1),BOX(4,2)),DXB )

          ! DELTA Y
          DYT = (YIN(BOX(2,1),BOX(2,2)) -  YIN(BOX(1,1),BOX(1,2)))*TPI

          DYB = (YIN(BOX(4,1),BOX(4,2)) -  YIN(BOX(3,1),BOX(3,2)))*TPI

# else
          ! DELTA X
          DXT = (XIN(BOX(2,1),BOX(2,2)) -  XIN(BOX(1,1),BOX(1,2)))

          DXB = (XIN(BOX(4,1),BOX(4,2)) -  XIN(BOX(3,1),BOX(3,2)))

          ! DELTA Y
          DYT = (YIN(BOX(2,1),BOX(2,2)) -  YIN(BOX(1,1),BOX(1,2)))

          DYB = (YIN(BOX(4,1),BOX(4,2)) -  YIN(BOX(3,1),BOX(3,2)))

!---------------------------------------------------------------
# endif
!---------------------------------------------------------------


          ! CHECK FOR REASONABLE SLOPES!
          IF(ABS(DYT) < TB_TOL/2._SP *ABS(DXT) .and. ABS(DYB) < TB_TOL/2._SP *ABS(DXB)) THEN
             
             
             XBOX(1) = XIN(BOX(1,1),BOX(1,2))
             XBOX(2) = XIN(BOX(2,1),BOX(2,2))
             XBOX(3) = XIN(BOX(3,1),BOX(3,2))
             XBOX(4) = XIN(BOX(4,1),BOX(4,2))
             
             YBOX(1) = YIN(BOX(1,1),BOX(1,2))
             YBOX(2) = YIN(BOX(2,1),BOX(2,2))
             YBOX(3) = YIN(BOX(3,1),BOX(3,2))
             YBOX(4) = YIN(BOX(4,1),BOX(4,2))
             
             
             
             CALL INTERP_0(Xbox,Ybox,X,Y,wghts)
             
             WEIGHTS%PTW(i)%WEIGHTS(1,1) = wghts(1)
             WEIGHTS%PTW(i)%WEIGHTS(2,1) = wghts(2)
             WEIGHTS%PTW(i)%WEIGHTS(3,1) = wghts(3)
             WEIGHTS%PTW(i)%WEIGHTS(4,1) = wghts(4)
             
          ELSE ! TRY ROTATION!

              ! CHECK FOR REASONABLE SLOPES!
             IF(ABS(DYT) > TB_TOL *ABS(DXT)) THEN
                CALL FATAL_ERROR&
                     &("THIS CURVALINEAR MESH HAS CELLS WHICH ARE TOO EXTREME FOR THIS METHOD",&
                     & "PLEASE EXAMINE TOP_DIR/testing/interp",&
                     & "IF YOU WOULD LIKE TO ADD SUPPORT FOR EXTREME CURVALINEAR INTERPOLATION")
             END IF


             ! ROTATION METHOD DOES NOT APPEAR TO WORK PROPERLY -
             ! LEAVE AS BASIS FOR FUTURE WORK
             
!!$             write(ipt,*) "BOX-pre",BOX
!!$             BOX = cshift(BOX,1)
!!$             write(ipt,*) "BOX-post",BOX
!!$
!!$
!!$             ! DELTA X
!!$             DXT = (XIN(BOX(2,1),BOX(2,2)) -  XIN(BOX(1,1),BOX(1,2)))
!!$             
!!$             DXB = (XIN(BOX(4,1),BOX(4,2)) -  XIN(BOX(3,1),BOX(3,2)))
!!$             
!!$             ! DELTA Y
!!$             DYT = (YIN(BOX(2,1),BOX(2,2)) -  YIN(BOX(1,1),BOX(1,2)))
!!$             
!!$             DYB = (YIN(BOX(4,1),BOX(4,2)) -  YIN(BOX(3,1),BOX(3,2)))
!!$             
!!$             
!!$             
!!$             
!!$             ! CHECK FOR REASONABLE SLOPES!
!!$             IF(ABS(DYT) > TB_TOL *ABS(DXT)) THEN
!!$                CALL FATAL_ERROR("CAN NOT FIND A GOOD ROTATION TO INTE&
!!$                     &RPOLATE THIS CELL!")
!!$             END IF
!!$             
!!$             IF(ABS(DYB) > TB_TOL *ABS(DXB)) THEN
!!$                CALL FATAL_ERROR("CAN NOT FIND A GOOD ROTATION TO INTE&
!!$                     &RPOLATE THIS CELL!")
!!$             END IF
!!$             
!!$
!!$             XBOX(1) = XIN(BOX(1,1),BOX(1,2))
!!$             XBOX(2) = XIN(BOX(2,1),BOX(2,2))
!!$             XBOX(3) = XIN(BOX(3,1),BOX(3,2))
!!$             XBOX(4) = XIN(BOX(4,1),BOX(4,2))
!!$             
!!$             YBOX(1) = YIN(BOX(1,1),BOX(1,2))
!!$             YBOX(2) = YIN(BOX(2,1),BOX(2,2))
!!$             YBOX(3) = YIN(BOX(3,1),BOX(3,2))
!!$             YBOX(4) = YIN(BOX(4,1),BOX(4,2))
!!$             
!!$             
!!$             
!!$             CALL INTERP_0(Xbox,Ybox,X,Y,wghts)
!!$             
!!$             wghts = cshift(wghts,-1)
!!$             
!!$
!!$             WEIGHTS%PTW(i)%WEIGHTS(1,1) = wghts(1)
!!$             WEIGHTS%PTW(i)%WEIGHTS(2,1) = wghts(2)
!!$             WEIGHTS%PTW(i)%WEIGHTS(3,1) = wghts(3)
!!$             WEIGHTS%PTW(i)%WEIGHTS(4,1) = wghts(4)
!!$             
          END IF

          


       CASE(1)

          WEIGHTS%PTW(i)%WEIGHTS(2,1) = 1.0_SP

       CASE(2)

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------

          call arcx(X,Y,XIN(BOX(2,1),BOX(2,2)),YIN(BOX(2,1),BOX(2,2)),arcx1)
          call arcx(XIN(BOX(1,1),BOX(1,2)),YIN(BOX(1,1),BOX(1,2)),X,Y,arcx2)
             
          arcx1 = abs(arcx1)
          arcx2 = abs(arcx2)
          arcx3 = arcx1+arcx2
          
          IF(ARCX3 > SMALL) THEN
!          WEIGHTS%PTW(i)%WEIGHTS(1,1)= (XIN(BOX(2,1),BOX(2,2)) - X)/(XIN(BOX(2,1),BOX(2,2)) - XIN(BOX(1,1),BOX(1,2)))
             WEIGHTS%PTW(i)%WEIGHTS(1,1)= arcx1/arcx3

!          WEIGHTS%PTW(i)%WEIGHTS(2,1)= (X - XIN(BOX(1,1),BOX(1,2)))/(XIN(BOX(2,1),BOX(2,2)) - XIN(BOX(1,1),BOX(1,2)))
             WEIGHTS%PTW(i)%WEIGHTS(2,1)= arcx2/arcx3
          ELSE
             WEIGHTS%PTW(i)%WEIGHTS(1,1)= 0.5_SP
             WEIGHTS%PTW(i)%WEIGHTS(2,1)= 0.5_SP
          END IF

# else
          WEIGHTS%PTW(i)%WEIGHTS(1,1)= (XIN(BOX(2,1),BOX(2,2)) - X)/(XIN(BOX(2,1),BOX(2,2)) - XIN(BOX(1,1),BOX(1,2)))

          WEIGHTS%PTW(i)%WEIGHTS(2,1)= (X - XIN(BOX(1,1),BOX(1,2)))/(XIN(BOX(2,1),BOX(2,2)) - XIN(BOX(1,1),BOX(1,2)))
!--------------------------------------------------------------------
# endif
!--------------------------------------------------------------------

       CASE(3)

          WEIGHTS%PTW(i)%WEIGHTS(1,1) = 1.0_SP

       CASE(4)

          ! RATIO OF THE DIFFERENCE IN LATITUDE IS THE SAME IN SPHERICAL
          WEIGHTS%PTW(i)%WEIGHTS(1,1)= (YIN(BOX(4,1),BOX(4,2)) - Y)/(YIN(BOX(4,1),BOX(4,2)) - YIN(BOX(1,1),BOX(1,2)))

          WEIGHTS%PTW(i)%WEIGHTS(4,1)= (Y - YIN(BOX(1,1),BOX(1,2)))/(YIN(BOX(4,1),BOX(4,2)) - YIN(BOX(1,1),BOX(1,2)))


       CASE(5)

          WEIGHTS%PTW(i)%WEIGHTS(4,1) = 1.0_SP

       CASE(6)
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
          call arcx(X,Y,XIN(BOX(4,1),BOX(4,2)),YIN(BOX(4,1),BOX(4,2)),arcx1)
          call arcx(XIN(BOX(3,1),BOX(3,2)),YIN(BOX(3,1),BOX(3,2)),X,Y,arcx2)

          arcx1 = abs(arcx1)
          arcx2 = abs(arcx2)
          arcx3 = arcx1+arcx2
           

          IF(ARCX3>small) THEN
!          WEIGHTS%PTW(i)%WEIGHTS(3,1)= (XIN(BOX(4,1),BOX(4,2)) - X)/(XIN(BOX(4,1),BOX(4,2)) - XIN(BOX(3,1),BOX(3,2)))
             WEIGHTS%PTW(i)%WEIGHTS(3,1)= arcx1/arcx3
          

!          WEIGHTS%PTW(i)%WEIGHTS(4,1)= (X - XIN(BOX(3,1),BOX(3,2)))/(XIN(BOX(4,1),BOX(4,2)) - XIN(BOX(3,1),BOX(3,2)))
             WEIGHTS%PTW(i)%WEIGHTS(4,1)=arcx2/arcx3
          ELSE
             WEIGHTS%PTW(i)%WEIGHTS(3,1)= 0.5_SP
             WEIGHTS%PTW(i)%WEIGHTS(4,1)= 0.5_SP
          END IF
# else 
          WEIGHTS%PTW(i)%WEIGHTS(3,1)= (XIN(BOX(4,1),BOX(4,2)) - X)/(XIN(BOX(4,1),BOX(4,2)) - XIN(BOX(3,1),BOX(3,2)))

          WEIGHTS%PTW(i)%WEIGHTS(4,1)= (X - XIN(BOX(3,1),BOX(3,2)))/(XIN(BOX(4,1),BOX(4,2)) - XIN(BOX(3,1),BOX(3,2)))
!--------------------------------------------------------------------
# endif
!--------------------------------------------------------------------
       CASE(7)

          WEIGHTS%PTW(i)%WEIGHTS(3,1) = 1.0_SP

       CASE(8)

          ! RATIO OF THE DIFFERENCE IN LATITUDE IS THE SAME IN SPHERICAL
          WEIGHTS%PTW(i)%WEIGHTS(2,1)= (Y - YIN(BOX(3,1),BOX(3,2)))/(YIN(BOX(2,1),BOX(2,2)) - YIN(BOX(3,1),BOX(3,2)))

          WEIGHTS%PTW(i)%WEIGHTS(3,1)= (YIN(BOX(2,1),BOX(2,2)) - Y)/(YIN(BOX(2,1),BOX(2,2)) - YIN(BOX(3,1),BOX(3,2)))

       CASE(-1)

          ! BOTTOM LEFT IS MISSING
          ! FOR THE TWO EDGES WITH DATA
 
          XBOX(1) = X_BND(BOX(1,1),BOX(1,2))
          XBOX(2) = X_BND(BOX(2,1),BOX(2,2))
          XBOX(3) = X_BND(BOX(3,1),BOX(3,2))
          XBOX(4) = X_BND(BOX(4,1),BOX(4,2))
          
          YBOX(1) = Y_BND(BOX(1,1),BOX(1,2))
          YBOX(2) = Y_BND(BOX(2,1),BOX(2,2))
          YBOX(3) = Y_BND(BOX(3,1),BOX(3,2))
          YBOX(4) = Y_BND(BOX(4,1),BOX(4,2))
          
          CALL INTERP_NEG(Xbox,Ybox,X,Y,wghts)
          
          WEIGHTS%PTW(i)%WEIGHTS(1,1) = wghts(1)
          WEIGHTS%PTW(i)%WEIGHTS(2,1) = wghts(2)
          WEIGHTS%PTW(i)%WEIGHTS(3,1) = wghts(3)
          WEIGHTS%PTW(i)%WEIGHTS(4,1) = wghts(4)

       CASE(-2)

          ! BOTTOM RIGHT IS MISSING
          BOX = CSHIFT(BOX,-1)

          XBOX(1) = X_BND(BOX(1,1),BOX(1,2))
          XBOX(2) = X_BND(BOX(2,1),BOX(2,2))
          XBOX(3) = X_BND(BOX(3,1),BOX(3,2))
          XBOX(4) = X_BND(BOX(4,1),BOX(4,2))
          
          YBOX(1) = Y_BND(BOX(1,1),BOX(1,2))
          YBOX(2) = Y_BND(BOX(2,1),BOX(2,2))
          YBOX(3) = Y_BND(BOX(3,1),BOX(3,2))
          YBOX(4) = Y_BND(BOX(4,1),BOX(4,2))
          
          CALL INTERP_NEG(Xbox,Ybox,X,Y,wghts)
          
          WGHTS = CSHIFT(wghts,1)

          WEIGHTS%PTW(i)%WEIGHTS(1,1) = wghts(1)
          WEIGHTS%PTW(i)%WEIGHTS(2,1) = wghts(2)
          WEIGHTS%PTW(i)%WEIGHTS(3,1) = wghts(3)
          WEIGHTS%PTW(i)%WEIGHTS(4,1) = wghts(4)
          
 
      CASE(-3)

          ! TOP RIGHT IS MISSING

          BOX = CSHIFT(BOX,-2)

          XBOX(1) = X_BND(BOX(1,1),BOX(1,2))
          XBOX(2) = X_BND(BOX(2,1),BOX(2,2))
          XBOX(3) = X_BND(BOX(3,1),BOX(3,2))
          XBOX(4) = X_BND(BOX(4,1),BOX(4,2))
          
          YBOX(1) = Y_BND(BOX(1,1),BOX(1,2))
          YBOX(2) = Y_BND(BOX(2,1),BOX(2,2))
          YBOX(3) = Y_BND(BOX(3,1),BOX(3,2))
          YBOX(4) = Y_BND(BOX(4,1),BOX(4,2))
          
          CALL INTERP_NEG(Xbox,Ybox,X,Y,wghts)
          
          WGHTS = CSHIFT(wghts,2)

          WEIGHTS%PTW(i)%WEIGHTS(1,1) = wghts(1)
          WEIGHTS%PTW(i)%WEIGHTS(2,1) = wghts(2)
          WEIGHTS%PTW(i)%WEIGHTS(3,1) = wghts(3)
          WEIGHTS%PTW(i)%WEIGHTS(4,1) = wghts(4)
          

      CASE(-4)

         ! TOP LEFT IS MISSING
 
         BOX = CSHIFT(BOX,-3)

 
          XBOX(1) = X_BND(BOX(1,1),BOX(1,2))
          XBOX(2) = X_BND(BOX(2,1),BOX(2,2))
          XBOX(3) = X_BND(BOX(3,1),BOX(3,2))
          XBOX(4) = X_BND(BOX(4,1),BOX(4,2))
          
          YBOX(1) = Y_BND(BOX(1,1),BOX(1,2))
          YBOX(2) = Y_BND(BOX(2,1),BOX(2,2))
          YBOX(3) = Y_BND(BOX(3,1),BOX(3,2))
          YBOX(4) = Y_BND(BOX(4,1),BOX(4,2))
          
          CALL INTERP_NEG(Xbox,Ybox,X,Y,wghts)
          
          WGHTS = CSHIFT(wghts,3)

          WEIGHTS%PTW(i)%WEIGHTS(1,1) = wghts(1)
          WEIGHTS%PTW(i)%WEIGHTS(2,1) = wghts(2)
          WEIGHTS%PTW(i)%WEIGHTS(3,1) = wghts(3)
          WEIGHTS%PTW(i)%WEIGHTS(4,1) = wghts(4)
          

       CASE(-5) 

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
          CALL ARC(X,Y,X_BND(BOX(2,1),BOX(2,2)),Y_BND(BOX(2,1),BOX(2,2)),DRA)
          CALL ARC(X,Y,X_BND(BOX(4,1),BOX(4,2)),Y_BND(BOX(4,1),BOX(4,2)),DRB)
# else
          DRA = SQRT( (X_BND(BOX(2,1),BOX(2,2)) - X)**2 + (Y_BND(BOX(2,1),BOX(2,2)) - Y)**2)
          DRB = SQRT( (X_BND(BOX(4,1),BOX(4,2)) - X)**2 + (Y_BND(BOX(4,1),BOX(4,2)) - Y)**2)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

          DRT = (DRA +DRB)
          ! NORMALIZE TO ONE
          WEIGHTS%PTW(i)%WEIGHTS(2,1) = DRB/DRT

          WEIGHTS%PTW(i)%WEIGHTS(4,1) = DRA/DRT

       CASE(-6) 

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
          CALL ARC(X,Y,X_BND(BOX(1,1),BOX(1,2)),Y_BND(BOX(1,1),BOX(1,2)),DRA)
          CALL ARC(X,Y,X_BND(BOX(3,1),BOX(3,2)),Y_BND(BOX(3,1),BOX(3,2)),DRB)
          
# else
          DRA = SQRT( (X_BND(BOX(1,1),BOX(1,2)) - X)**2 + (Y_BND(BOX(1,1),BOX(1,2)) - Y)**2)
          DRB = SQRT( (X_BND(BOX(3,1),BOX(3,2)) - X)**2 + (Y_BND(BOX(3,1),BOX(3,2)) - Y)**2)

!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

          DRT = (DRA +DRB)
          ! NORMALIZE TO ONE
          WEIGHTS%PTW(i)%WEIGHTS(1,1) = DRB/DRT

          WEIGHTS%PTW(i)%WEIGHTS(3,1) = DRA/DRT


       CASE(-7) 

          ! NOTE: THIS METHOD IS DESIGNED TO SEARCH FOR THE NEAREST
          ! MESH NEIGHBOR THAT HAS AN UNMASKED VALUES AND TAKE THE
          ! AVERAGE OF THESE VALUES, PROPAGATING SEQUENTIALLY INTO
          ! THE MASKED REGION. THIS CAN BE SUCCESSFUL IN SMALL
          ! REGIONS OF MASKED DATA LIKE A RIVER. IT CAN NOT MAKE UP
          ! PHYSICS OR DATA AND SHOULD GENERALLY NOT BE TRUSTED!

          IF(N_TYPE==TYPE_NODE) THEN
             N_FOUND(i) = -1
             WEIGHTS%PTW(i)%QUAD_NUM = 0
             IF(I>M) CALL FATAL_ERROR&
                  &("INTERPOLATION WITH MASK: ERORR!",&
                  & "NEARST MESH NEIGHBOR SEARCH CAN NOT WORK WHEN THE MASK COVERS A DOMAIN BOUNDARY!")
             
          ELSEIF( N_TYPE==TYPE_ELEMENT) THEN
             N_FOUND(i) = -1
             WEIGHTS%PTW(i)%QUAD_NUM = 0
             IF(I>N) CALL FATAL_ERROR&
                  &("INTERPOLATION WITH MASK: ERORR!",&
                  & "NEARST MESH NEIGHBOR SEARCH CAN NOT WORK WHEN THE MASK COVERS A DOMAIN BOUNDARY!")

          ELSE
             ! IF NOT USING FVCOM GRID USE NEAREST VALUE BY DISTANCE
             WEIGHTS%PTW(i)%WEIGHTS(1,1) = 1.0_SP
             WEIGHTS%PTW(i)%QUAD_NUM(2:4,:) = 0
             
          END IF



       CASE DEFAULT

          CALL FATAL_ERROR("SETUP_INTERP_BILINEAR: INVALID CONDITION" )
       END SELECT

       CALL WATCH_STOP_LAP(CASE_W)

    END do

    
    ! TAKE CARE OF ANY MARKED AS MASKED
    IF(N_TYPE==TYPE_NODE) THEN
       IF(ANY(N_FOUND < 0)) THEN
          
          CALL WARNING&
               &("FVCOM ATTEMPTING TO GUESS VALUES FOR A REGION WHERE ALL DATA ARE MASKED",&
               & "PLEASE CHECK RESULTS CAREFULLY AND GET BETTER DATA IF POSSIBLE")
          
          mk_cnt = abs(sum(n_found))


          ALLOCATE(N_ORDER(mk_cnt))

          ! IMPORTANT INTIALIZATION TO PREVENT COUNTING BOUNDARY AS
          ! FOUND CELL!
          N_FOUND(0) = -2

          CALL GRID_NEIGHBOR_INDEX(N_FOUND,N_INDEX,N_CNT,N_ORDER)

          WEIGHTS%N_INDEX => N_INDEX
          WEIGHTS%N_CNT => N_CNT
          WEIGHTS%N_ORDER => N_ORDER
          WEIGHTS%N_FOUND => N_FOUND

          nullify(N_INDEX)
          nullify(N_CNT)
          nullify(N_ORDER)
          nullify(N_FOUND)

       ELSE

          DEALLOCATE(N_FOUND)
          DEALLOCATE(N_INDEX)
          DEALLOCATE(N_CNT)

       END IF

    ELSEIF( N_TYPE==TYPE_ELEMENT) THEN
       IF(ANY(N_FOUND < 0)) THEN
          CALL WARNING&
            &("FVCOM ATTEMPTING TO GUESS VALUES FOR A REGION WHERE ALL DATA ARE MASKED",&
            & "PLEASE CHECK RESULTS CAREFULLY AND GET BETTER DATA IF POSSIBLE")
       
          mk_cnt = abs(sum(n_found))
          

          ALLOCATE(N_ORDER(mk_cnt))
          
          ! IMPORTANT INTIALIZATION TO PREVENT COUNTING BOUNDARY AS
          ! FOUND CELL!
          N_FOUND(0) = -2

          CALL GRID_NEIGHBOR_INDEX(N_FOUND,N_INDEX,N_CNT,N_ORDER)

          WEIGHTS%N_INDEX => N_INDEX
          WEIGHTS%N_CNT => N_CNT
          WEIGHTS%N_ORDER => N_ORDER
          WEIGHTS%N_FOUND => N_FOUND

          nullify(N_INDEX)
          nullify(N_CNT)
          nullify(N_ORDER)
          nullify(N_FOUND)

       ELSE

          DEALLOCATE(N_FOUND)
          DEALLOCATE(N_INDEX)
          DEALLOCATE(N_CNT)


       END IF

    END IF


    IF(DBG_SET(DBG_SBR)) THEN
       WRITE(IPT,*) "///////////////////////////////////////////////////"
       WRITE(IPT,*) "//// TIMING REPORT FROM BILINEAR INTERP ///////////"
       WRITE(IPT,*) "///////////////////////////////////////////////////"
       CALL WATCH_REPORT(MIN_W,IPT,"MIN LOCATION")

       CALL WATCH_REPORT(BOX_W,IPT,"GET BOX")

       CALL WATCH_REPORT(COND_W,IPT,"CONDITION")

       CALL WATCH_REPORT(CASE_W,IPT,"WEIGHTS")

       CALL WATCH_LAP(TOT_W,IPT,"TOTAL")
       WRITE(IPT,*) "///////////////////////////////////////////////////"
    END IF

    DEALLOCATE(X_BND)
    DEALLOCATE(Y_BND)
    DEALLOCATE(RSQ)


    IF(BOUNDS_FLAG) CALL WARNING&
         &("FVCOM GRID IS OUT OF BOUNDS FOR FORCING",&
         & "USING CONSTANT VALUE OUTSIDE OF FORCING REGION")

  END SUBROUTINE SETUP_INTERP_BILINEAR_P

  ! CALCULATE WEIGHTS FOR CURVILINEAR INTERP
  SUBROUTINE INTERP_0(Xbox,Ybox,X,Y,wghts)
    IMPLICIT NONE 

    REAL(SP), INTENT(IN) :: Xbox(4)
    REAL(SP), INTENT(IN) :: Ybox(4)
    REAL(SP), INTENT(IN) :: X
    REAL(SP), INTENT(IN) :: Y
    REAL(SP), INTENT(OUT) :: wghts(4)


    ! FOR GENERAL CASE
    REAL(SP) :: AB,BB,AT,BT ! SLOPE AND INTERCEPT, TOP and BOTTOM
    REAL(SP) :: AL,BL,AR,BR ! SLOPE AND INTERCEPT, TOP and BOTTOM

    ! DELTA X/Y for each edge
    REAL(SP) :: DXL, DXR, DXT, DXB
    REAL(SP) :: DYL, DYR, DYT, DYB

    ! X/Y INTERCEPT ON TOP AND BOTTOM EDGES
    REAL(sp) :: XT,YT,XB,YB

    ! SLOPE, INTERPCEPT and DELTA's FOR AN 'AVERAGE' LINE BETWEEN TOP
    ! AND BOTTOM 
    REAL(SP) :: XLA,XRA
    REAL(SP) :: DYM, DXM
    REAL(SP) :: AM,BM
    REAL(SP) :: RADL,RADR,RADM

    ! DISTANCES
    REAL(SP) :: DR1, DR2,DR3,DR4,DR5,DR6
    REAL(SP) :: DRM,DRT,DRB
    
    ! FOR SPHERICAL MODEL, DO INTERCEPT CALCULATION IN LATITUDE AND LONGITUDE!

    ! DELTA X
    DXT = XBOX(2) -  XBOX(1)
    
    DXB = XBOX(4) -  XBOX(3)
    
    DXL = XBOX(1) -  XBOX(4)
    
    DXR = XBOX(2) -  XBOX(3)
    
    ! DELTA Y
    DYT = YBOX(2) -  YBOX(1)
    
    DYB = YBOX(4) -  YBOX(3)
    
    DYL = YBOX(1) -  YBOX(4)
    
    DYR = YBOX(2) -  YBOX(3)
    

    !SLOPE AND INTERCEPT FOR TOP AND BOTTOM OF QUAD
    AT = DYT/DXT
    BT = YBOX(2) - AT*XBOX(2)

    !       write(ipt,*) "AT=",AT,"; BT=",BT


    AB = DYB/DXB
    BB = YBOX(4) - AB*XBOX(4)

    !       write(ipt,*) "AB=",AB,"; BB=",BB


    !GUESS AT X INTERCEPT ON TOP AND BOTTOM LINES
    XT = X

    XB = X


    ! IF THE SIDES ARE NOT NEARLY VERTICAL GET A NEW ESTIMATE FOR XT, XB
    IF(ABS(DYL) < LR_TOL * ABS(DXL) .OR. ABS(DYR) < LR_TOL*ABS(DXR) ) THEN

       !XLA = (XBOX(BOX(1,1),BOX(1,2)) + XBOX(BOX(4,1),BOX(4,2)) ) / 2.0_SP

       IF(ABS(DYL) > LR_TOL * ABS(DXL)) THEN ! IF THE LEFT SIDE IS HORIZONTAL
          XLA = MAX(XBOX(1),XBOX(4))
       ELSE
          ! USE INTERCEPT AT Y!
          AL = DYL/DXL
          BL = YBOX(4) - AL*XBOX(4)

          XLA = (Y -BL)/AL
       END IF


       IF(ABS(DYR) > LR_TOL * ABS(DXR)) THEN ! IF THE LEFT SIDE IS HORIZONTAL
          XRA = MIN(XBOX(2),XBOX(3))
       ELSE
          ! USE INTERCEPT AT Y!
          AR = DYR/DXR
          BR = YBOX(2) - AR*XBOX(2)

          XRA = (Y -BR)/AR
       END IF


!       WRITE(IPT,*) "XLA=",XLA,"; XRA=",XRA

       !             write(ipt,*) "Left slope=",(DYL)/(DXL)
       !             write(ipt,*) "Right slope=",(DYR)/(DXR)


       RADL = ATAN2(DYL,DXL)
       RADR = ATAN2(DYR,DXR)

       !             write(ipt,*) "Left RAD=",RADL
       !             write(ipt,*) "Right RAD=",RADR


       ! GET AN AVERAGE DX and DY
       RADM  = (XLA -XT)*(RADR)/(XLA-XRA) + (XT - XRA)*RADL/(XLA-XRA)

       !             write(ipt,*) "RADM=",RADM

       AM = TAN(RADM)

       ! IF THE RESULTING SLOPE IS NOT VERTICAL
       IF(AM < LR_TOL )THEN
          BM = Y - AM *X

          XT = (BT-BM)/(AM - AT)

          XB = (BB-BM)/(AM - AB)

          !                write(ipt,*) "Middle! AM=",AM


       END IF

    END IF

    !GET YT, YB FROM XB,XT
    YB = AB*XB+BB 

    YT = AT*XT+BT 


!    WRITE(IPT,*) "XT=",XT,"; YT=",YT

!    WRITE(IPT,*) "XB=",XB,"; YB=",YB



!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------

    !  NOW CONVERT TO DISTANCES IN SPHERICAL COORDINATES
    ! DELTA X

    ! LENGTH OF TOP AND BOTTOM
    CALL ARC(XBOX(1),YBOX(1),XBOX(2),YBOX(2),DRT)
    CALL ARC(XBOX(3),YBOX(3),XBOX(4),YBOX(4),DRB)
    ! MIDDLE LINE LENGTH
    CALL ARC(XB,YB,XT,YT,DRM)
    
    ! TOP SEGMENTS
    IF(DRT > SMALL) THEN
       CALL ARC(XBOX(1),YBOX(1),XT,YT,DR1)
       DR2 = DRT - DR1
    ELSE
       DRT = 1.0_SP
       DR1 = 0.5_SP
       DR2 = 0.5_SP
    END IF
    
    ! BOTTOM SEGMENTS
    IF(DRB> SMALL) THEN
       CALL ARC(XBOX(3),YBOX(3),XB,YB,DR3)
       DR4 = DRB - DR3
    ELSE
       DRB = 1.0_SP
       DR3 = 0.5_SP
       DR4 = 0.5_SP
    END IF


    ! MIDDLE SEGMENTS
    IF(DRM > SMALL) THEN
       CALL ARC(X,Y,XT,YT,DR5)
       DR6 = DRM-DR5
    ELSE
       DRM = 1.0_SP
       DR5 = 0.5_SP
       DR6 = 0.5_SP
    END IF

# else

    ! TOP AND BOTTOM LENGTH
    DRT = SQRT(DXT**2+ DYT**2)
    DRB = SQRT(DXB**2+ DYB**2)

    ! TOP SEGMENTS
    DR1 = SQRT( (XBOX(1) - XT)**2 + (YBOX(1) - YT)**2)
    DR2 = DRT - DR1

    ! BOTTOM SEGMENTS
    DR3 = SQRT( (XBOX(3) - XB)**2 + (YBOX(3) - YB)**2)
    DR4 = DRB - DR3

    ! MIDDLE LINE LENGTH
    DRM = SQRT((YB-YT)**2 + (XB-XT)**2)

    ! MIDDLE SEGMENTS
    DR5 = SQRT((Y-YT)**2 + (X-XT)**2)
    DR6 = DRM-DR5

!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

!    WRITE(IPT,*) "DRT",DRT
!    WRITE(IPT,*) "DRB",DRB
!    WRITE(IPT,*) "DRM",DRM
!    WRITE(IPT,*) "DR1",DR1
!    WRITE(IPT,*) "DR2",DR2
!    WRITE(IPT,*) "DR3",DR3
!    WRITE(IPT,*) "DR4",DR4
!    WRITE(IPT,*) "DR5",DR5
!    WRITE(IPT,*) "DR6",DR6




    wghts(1) =(DR6*DR2) / (DRT*DRM)
    wghts(2) =(DR6*DR1) / (DRT*DRM)
    wghts(3) =(DR5*DR4) / (DRB*DRM)
    wghts(4) =(DR5*DR3) / (DRB*DRM)

  END SUBROUTINE INTERP_0

  SUBROUTINE INTERP_NEG(Xbox,Ybox,X,Y,wghts)
    IMPLICIT NONE 

    REAL(SP), INTENT(IN) :: Xbox(4)
    REAL(SP), INTENT(IN) :: Ybox(4)
    REAL(SP), INTENT(IN) :: X
    REAL(SP), INTENT(IN) :: Y
    REAL(SP), INTENT(OUT) :: wghts(4)


    ! FOR GENERAL CASE
    REAL(SP) :: A1,B1,A2,B2,A3,B3 ! SLOPE AND INTERCEPT, TOP and BOTTOM

    ! DELTA X/Y for each edge
    REAL(SP) :: DX1, DX2, DX3
    REAL(SP) :: DY1, DY2, DY3

    ! X/Y INTERCEPT ON EDGES
    REAL(sp) :: X1,Y1,X2,Y2

    ! SLOPE, INTERPCEPT and DELTA's FOR AN 'AVERAGE' LINE BETWEEN TOP
    ! AND BOTTOM 

    ! DISTANCES
    REAL(SP) :: DRA, DRB,DRC,DRD,DRE,DRF
    REAL(SP) :: DR1,DR2,DR3


    REAL(SP) :: Xt(3),Yt(3)


    XT = XBOX(1:3)
    YT = YBOX(1:3)

    IF(ISINTRI(X,Y,XT,YT))THEN
       ! BOTTOM LEFT IS MISSING
       ! FOR THE TWO EDGES WITH DATA
       ! DELTA X
       DX1 = XBOX(2) -  XBOX(1) !SIDE ONE
       DX2 = XBOX(2) -  XBOX(3) !SIDE TWO

       ! DELTA Y 
       DY1 = YBOX(2) -  YBOX(1)
       DY2 = YBOX(2) -  YBOX(3)


       ! ACCROSS THE CENTER - THE TWO VALID POINTS!
       DX3 = XBOX(3) -  XBOX(1) 
       DY3 = YBOX(3) -  YBOX(1)


       ! USING THE SAME SLOPE AS THE MIDSECTION, GET THE INTERCEPT OF
       ! THE POINT X/Y ON THE VALID EDGES

       ! IF THE MID SECTION IS NEARLY VERTICAL - GREAT
       IF(ABS(DY3) > TB_TOL *ABS(DX3)) THEN
          X1 = X
          X2 = X

          ! assume the other sides are not also vertical!
          A1 = DX1/DX1
          B1 = YBOX(2) - A1*XBOX(2)

          ! GET YT 
          Y1 = A1*X1+B1 


          A2 = DY2/DX2
          B2 = YBOX(2) - A2*XBOX(2)
          ! GET YB
          Y2 = A2*X2+B2 

       ELSE
          ! GET THE SLOPE OF THE MIDDLE LINE, AND THE EQ FOR THE LINE
          ! THROUGH THE POINT X/Y WITH THAT SLOPE
          A3 = DY3 / DX3
          B3 = Y - A3 *X

!          WRITE(IPT,*) "A3=",A3,"; BM=",B3

          ! SIDE ONE
          IF(ABS(DY1) > TB_TOL *ABS(DX1)) THEN
             X1 = XBOX(2)
             Y1 = A3*X1+B3
          ELSE
             ! ASSUME THE LINES ARE NOT NEARLY PARALLEL
             A1 = DY1/DX1
             B1 = YBOX(2) - A1*XBOX(2)

             X1 = (B1-B3)/(A3 - A1)

             Y1 = A1*X1+B1 
          END IF


          IF(ABS(DY2) > TB_TOL *ABS(DX2)) THEN
             X2 = XBOX(2)
             Y2 = A3*X2+B3
          ELSE
             ! ASSUME THE LINES ARE NOT NEARLY PARALLEL

             A2 = DY2/DX2
             B2 = YBOX(2) - A2*XBOX(2)
             X2 = (B2-B3)/(A3 - A2)

             Y2 = A2*X2+B2 
          END IF

       END IF

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------       

       ! LENGTH OF SIDES
       CALL ARC(XBOX(1),YBOX(1),XBOX(2),YBOX(2),DR1)
       CALL ARC(XBOX(2),YBOX(2),XBOX(3),YBOX(3),DR2)
       CALL ARC(X1,Y1,X2,Y2,DR3)

       ! SIDE ONE SEGMENTS
       IF (DR1 > SMALL) THEN
          CALL ARC(XBOX(1),YBOX(1), X1,Y1, DRA)
          DRB = DR1 - DRA
       ELSE
          DR1 = 1.0_SP
          DRA = 0.5_SP
          DRB = 0.5_SP
       END IF

       ! SIDE TWO SEGMENTS
       IF(DR2 > SMALL) THEN
          CALL ARC(XBOX(2),YBOX(2), X2,Y2, DRC)
          DRD = DR2 - DRC
       ELSE
          DR2 = 1.0_SP
          DRC = 0.5_SP
          DRD = 0.5_SP
       END IF

       IF(DR3 > SMALL) THEN
          CALL ARC(X,Y,X1,Y1,DRE)
          DRF = DR3-DRE
       ELSE
          DR3 = 1.0_SP
          DRE = 0.5_SP
          DRF = 0.5_SP
       END IF

# else
       ! LENGTH OF SIDES
       DR1 = SQRT(DX1**2+ DY1**2) !SIDE ONE
       DR2 = SQRT(DX2**2+ DY2**2) !SIDE TWO

       ! SIDE ONE SEGMENTS
       DRA = SQRT( (XBOX(1) - X1)**2 + (YBOX(1) - Y1)**2)
       DRB = DR1 - DRA

       ! SIDE TWO SEGMENTS
       DRC = SQRT( (XBOX(2) - X2)**2 + (YBOX(2) - Y2)**2)
       DRD = DR2 - DRC


       DR3 = SQRT((Y2-Y1)**2 + (X2-X1)**2)
       DRE = SQRT((Y-Y1)**2 + (X-X1)**2)
       DRF = DR3-DRE
# endif

!       write(ipt,*) "Middle! A3=",A3

!       WRITE(IPT,*) "X1=",X1,"; YT=",Y1


!       WRITE(IPT,*) "XB=",X2,"; YB=",Y2


!       write(ipt,*) "DR1=",DR1

!       write(ipt,*) "DR2=",DR2
!       write(ipt,*) "DR3=",DR3

       wghts(1) = (DRF*DRB) / (DR1*DR3)
       wghts(2) = (DRF*DRA) / (DR1*DR3) + (DRE*DRD) / (DR2*DR3)
       wghts(3) = (DRE*DRC) / (DR2*DR3)
       wghts(4) = 0.0_SP

    ELSE

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------       

       CALL ARC(XBOX(1),YBOX(1),X,Y,DR3)
       CALL ARC(XBOX(3),YBOX(3),X,Y,DR1)

# else
       DR3 = SQRT( (XBOX(1) - X)**2 + (YBOX(1) - Y)**2)

       DR1 = SQRT( (XBOX(3) - X)**2 + (YBOX(3) - Y)**2)

# endif
       DRA = (DR1 +DR3)
       ! NORMALIZE TO ONE

       wghts(1) = DR1/DRA
       wghts(2) = 0.0_SP
       wghts(3) = DR3/DRA
       wghts(4) = 0.0_SP


    END IF


  END SUBROUTINE INTERP_NEG
    

  subroutine interp_bilinear_A(Zin,WEIGHTS,zout)
    implicit none
    TYPE(INTERP_WEIGHTS), INTENT(IN) :: WEIGHTS
    real(SP), ALLOCATABLE ,TARGET, INTENT(IN):: zin(:,:)
    real(SP), ALLOCATABLE,TARGET, INTENT(INOUT) :: Zout(:)


    REAL(SP), POINTER :: ZinP(:,:),ZoutP(:)

    IF (allocated(Zin)) THEN
       ZinP => Zin
    ELSE
       CALL FATAL_ERROR &
            &("INTERP: ZIN IS NOT ALLOCATED")
    END IF

    IF (allocated(Zout)) THEN
       ZoutP => Zout   !(lbound(zout):ubound(zout))
    ELSE
       CALL FATAL_ERROR &
            &("INTERP: ZOUT IS NOT ALLOCATED")
    END IF

    CALL INTERP_bilinear_P(zinp,Weights,zoutp)


  end subroutine interp_bilinear_A

  subroutine interp_bilinear_P(Zin,WEIGHTS,zout)
    implicit none
    TYPE(INTERP_WEIGHTS), INTENT(IN) :: WEIGHTS
    real(SP), POINTER , INTENT(IN):: zin(:,:)
    real(SP), POINTER, INTENT(INOUT) :: Zout(:)
    TYPE(R2PTS) :: PTW
    integer  :: i,j,x,y, kk
    real(SP) :: tt

    ! WEIGHTS ARE NORMALIZED TO ONE ALREADY: DO NOT DIVIDE BY THE SUM!

    DO J = 1,ubound(Zout,1)

       PTW=WEIGHTS%PTW(J)

       tt=0.0_SP
       do i=1,4
          x=PTW%quad_num(i,1)
          y=PTW%quad_num(i,2)
          IF(X*Y == 0) CYCLE
          tt=Zin(x,y)*PTW%weights(i,1)+tt
       end do

       Zout(J)=tt

    END DO


    IF(ASSOCIATED(WEIGHTS%N_INDEX)) THEN
       ! THERE IS MASKED DATA - USE THE NEIGHBOR INDEX!
       DO I = 1,ubound(WEIGHTS%N_ORDER,1)
          
          J = WEIGHTS%N_ORDER(I)
          tt=0.0_SP
          
          DO kk = 1, WEIGHTS%N_CNT(J)
             tt = tt + Zout(WEIGHTS%N_INDEX(J,KK))
          END DO

          Zout(J) = tt / real(WEIGHTS%N_CNT(J),SP)
       END DO
    END IF


  end subroutine interp_bilinear_P









 
  SUBROUTINE SETUP_INTERP_QUAD_A(Xin,Yin,Xout,Yout,WEIGHTS, rzero)
    IMPLICIT NONE 
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Xin(:)
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Yin(:)
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Xout(:)
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Yout(:)
    REAL(SP), OPTIONAL :: rzero

    TYPE(INTERP_WEIGHTS), INTENT(OUT) :: WEIGHTS

    REAL(SP), POINTER :: XinP(:)
    REAL(SP), POINTER :: YinP(:)
    REAL(SP), POINTER :: XoutP(:)
    REAL(SP), POINTER :: YoutP(:)

    NULLIFY(XinP,YinP,XoutP,YoutP)


    IF(ALLOCATED(Xin)) XinP => Xin

    IF(ALLOCATED(Yin)) YinP => Yin

    IF(ALLOCATED(Xout)) XoutP => Xout

    IF(ALLOCATED(Yout)) YoutP => Yout

    IF(PRESENT(rzero)) THEN
       CALL SETUP_INTERP_QUAD_P(XinP,YinP,XoutP,YoutP,WEIGHTS,rzero)
    ELSE
       CALL SETUP_INTERP_QUAD_P(XinP,YinP,XoutP,YoutP,WEIGHTS)
    END IF
  END SUBROUTINE SETUP_INTERP_QUAD_A



  SUBROUTINE SETUP_INTERP_QUAD_P(Xin,Yin,Xout,Yout,WEIGHTS,rzero)
    IMPLICIT NONE 
    REAL(SP), POINTER, INTENT(IN) :: Xin(:)
    REAL(SP), POINTER, INTENT(IN) :: Yin(:)
    REAL(SP), POINTER, INTENT(IN) :: Xout(:)
    REAL(SP), POINTER, INTENT(IN) :: Yout(:)
    REAL(SP), OPTIONAL :: rzero

    REAL(SP), ALLOCATABLE :: Xvec(:),Yvec(:)

    TYPE(INTERP_WEIGHTS), INTENT(OUT) :: WEIGHTS

    INTEGER :: I, status, lb, ub

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
    CALL FATAL_ERROR("SETUP_INTERP_QUAD - does not run in spherical!" )
# endif

    ! CHECK INPUT ALLOCATION STATUS
    if(.not.associated(Xin)) CALL FATAL_ERROR&
         &("SETUP_INTERP: INPUT ARGUMENTS MUST BE ALLOCATED!")
    if(.not.associated(Yin)) CALL FATAL_ERROR&
         &("SETUP_INTERP: INPUT ARGUMENTS MUST BE ALLOCATED!")
    if(.not.associated(Xout)) CALL FATAL_ERROR&
         &("SETUP_INTERP: INPUT ARGUMENTS MUST BE ALLOCATED!")
    if(.not.associated(Yout)) CALL FATAL_ERROR&
         &("SETUP_INTERP: INPUT ARGUMENTS MUST BE ALLOCATED!")
    

    ! GET THE DIMENSIONS
 
    WEIGHTS%Nout = ubound(Xout,1)
    WEIGHTS%Nin  = ubound(Xin,1)

    ALLOCATE(WEIGHTS%PTW(WEIGHTS%Nout), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR("SETUP_INTERP: COULD NOT ALLOCATE SPACE")

    ALLOCATE(WEIGHTS%INDEX(WEIGHTS%Nin), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR("SETUP_INTERP: COULD NOT ALLOCATE SPACE")


    CALL SORTRX(WEIGHTS%Nin,XIN,WEIGHTS%INDEX)


    ALLOCATE(XVEC(WEIGHTS%Nin), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR&
         & ("SETUP_INTERP: COULD NOT ALLOCATE SPACE")


    ALLOCATE(YVEC(WEIGHTS%Nin), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR&
         & ("SETUP_INTERP: COULD NOT ALLOCATE SPACE")
    

    DO I = 1,WEIGHTS%Nin
       Xvec(I) = XIN(WEIGHTS%INDEX(I))
       Yvec(I) = YIN(WEIGHTS%INDEX(I))
    END DO
    

    IF(PRESENT(RZERO)) THEN
       DO I = 1, WEIGHTS%Nout
          CALL GEN_WTS(Xvec,Yvec,Xout(I),Yout(I),WEIGHTS%INDEX,WEIGHTS&
               &%PTW(I),rzero)
       END DO
    ELSE
       DO I = 1, WEIGHTS%Nout
          CALL GEN_WTS(Xvec,Yvec,Xout(I),Yout(I),WEIGHTS%INDEX,WEIGHTS&
               &%PTW(I))
       END DO
    END IF
       

  END SUBROUTINE SETUP_INTERP_QUAD_P


  ! THIS CODE IS SETUP TO WORK FOR REGULAR GRIDDED DATA
  ! TO USE IT ON OBSERVATION DATA IT WOULD NEED TO BE MODIFIED TO
  ! HANDLE STRANGE CASES WHERE THERE ARE CLUSTERS OF DATA POINTS
  SUBROUTINE GEN_WTS(Xvec,Yvec,Xout,Yout,INDX,PTW,rzero)
    IMPLICIT NONE
    real(SP), INTENT(in) :: Xout,Yout
    real(SP), ALLOCATABLE, INTENT(IN) :: Xvec(:),Yvec(:)
    INTEGER, ALLOCATABLE, INTENT(IN) :: INDX(:)
    TYPE(R2PTS), INTENT(OUT) :: PTW
    REAL(SP), OPTIONAL :: rzero

    INTEGER :: I,Q,K,NB,J
    
    Real(SP) :: DELX, DELY, D, SUMWGHT

    REAL(SP) :: D_MAX(4)
    REAL(SP) :: DIST(4,2)
    INTEGER :: NDM

    NB = ubound(Xvec,1)

    PTW%WEIGHTS=0
    PTW%QUAD_NUM=0

    D_MAX = HUGE(D)
    DIST  = HUGE(D)

    NDM=1

    k = 0
    DO I=1,NB

       DELX = Xvec(I) - Xout

       ! SINCE XVEC IS IN INCREASING ORDER, ONCE DELX IS GREATER THAN
       ! SEARCH WE ARE DONE
       If (DELX.GT.SEARCH) EXIT

       If (ABS(DELX).LE.SEARCH) Then
          DELY = Yvec(I) - Yout
          D = SQRT(DELX**2+DELY**2)
          
          ! WE FOUND A POINT IN THE SEARCH RADIUS
          If (D.LE.SEARCH) Then
             
             k = k +1
             ! DEAL WITH D == 0.0 here?
             Q=QUADRANT(DELX,DELY)

             ! REPLACE THAT VALUE IF IT IS CLOSER THAN THE LAST ONE
             IF(D .LT. D_MAX(Q))THEN
                ! RECORD THE DISTANCE AND UPDATE THE MAX FOR THAT QUADRANT
                NDM= MAXLOC(DIST(Q,:),1)
                DIST(Q,NDM) = D
                D_MAX(Q) = MAXVAL(DIST(Q,:))
                
                ! RECORD THE INDEX IN QUAD_NUM
                PTW%QUAD_NUM(Q,NDM) = I
             END IF
          End If
       End If
    END DO

    if (k == 0) call warning("FOUND NO Values in interpolation serach radius!")

    IF (PRESENT(RZERO)) THEN
       DIST = Dist/rzero
    END IF


    DO I = 1,4
       DO K =1 , 2
          
          IF (DIST(I,K) < HUGE(D)) THEN
             PTW%WEIGHTS(I,K) = 1/(.5_SP+DIST(I,K)**2)
          END IF
       END DO
    END DO
    
    ! NORMALIZE THE TOTAL WEIGHT TO ONE
    SUMWGHT = SUM(PTW%WEIGHTS)
    PTW%WEIGHTS = PTW%WEIGHTS / SUMWGHT

  END SUBROUTINE GEN_WTS


  INTEGER FUNCTION QUADRANT(DX,DY)
    IMPLICIT NONE
    REAL(SP), INTENT(IN):: DX,DY
    INTEGER :: IBIT1, IBIT2
    
    IBIT1 = (INT(SIGN(REAL(1.0,SP),DX))+1) / 2
    IBIT2 = (INT(SIGN(REAL(1.0,SP),DY))+1) / 2
    !    quadrant number:
    QUADRANT = 2 * IBIT1 + IBIT2 + 1
  END FUNCTION QUADRANT
  

!*****************************************************************
!C Here's a hybrid QuickSort I wrote a number of years ago.  It's
!C based on suggestions in Knuth, Volume 3, and performs much better
!C than a pure QuickSort on short or partially ordered input arrays.  

  SUBROUTINE SORTRX(N,DATA,INDX)
!C===================================================================
!C
!C     SORTRX -- SORT, Real input, indx output
!C
!C
!C     Input:  N     INTEGER
!C             DATA  REAL
!C
!C     Output: INDX INTEGER (DIMENSION N)
!C
!C This routine performs an in-memory sort of the first N elements of
!C array DATA, returning into array INDX the indices of elements of
!C DATA arranged in ascending order.  Thus,
!C
!C    DATA(INDX(1)) will be the smallest number in array DATA;
!C    DATA(INDX(N)) will be the largest number in DATA.
!C
!C The original data is not physically rearranged.  The original order
!C of equal input values is not necessarily preserved.
!C
!C===================================================================
!C
!C SORTRX uses a hybrid QuickSort algorithm, based on several
!C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
!C "pivot key" [my term] for dividing each subsequence is chosen to be
!C the median of the first, last, and middle values of the subsequence;
!C and the QuickSort is cut off when a subsequence has 9 or fewer
!C elements, and a straight insertion sort of the entire array is done
!C at the end.  The result is comparable to a pure insertion sort for
!C very short arrays, and very fast for very large arrays (of order 12
!C micro-sec/element on the 3081K for arrays of 10K elements).  It is
!C also not subject to the poor performance of the pure QuickSort on
!C partially ordered data.
!C
!C Created:  15 Jul 1986  Len Moss
!C
!C===================================================================
    
    INTEGER, INTENT(IN)::  N
    REAL(SP)    :: DATA(N)
    INTEGER :: INDX(N)
    
    INTEGER  :: LSTK(31),RSTK(31),ISTK
    INTEGER  :: L,R,I,J,P,INDXP,INDXT
    REAL     :: DATAP
    
!C     QuickSort Cutoff
!C
!C     Quit QuickSort-ing when a subsequence contains M or fewer
!C     elements and finish off at end with straight insertion sort.
!C     According to Knuth, V.3, the optimum value of M is around 9.
 
    INTEGER, PARAMETER ::   M = 9
    
!C===================================================================
!C
!C     Make initial guess for INDX
 
    DO 50 I=1,N
       INDX(I)=I
50     CONTINUE
 
!C     If array is short, skip QuickSort and go directly to
!C     the straight insertion sort.
 
       IF (N.LE.M) GOTO 900
       
!C===================================================================
!C
!C     QuickSort
!C
!C     The "Qn:"s correspond roughly to steps in Algorithm Q,
!C     Knuth, V.3, PP.116-117, modified to select the median
!C     of the first, last, and middle elements as the "pivot
!C     key" (in Knuth's notation, "K").  Also modified to leave
!C     data in place and produce an INDX array.  To simplify
!C     comments, let DATA[I]=DATA(INDX(I)).
 
!C Q1: Initialize
       ISTK=0
       L=1
       R=N
       
200    CONTINUE
 
!C Q2: Sort the subsequence DATA[L]..DATA[R].
!C
!C     At this point, DATA[l] &lt;= DATA[m] &lt;= DATA[r] for all l &lt; L,
!C     r &gt; R, and L &lt;= m &lt;= R.  (First time through, there is no
!C     DATA for l &lt; L or r &gt; R.)
 
       I=L
       J=R
 
!C Q2.5: Select pivot key
!C
!C     Let the pivot, P, be the midpoint of this subsequence,
!C     P=(L+R)/2; then rearrange INDX(L), INDX(P), and INDX(R)
!C     so the corresponding DATA values are in increasing order.
!C     The pivot key, DATAP, is then DATA[P].
 
       P=(L+R)/2
       INDXP=INDX(P)
       DATAP=DATA(INDXP)
       
       IF (DATA(INDX(L)) .GT. DATAP) THEN
          INDX(P)=INDX(L)
          INDX(L)=INDXP
          INDXP=INDX(P)
          DATAP=DATA(INDXP)
       ENDIF
       
       IF (DATAP .GT. DATA(INDX(R))) THEN
          IF (DATA(INDX(L)) .GT. DATA(INDX(R))) THEN
             INDX(P)=INDX(L)
             INDX(L)=INDX(R)
          ELSE
             INDX(P)=INDX(R)
          ENDIF
          INDX(R)=INDXP
          INDXP=INDX(P)
          DATAP=DATA(INDXP)
       ENDIF
 
!C     Now we swap values between the right and left sides and/or
!C     move DATAP until all smaller values are on the left and all
!C     larger values are on the right.  Neither the left or right
!C     side will be internally ordered yet; however, DATAP will be
!C     in its final position.
 
300    CONTINUE
 
!C Q3: Search for datum on left &gt;= DATAP
!C
!C     At this point, DATA[L] &lt;= DATAP.  We can therefore start scanning
!C     up from L, looking for a value &gt;= DATAP (this scan is guaranteed
!C     to terminate since we initially placed DATAP near the middle of
!C     the subsequence).
 
       I=I+1
       IF (DATA(INDX(I)).LT.DATAP) GOTO 300
       
400    CONTINUE
       
!C Q4: Search for datum on right &lt;= DATAP
!C
!C     At this point, DATA[R] &gt;= DATAP.  We can therefore start scanning
!C     down from R, looking for a value &lt;= DATAP (this scan is guaranteed
!C     to terminate since we initially placed DATAP near the middle of
!C     the subsequence).
 
       J=J-1
       IF (DATA(INDX(J)).GT.DATAP) GOTO 400
 
!C Q5: Have the two scans collided?
 
       IF (I.LT.J) THEN
 
!C Q6: No, interchange DATA[I] &lt;--&gt; DATA[J] and continue
          
          INDXT=INDX(I)
          INDX(I)=INDX(J)
          INDX(J)=INDXT
          GOTO 300
       ELSE
          
!C Q7: Yes, select next subsequence to sort
!C
!C    At this point, I &gt;= J and DATA[l] &lt;= DATA[I] == DATAP &lt;= DATA[r],
!C     for all L &lt;= l &lt; I and J &lt; r &lt;= R.  If both subsequences are
!C     more than M elements long, push the longer one on the stack and
!C     go back to QuickSort the shorter; if only one is more than M
!C     elements long, go back and QuickSort it; otherwise, pop a
!C     subsequence off the stack and QuickSort it.
 
          IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
             ISTK=ISTK+1
             LSTK(ISTK)=J+1
             RSTK(ISTK)=R
             R=I-1
          ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
             ISTK=ISTK+1
             LSTK(ISTK)=L
             RSTK(ISTK)=I-1
             L=J+1
          ELSE IF (R-J .GT. M) THEN
             L=J+1
          ELSE IF (I-L .GT. M) THEN
             R=I-1
          ELSE
!C Q8: Pop the stack, or terminate QuickSort if empty
             IF (ISTK.LT.1) GOTO 900
             L=LSTK(ISTK)
             R=RSTK(ISTK)
             ISTK=ISTK-1
          ENDIF
          GOTO 200
       ENDIF
 
900    CONTINUE
 
!C===================================================================
!C
!C Q9: Straight Insertion sort
 
       DO 950 I=2,N
          !         if(mod(I,10000).eq.0) print*,'i=',i
          IF (DATA(INDX(I-1)) .GT. DATA(INDX(I))) THEN
             INDXP=INDX(I)
             DATAP=DATA(INDXP)
             P=I-1
920          CONTINUE
             INDX(P+1) = INDX(P)
             P=P-1
             IF (P.GT.0) THEN
                IF (DATA(INDX(P)).GT.DATAP) GOTO 920
             ENDIF
             INDX(P+1) = INDXP
          ENDIF
950       CONTINUE
          
!C===================================================================
!C
!C     All done
          RETURN
!       END DO
!    END DO
  END SUBROUTINE SORTRX
!******************************


  SUBROUTINE INTERP_QUAD_A(zin,Weights,zout)
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Zin(:)
    TYPE(INTERP_WEIGHTS), INTENT(in) :: WEIGHTS
    REAL(SP), ALLOCATABLE, TARGET :: Zout(:)

    REAL(SP), POINTER :: ZinP(:),ZoutP(:)

    IF (allocated(Zin)) THEN
       ZinP => Zin
    ELSE
       CALL FATAL_ERROR &
          &("INTERP: ZIN IS NOT ALLOCATED")
    END IF

    IF (allocated(Zout)) THEN
       ZoutP => Zout   !(lbound(zout):ubound(zout))
    ELSE
       CALL FATAL_ERROR &
          &("INTERP: ZOUT IS NOT ALLOCATED")
    END IF
 
    CALL INTERP_QUAD_P(zinp,Weights,zoutp)

  END SUBROUTINE INTERP_QUAD_A
  
  SUBROUTINE INTERP_QUAD_P(zin,Weights,zout)
     REAL(SP), POINTER :: Zin(:)
     TYPE(INTERP_WEIGHTS), INTENT(in) :: WEIGHTS
     REAL(SP), POINTER :: Zout(:)

     REAL(SP), ALLOCATABLE :: Zvec(:)
     
     integer :: I,status, msze,nsze, lb, ub

     IF (.NOT. ASSOCIATED(ZIN))CALL FATAL_ERROR &
          &("INTERP: ZIN IS NOT ALLOCATED")

     ! CAN'T CHECK ALLOCATION STATUS OF INTENT OUT VARIABLE!

     
     IF (Weights%Nin /= ubound(zin,1)) CALL FATAL_ERROR &
          &("INTERP: THE DATA INPUT SIZE DOES NOT MATCH THE WEIGHT!")

     ! MAKE VECTOR IN RIGHT ORDER...
     ALLOCATE(Zvec(0:Weights%Nin),stat=status)
     if(status /= 0) CALL FATAL_ERROR("INPTERP: COULD NOT ALLOCATE SPACE")
     Zvec = 0.0_SP

     DO I = 1,WEIGHTS%Nin
        Zvec(I) = Zin(WEIGHTS%INDEX(I))
     END DO

     DO I = 1,WEIGHTS%Nout

        CALL INTERP_WEIGH(ZVEC,WEIGHTS%PTW(I),Zout(I))

     END DO

     !DEALLOCATE IS AUTOMATIC FOR ALLOCATABLES

   end subroutine INTERP_QUAD_P

! VECTORIZE THIS - THERE IS NO REASON TO USE AN ARRAY TO STORE THE WEIGHTS!
   subroutine interp_weigh(Zvec,PTW,zval)
     implicit none
     real(SP), INTENT(OUT) :: ZVAL
     TYPE(R2PTS), INTENT(IN) :: PTW
     real(SP), ALLOCATABLE , INTENT(IN):: zvec(:)
     integer  :: i,j,n

     ! WEIGHTS ARE NORMALIZED TO ONE ALREADY: DO NOT DIVIDE BY THE SUM!

     zval=0.0_SP
     do i=1,4
        do j=1,2
           n=PTW%quad_num(i,j)
           zval=Zvec(n)*PTW%weights(i,j)+zval
        end do
     end do

   end subroutine interp_weigh
   

END MODULE MOD_INTERP
