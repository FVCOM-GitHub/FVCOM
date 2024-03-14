  SUBROUTINE SET_W_OBC
# if defined(NH)
  USE MOD_NESTING
  USE NON_HYDRO
  USE MOD_OBCS
  IMPLICIT NONE
  INTEGER I, J, K
  INTEGER JN, II, I1, I2

!---> Modified by Dr. Lai 2021-01-15
! The original version: use hyw_nest from nesting file as w4zt (with nesting)
!                       use wts and elevation to calculate w4zt
! The update version: directly use w4zt 
!%  IF(NESTING_ON)THEN
    CALL SET_VAR(intTime,HYW=W4ZT)

!    DO I=1, M
!      IF(ISONB(I)==2) THEN
!        W4ZT(I,1:KB) = HYW_NEST(I,1:KB)
!      ENDIF
!    ENDDO

!%  ELSE

!%    DO J = 1, IBCN(1)
!%      JN = OBC_LST(1,J)
!%      II = I_OBC_N(JN)

!%      DO K=1, KB
!%        W4ZT(II,K) = WTS(II,K) + (1.0_SP+Z(II,K))*(EL(II)-ET(II))/DTI
!%      ENDDO
!%    END DO
!%    DO J = 1,IBCN(4)
!%      JN = OBC_LST(4,J)
!%      I1 = I_OBC_N(JN)
!%      I2 = NEXT_OBC(JN)

!%      DO K=1, KB
!%        W4ZT(I1,K) = W4ZT(I2,K)
!%      ENDDO
!%    END DO

!%  ENDIF

!<--- Modified by Dr. Lai 2021-01-15

  RETURN
# endif
  END SUBROUTINE SET_W_OBC

  SUBROUTINE SET_UV_OBC
# if defined(NH)
  USE MOD_NESTING
  USE NON_HYDRO
  IMPLICIT NONE
  INTEGER I

  IF(NESTING_ON)THEN
    CALL SET_VAR(intTime,U=UF)
    CALL SET_VAR(intTime,V=VF)
  END IF

  RETURN
# endif
  END SUBROUTINE SET_UV_OBC

  SUBROUTINE SET_NH_OBC
# if defined(NH)
  USE MOD_NESTING
  USE NON_HYDRO
  IMPLICIT NONE
  INTEGER I

!  IF(NESTING_ON)THEN
!    DO I=1,nt_nest
!      NHQDRX(NID_NEST(I),:) = 0.0_SP
!      NHQDRY(NID_NEST(I),:) = 0.0_SP
!      NHQ2DX(NID_NEST(I))   = 0.0_SP
!      NHQ2DY(NID_NEST(I))   = 0.0_SP
!    END DO
!  END IF

  RETURN
# endif
  END SUBROUTINE SET_NH_OBC
