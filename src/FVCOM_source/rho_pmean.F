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
SUBROUTINE RHO_PMEAN
  !==============================================================================|
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_PAR
  USE SINTER
# if defined (WET_DRY)
  USE MOD_WD
# endif
  IMPLICIT NONE

  INTEGER, PARAMETER :: KBB=601
  INTEGER, PARAMETER :: KBBM1=KBB-1

  REAL(SP) DELTZ
  REAL(SP), DIMENSION(KBB)   :: PHY_Z   !Depth(m) in every standary Z levels 
  REAL(SP), DIMENSION(KBB)   :: RHOZ    !density in standary Z levels 

  REAL(SP), DIMENSION(KBB)   :: RHOA    !density mean in standary Z levals

  REAL(SP), DIMENSION(KBM1)  :: ZM      !Depth (m) in every sigma levels for giving node
  REAL(SP), DIMENSION(KBM1)  :: RHOS      !DENS AT SIGMA LEVELS

  REAL(SP), DIMENSION(KBB)   :: FCOUNT

  INTEGER :: I,K,status,IERR,IB
  REAL(SP) :: Z_TOP, Z_BOTTOM,SBUF

  !========================================
  ! ONLY USED FOR PAR
  REAL(SP), DIMENSION(KBB,NPROCS) :: RHOA_RCV, FCOUNT_RCV
  !========================================

  !--CALCULATE Z-LEVELS TO MAX DEPTH---------------------------------------------|

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: RHO_PMEAN"

  IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! Recalculating mean density in pressure coordinates"

# if defined (WET_DRY)

  IF(SERIAL) THEN
     IF(WETTING_DRYING_ON)THEN
       Z_TOP = maxval(EL(1:M)*ISWETN(1:M))
     ELSE
       Z_TOP = maxval(EL(1:M))
     END IF
     Z_BOTTOM = HMAX+spacing(hmax)
  ELSE
#    if defined(MULTIPROCESSOR)
     IF(WETTING_DRYING_ON)THEN
       SBUF = maxval(EL(1:M)*ISWETN(1:M))
     ELSE
       SBUF = maxval(EL(1:M))
     END IF
     CALL MPI_ALLREDUCE(SBUF,Z_TOP,1,MPI_F,MPI_MAX,MPI_FVCOM_GROUP,IERR)
     Z_BOTTOM = HMAX ! ALREADY GLOBALLY COLLECTED!
# endif
  END IF

# else

  IF(SERIAL) THEN
#   if defined (ICE_EMBEDDING)
     Z_TOP = maxval(EL(1:M))-minval(QTHICE(1:M))   !!yding
#   else
     Z_TOP = maxval(EL(1:M))
#   endif
     Z_BOTTOM = HMAX+spacing(hmax)
  ELSE
#    if defined(MULTIPROCESSOR)
#   if defined (ICE_EMBEDDING)
     SBUF = maxval(EL(1:M))-minval(QTHICE(1:M))   !!yding
#   else
     SBUF = maxval(EL(1:M))
#   endif
     CALL MPI_ALLREDUCE(SBUF,Z_TOP,1,MPI_F,MPI_MAX,MPI_FVCOM_GROUP,IERR)
     Z_BOTTOM = HMAX ! ALREADY GLOBALLY COLLECTED!
# endif
  END IF

# endif  

!  DELTZ=HMAX/FLOAT(KBBM1)
  DELTZ=(Z_BOTTOM + Z_TOP)/FLOAT(KBBM1)

  DO K=1,KBB
     PHY_Z(K)= Z_TOP - FLOAT(K-1)*DELTZ
  END DO

!  PHY_Z(1)=PHY_Z(1)-spacing(PHY_Z(1))
!  PHY_Z(KBB)=PHY_Z(KBB)+spacing(PHY_Z(KBB))

  PHY_Z(1)=nearest(Z_TOP,-1.0_SP) ! Nearest number smaller than Z_TOP

  PHY_Z(KBB)=nearest(-Z_BOTTOM,1.0_SP) ! Nearest number larger than Z_BOTTOM

  
  !--DO THE AVERAGE OVER Z_levels 


  IF(SERIAL) THEN
     RHOA=0.0_SP
     FCOUNT=0.0_SP
     DO I=1,M

        !--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT Z LEVELS-------------------|
#   if defined (ICE_EMBEDDING)
        ZM(1:kbm1)=ZZ(I,1:kbm1)*D(I)+EL(I)-QTHICE(I)    !!yding
#   else        
        ZM(1:kbm1)=ZZ(I,1:kbm1)*D(I)+EL(I)
#   endif
        RHOS = RHO1(I,1:KBM1)
        CALL SINTER_EXTRP_NONE(ZM,RHOS,PHY_Z,RHOZ,KBM1,KBB)

        !--SUM THE DENSITY ACROSS ALL THE NODES AT ZLEVELS ----------!
        DO K=1,KBB       
#   if defined (ICE_EMBEDDING)
           IF(-H(I).LE.PHY_Z(K) .AND. EL(I)-QTHICE(I) .GE. PHY_Z(K)) THEN  !!yding
#   else 
           IF(-H(I).LE.PHY_Z(K) .AND. EL(I) .GE. PHY_Z(K)) THEN
#   endif
              FCOUNT(K) = FCOUNT(K) + 1.0_SP
              RHOA(K)=RHOA(K)+RHOZ(K)
           END IF

        END DO

     END DO
        
     ! TAKE THE AVERAGE
     IF(MINVAL(FCOUNT) .LT. 1.0_SP)THEN
        IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND NO DATA AT DEPTH:",PHY_Z(MINLOC(FCOUNT,1))
        CALL FATAL_ERROR("RHO_PMEAN: In Serial case, found fcount LT 0.0!")
     END IF
     RHOA = RHOA / FCOUNT

  ELSE

#    if defined(MULTIPROCESSOR)

     FCOUNT=0.0_SP
     RHOA=0.0_SP
     RHOA_RCV=0.0_SP
     FCOUNT_RCV=0.0_SP

     DO I=1,M        
        !--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT Z LEVELS-------------------|
#   if defined (ICE_EMBEDDING)
        ZM(1:kbm1)=ZZ(I,1:kbm1)*D(I)+EL(I)-QTHICE(I)   !!yding
#   else
        ZM(1:kbm1)=ZZ(I,1:kbm1)*D(I)+EL(I)
#   endif
        RHOS = RHO1(I,1:KBM1)
        CALL SINTER_EXTRP_NONE(ZM,RHOS,PHY_Z,RHOZ,KBM1,KBB) 

        !--SUM THE DENSITY ACROSS ALL THE NODES AT ZLEVELS ----------!
        DO K=1,KBB       
#   if defined (ICE_EMBEDDING) 
           IF(-H(I).LE.PHY_Z(K) .AND. EL(I)-QTHICE(I) .GE. PHY_Z(K)) THEN   !!yding
#   else
           IF(-H(I).LE.PHY_Z(K) .AND. EL(I) .GE. PHY_Z(K)) THEN
#   endif
              IF(NDE_ID(I) == 0)THEN !!INTERNAL NODE
                 FCOUNT(K) = FCOUNT(K) + 1.0_SP
                 RHOA(K)=RHOA(K)+RHOZ(K)
              ELSE  !!BOUNDARY NODE, ACCUMULATE FRACTION ONLY
                 DO IB = 1,NBN
                    IF(BN_LOC(IB) == I)THEN
                       FCOUNT(K) = FCOUNT(K) + 1.0_SP/FLOAT(BN_MLT(IB))
                       RHOA(K)=RHOA(K)+RHOZ(K)
                    END IF
                 END DO
              END IF
           END IF
        END DO

     END DO

     CALL MPI_ALLGATHER(RHOA,KBB,MPI_F,RHOA_RCV,KBB,MPI_F,MPI_FVCOM_GROUP,IERR)
     CALL MPI_ALLGATHER(FCOUNT,KBB,MPI_F,FCOUNT_RCV,KBB,MPI_F,MPI_FVCOM_GROUP,IERR)
     RHOA = SUM(RHOA_RCV,2)
     FCOUNT = SUM(FCOUNT_RCV,2)
        
!     IF(MINVAL(FCOUNT) .LT. 1.0_SP) THEN
# if !defined (DOUBLE_PRECISION)
     IF(NINT(MINVAL(FCOUNT)) .LT. 1) THEN
# else
     IF(IDNINT(MINVAL(FCOUNT)) .LT. 1) THEN
# endif
        IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND NO DATA AT DEPTH:",PHY_Z(MINLOC(FCOUNT,1))
        CALL FATAL_ERROR("RHO_PMEAN: In Parallel case, found fcount LT 0.0!")
     END IF

     RHOA = RHOA/FCOUNT

#    endif
  END IF


  !--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT SIGMA LEVELS-------------------|


  DO I=1,M
# if defined (WET_DRY)
     IF(WETTING_DRYING_ON.and.ISWETN(i)==0)cycle
# endif
     DO K=1,KBM1
#   if defined (ICE_EMBEDDING)
        ZM(K)=ZZ(I,K)*D(I)+EL(I)-QTHICE(I)       !!yding
#   else
        ZM(K)=ZZ(I,K)*D(I)+EL(I)
#   endif
     END DO
     CALL SINTER_EXTRP_NONE(PHY_Z,RHOA,ZM,RHOS,KBB,KBM1)
     RMEAN1(I,1:KBM1) = RHOS
  END DO

  RMEAN1(:,KB)=0.0_SP


  CALL N2E3D(RMEAN1,RMEAN)


  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: RHO_PMEAN"

  RETURN
END SUBROUTINE RHO_PMEAN
