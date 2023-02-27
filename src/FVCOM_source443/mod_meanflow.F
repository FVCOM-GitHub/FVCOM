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

MODULE MOD_MEANFLOW
   USE ALL_VARS
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   USE MOD_PREC
   USE MOD_TYPES
   IMPLICIT NONE
   SAVE
   INTEGER              :: MF_RST_STCNT
   INTEGER              :: INMF,INTCELL,INTNODE,INTELEL,INTUV
   INTEGER              :: nmfcell_GL, nmfcell, nmfcell_i
   INTEGER, ALLOCATABLE :: MF_GL2LOC(:)
   INTEGER, ALLOCATABLE :: I_MFCELL_GL(:),I_MFCELL_N(:)
   REAL(SP),ALLOCATABLE :: DMFQDIS(:,:),MFQDIS(:),MFDIST(:,:)
   REAL(SP),ALLOCATABLE :: ANGLEMF(:),MFAREA(:),VLCTYMF(:)
   TYPE(BC)             :: MF_TM           !!TIME MAP FOR MEAN FLOW DATA
   REAL(SP),ALLOCATABLE :: RDISMF(:,:)
   INTEGER ,ALLOCATABLE :: NODE_MFCELL(:,:)

   CONTAINS

! we still need to consider the case in which MEAN FLOW bring in/take out T & S
!==============================================================================|
!  READ IN MEAN FLOW OPEN BOUNDARY FLUX (m^3/s^1) TIME SERIES                  |
!==============================================================================|

   SUBROUTINE READ_MEANFLOW

!------------------------------------------------------------------------------!
     INTEGER              :: k,i,j,i1,i2,i3,ii,NCNT,itemp,IERR
     INTEGER, ALLOCATABLE :: temp1(:),temp2(:)
     REAL(SP),ALLOCATABLE :: RTEMP1(:,:),RTEMP2(:,:)
     REAL(SP)             :: ttemp

     REWIND(INMF)
     READ(INMF,*) nmfcell_GL

     nmfcell_i = 0
     nmfcell   = 0
  IF (nmfcell_GL > 0) THEN

     ALLOCATE(I_MFCELL_GL(nmfcell_GL))
     DO I=1,nmfcell_GL
        READ(INMF,*)I_MFCELL_GL(I)
     ENDDO

!----Read in Mean Flow Flux Vertical Distribution---------------------
     ALLOCATE(RTEMP1(nmfcell_GL,KBM1))
     DO I = 1, nmfcell_GL
       READ(INMF,*) J,(RTEMP1(I,K),K = 1,KBM1)
     END DO

!----Read in Time Dependent DataSets ---------------------------------
       READ(INMF,*) itemp
       MF_TM%NTIMES = itemp
       MF_TM%LABEL  = "open boundary mean flow flux"
       ALLOCATE(MF_TM%TIMES(itemp))
       ALLOCATE(RTEMP2(nmfcell_GL,itemp))
       DO I = 1, itemp
         READ(INMF,*) ttemp
         MF_TM%TIMES(I) = ttemp
         READ(INMF,*) (RTEMP2(J,I),J = 1,nmfcell_GL)
!--------------------------------Jianzhong----------------------------
         IF(MSR)WRITE(IPT,*)MAXVAL(RTEMP2(1:NMFCELL_GL,I))&
              &,MAXLOC(RTEMP2(1:NMFCELL_GL,I)) ,MINVAL(RTEMP2(1:NMFCELL_GL,I))&
              &,MINLOC(RTEMP2(1:NMFCELL_GL,I)) 
!---------------------------------------------------------------------
 !        WRITE(IOPRT,*) ttemp
 !        WRITE(IOPRT,*) (RTEMP2(J,I),J = 1,nmfcell_GL)
       END DO
       CLOSE(INMF)

!
!---Map to Local Domain----------------------------------------

     IF(SERIAL)THEN
       nmfcell_i = nmfcell_GL
       nmfcell   = nmfcell_GL
       ALLOCATE(I_MFCELL_N(nmfcell))
       I_MFCELL_N = I_MFCELL_GL
       ALLOCATE(MFDIST(nmfcell,kbm1))
       MFDIST = RTEMP1
       ALLOCATE(DMFQDIS(nmfcell,MF_TM%NTIMES))
       DMFQDIS = RTEMP2
     END IF

#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       ALLOCATE(TEMP1(nmfcell_GL))
       ALLOCATE(TEMP2(nmfcell_GL))
       NCNT = 0
       DO I=1,nmfcell_GL
!         I1=ELID_X(I_MFCELL_GL(I))
         I1=ELID(I_MFCELL_GL(I))
	 IF(I1 /= 0)THEN
	   NCNT = NCNT + 1
	   TEMP1(NCNT) = I1
	   TEMP2(NCNT) = I
	 END IF
       END DO
       nmfcell_i = NCNT

       DO I=1,nmfcell_GL
         I1=ELID_X(I_MFCELL_GL(I))
         I2=ELID(I_MFCELL_GL(I))
	 IF(I1 /= 0 .and. I1 /= I2)THEN
	   NCNT = NCNT + 1
	   TEMP1(NCNT) = I1
	   TEMP2(NCNT) = I
	 END IF
       END DO
       nmfcell = NCNT
       IF(nmfcell > 0)THEN
         ALLOCATE(I_MFCELL_N(nmfcell),MF_GL2LOC(nmfcell))
	 I_MFCELL_N(1:nmfcell) = TEMP1(1:nmfcell)
         MF_GL2LOC (1:nmfcell) = TEMP2(1:nmfcell)
       END IF

!       do i = 1,nmfcell
!          write(ipt_p,*)I_MFCELL_N(I), MF_GL2LOC(I)
!       end do



       DEALLOCATE(TEMP1,TEMP2)

       IF(nmfcell > 0)THEN  
          ALLOCATE(MFDIST(nmfcell,kbm1))
          DO I=1,nmfcell
          DO K=1,KBM1
             MFDIST(I,K) = RTEMP1(MF_GL2LOC(I),K)
          END DO
          END DO
       
          ALLOCATE(DMFQDIS(nmfcell,MF_TM%NTIMES))
          DO I=1,MF_TM%NTIMES
             DMFQDIS(1:nmfcell,I)=RTEMP2(MF_GL2LOC(1:nmfcell),I)
          END DO
       END IF

     END IF
#   endif

     DEALLOCATE(RTEMP1,RTEMP2)

  ELSE  ! if statement end for nmfcell_GL > 0
    close(INMF)
  END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
#  endif
!--------------------------------------Jianzhong----------------------
   WRITE(IPT,*)'NMFCELL_I=',NMFCELL_I,'NMFCELL=',NMFCELL,'IN THREAD:',MYID
!---------------------------------------------------------------------

   RETURN
   END SUBROUTINE READ_MEANFLOW
!==============================================================================|


!==============================================================================|
!  SET METRICS FOR MEAN FLOW BOUNDARY CONDITIONS       			       |
!==============================================================================|

   SUBROUTINE SET_BNDRY_MEANFLOW     

!------------------------------------------------------------------------------!

   USE BCS
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
   USE MOD_OBCS

   IMPLICIT NONE
   REAL(DP)  DX12,DY12,ATMP1,HTMP
   INTEGER I,J,I1,I2,J1,J2,II,ITMP,JTMP
# if defined (SPHERICAL)
   REAL(DP) X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE
# endif
!------------------------------------------------------------------------------!

   IF(nmfcell > 0)THEN

     ALLOCATE(ANGLEMF(nmfcell),MFAREA(nmfcell),VLCTYMF(nmfcell),MFQDIS(nmfcell))
     ALLOCATE(NODE_MFCELL(nmfcell,2),RDISMF(nmfcell,2))

     DO I=1,nmfcell
       II=I_MFCELL_N(I)
       IF(I <= nmfcell_i .and. ISBCE(II) /= 2) THEN
          PRINT*, 'NO.',I,'MEAN FLOW CELL'
          PRINT*, 'IS NOT A OPEN BOUNDARY ONE'
          CALL PSTOP
       END IF
       ITMP=0
       DO J=1,3
         IF(NBE(II,J) == 0 .and. ISONB(nv(II,J)) /= 2) THEN
           JTMP=J
           ITMP=ITMP+1
         END IF
       END DO
       IF(ITMP /= 1) THEN
         PRINT*, 'NO OPEN BOUNDARY OR MORE THAN ONE OPEN BOUNDARY'
         PRINT*, 'IN NO.',I,'MEAN FLOW CELL'
         CALL PSTOP
       END IF
       J1=JTMP+1-INT((JTMP+1)/4)*3
       J2=JTMP+2-INT((JTMP+2)/4)*3
       I1=NV(II,J1)
       I2=NV(II,J2)
         
       NODE_MFCELL(I,1)=I1
       NODE_MFCELL(I,2)=I2

       HTMP=0.5_SP*(H(I1)+H(I2))     ! may be a problem here, should be replaced dy D
       DY12=VY(I1)-VY(I2)
#      if defined (SPHERICAL)
       X1_DP = VX(I2)
       Y1_DP = VY(I2)
       X2_DP = VX(I1)
       Y2_DP = VY(I1)
       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
       DX12 = SIDE

       DY12 = TPI*DY12
#      else
       DX12=VX(I1)-VX(I2)
#      endif
       ATMP1=ATAN2(DY12,DX12)
       MFAREA(I)=SQRT(DX12**2+DY12**2)*HTMP    ! for spherical coordinates is Phthagolean Theorem still valid?
       ANGLEMF(I)=ATMP1+3.1415927/2.0
       RDISMF(I,1)=ART1(I1)/(ART1(I1)+ART1(I2))
       RDISMF(I,2)=ART1(I2)/(ART1(I1)+ART1(I2))
     END DO
   END IF

   RETURN
   END SUBROUTINE SET_BNDRY_MEANFLOW
!==============================================================================|

!==============================================================================|
!  INTERPOLATION MEAN FLOW OPEN BOUNDARY FLUX (m^3/s^1) TIME SERIES            |
!==============================================================================|

   SUBROUTINE BCOND_MEANFLOW
!
!------------------------------------------------------------------------------!
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   
   INTEGER  L1,L2,IERR,II
   REAL(SP) :: FACT,UFACT
   REAL(SP) :: THOUR

   THOUR = DTI*FLOAT(IINT)/3600.0
   IF(nmfcell > 0)THEN
     CALL BRACKET(MF_TM,THOUR,L1,L2,FACT,UFACT,IERR)
     MFQDIS(:) = UFACT*DMFQDIS(:,L1) + FACT*DMFQDIS(:,L2)
     MFQDIS    = MFQDIS*RAMP
   END IF

   RETURN
   END SUBROUTINE BCOND_MEANFLOW

   SUBROUTINE BRACKET(TMAP,STIME,L1,L2,FACT,BACT,IERR)             
!==============================================================================|
!  DETERMINE DATA INTERVAL IN WHICH CURRENT TIME LIES                          |
!									       | 
!  L1:  DATA INTERVAL PROCEEDING TIME				               |
!  L2:  DATA INTERVAL AFTER TIME                                               |
!  FACT: LINEAR INTERPOLATION COEFFICIENT (0->1)                               |
!     FACT = .5  : STIME LIES EXACTLY BETWEEN TWO DATA TIMES                   |
!     FACT = 1.  : STIME OCCURS AT SECOND DATA TIME                            |
!  BACT  = 1.-FACT
!  IERR: RETURNS INTEGER ERROR                                                 |
!     IERR = 0   : NO ERROR, TIME IS BRACKETED BY DATA TIMES                   |
!     IERR =-1   : STIME PROCEEDS ALL DATA TIMES                               |
!     IERR = 1   : STIME IS GREATER THAN ALL DATA TIMES                        |
!                                                                              |
!  IF STIME PROCEEDS DATA, IERR IS SET TO -1, L1 TO 1, AND FACT TO 0.          !
!  IF STIME SUPERCEEDS DATA, IERR IS SET TO -1, L2 TO LMAX, AND FACT TO 1.     !
!==============================================================================|
   USE MOD_TYPES 
   IMPLICIT NONE
!------------------------------------------------------------------------------!
   TYPE(BC), INTENT(IN)  :: TMAP
   REAL(SP), INTENT(IN)  :: STIME
   INTEGER,  INTENT(OUT) :: L1,L2
   REAL(SP), INTENT(OUT) :: FACT,BACT
   INTEGER,  INTENT(OUT) :: IERR
!------------------------------------------------------------------------------!
   REAL(SP)  T1,T2
   INTEGER I,NTMAX
!==============================================================================|

   NTMAX = TMAP%NTIMES
   IF(STIME < TMAP%TIMES(1))THEN
     FACT = 0.0_SP
     BACT = 1.0_SP
     L1   = 1
     L2   = 1
     IERR = -1
     RETURN
   END IF

   IF(STIME > TMAP%TIMES(NTMAX))THEN
     FACT = 1.0_SP
     BACT = 0.0_SP
     L1   = NTMAX
     L2   = NTMAX
     IERR = 1
     RETURN
   END IF

   IF(NTMAX == 1)THEN
     FACT = 1.0_SP
     BACT = 0.0_SP
     L1   = 1
     L2   = 1
     IERR = 0
     RETURN
   END IF
   

   DO I=2,TMAP%NTIMES
     T1 = TMAP%TIMES(I-1)
     T2 = TMAP%TIMES(I)
     IF(STIME >= T1 .AND. STIME <= T2)THEN  
       L1   = I-1
       L2   = I
       IERR = 0
       FACT = (STIME-T1)/(T2-T1)
       BACT = 1.0_SP-FACT
     END IF
   END DO
   
     
   RETURN
   END SUBROUTINE BRACKET
!==============================================================================|


END MODULE MOD_MEANFLOW
