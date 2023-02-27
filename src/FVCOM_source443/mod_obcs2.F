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

MODULE MOD_OBCS2

# if !defined (SEMI_IMPLICIT)
   USE ALL_VARS
   USE MOD_PREC
   USE MOD_OBCS
   USE MOD_TYPES
   USE BCS

   USE MOD_MEANFLOW

   IMPLICIT NONE
   SAVE

   INTEGER              :: NOBE,NOBCV
   INTEGER,ALLOCATABLE  :: NOBEDGE_LST(:),COBEDGE_LST(:)
   INTEGER,ALLOCATABLE  :: I_OBC_CELL(:),I_OBC_NODE(:),I_OBC_CELL2(:),I_OBC_NODE2(:)
   REAL(SP),ALLOCATABLE :: UATTS(:,:),VATTS(:,:),UTTS(:,:,:),VTTS(:,:,:),ELTTS(:,:)

   REAL(SP),ALLOCATABLE :: ELT(:), ELTF(:), ELRKT(:), ELTDT(:)
   REAL(SP),ALLOCATABLE :: ELP(:), ELPF(:), ELRKP(:)
   REAL(SP),ALLOCATABLE :: UAT(:), VAT(:), UATF(:), VATF(:)
   REAL(SP),ALLOCATABLE :: UAP(:), VAP(:)

   REAL(SP),ALLOCATABLE :: UANT (:), VANT (:), UAN (:), VAN (:)
   REAL(SP),ALLOCATABLE :: UANP (:), VANP (:)

   REAL(SP),ALLOCATABLE :: UT(:,:),   VT(:,:)
   REAL(SP),ALLOCATABLE :: UNT(:,:),  VNT(:,:),  UN(:,:),  VN(:,:)

   REAL(SP),ALLOCATABLE :: UAPF  (:), VAPF  (:)
   REAL(SP),ALLOCATABLE :: UANTF (:), VANTF (:), UANF (:), VANF (:)
   REAL(SP),ALLOCATABLE :: UANPF (:), VANPF (:)
   REAL(SP),ALLOCATABLE :: UARKNT(:), VARKNT(:), UARKN(:), VARKN(:)
   REAL(SP),ALLOCATABLE :: UNTB(:,:), VNTB(:,:), UNB(:,:), VNB(:,:)

   INTEGER              :: ntidecell_GL, ntidecell, ntidecell_i
   INTEGER              :: ntidenode_GL, ntidenode, ntidenode_i
   INTEGER, ALLOCATABLE :: I_TIDENODE_GL(:),I_TIDENODE_N(:)
   INTEGER, ALLOCATABLE :: I_TIDECELL_GL(:),I_TIDECELL_N(:)

   REAL(SP),ALLOCATABLE :: DLTN(:)
   INTEGER :: IOS
   CONTAINS

!=========================================================================|
   SUBROUTINE ALLOC_OBC2_DATA

   IMPLICIT NONE

   ALLOCATE(ELT  (0:ntidenode));                 ELT           = ZERO
   ALLOCATE(ELTF (0:ntidenode));                 ELTF          = ZERO
   ALLOCATE(ELRKT(0:ntidenode));                 ELRKT         = ZERO
   ALLOCATE(ELTDT(0:ntidenode));                 ELTDT         = ZERO
   ALLOCATE(ELP  (0:ntidenode));                 ELP           = ZERO
   ALLOCATE(ELPF (0:ntidenode));                 ELPF          = ZERO
   ALLOCATE(ELRKP(0:ntidenode));                 ELRKP         = ZERO
   ALLOCATE(UAT  (0:ntidecell));                 UAT           = ZERO
   ALLOCATE(VAT  (0:ntidecell));                 VAT           = ZERO
   ALLOCATE(UATF (0:ntidecell));                 UATF          = ZERO
   ALLOCATE(VATF (0:ntidecell));                 VATF          = ZERO
   ALLOCATE(UAP  (0:ntidecell));                 UAP           = ZERO
   ALLOCATE(VAP  (0:ntidecell));                 VAP           = ZERO
   ALLOCATE(UANT (0:  nmfcell));                 UANT          = ZERO
   ALLOCATE(VANT (0:  nmfcell));                 VANT          = ZERO
   ALLOCATE(UAN  (0:  nmfcell));                 UAN           = ZERO
   ALLOCATE(VAN  (0:  nmfcell));                 VAN           = ZERO
   ALLOCATE(UANP (0:  nmfcell));                 UANP          = ZERO
   ALLOCATE(VANP (0:  nmfcell));                 VANP          = ZERO
   ALLOCATE(UT   (0:ntidecell,1:KBM1));          UT            = ZERO
   ALLOCATE(VT   (0:ntidecell,1:KBM1));          VT            = ZERO
   ALLOCATE(UNT  (0:  nmfcell,1:KBM1));          UNT           = ZERO
   ALLOCATE(VNT  (0:  nmfcell,1:KBM1));          VNT           = ZERO
   ALLOCATE(UN   (0:  nmfcell,1:KBM1));          UN            = ZERO
   ALLOCATE(VN   (0:  nmfcell,1:KBM1));          VN            = ZERO

   ALLOCATE(UAPF (0:ntidecell));                 UAPF          = ZERO
   ALLOCATE(VAPF (0:ntidecell));                 VAPF          = ZERO
   ALLOCATE(UANTF(0:  nmfcell));                 UANTF         = ZERO
   ALLOCATE(VANTF(0:  nmfcell));                 VANTF         = ZERO
   ALLOCATE(UANF (0:  nmfcell));                 UANF          = ZERO
   ALLOCATE(VANF (0:  nmfcell));                 VANF          = ZERO
   ALLOCATE(UANPF(0:  nmfcell));                 UANPF         = ZERO
   ALLOCATE(VANPF(0:  nmfcell));                 VANPF         = ZERO
   ALLOCATE(UARKNT (0:  nmfcell));               UARKNT        = ZERO
   ALLOCATE(VARKNT (0:  nmfcell));               VARKNT        = ZERO
   ALLOCATE(UARKN  (0:  nmfcell));               UARKN         = ZERO
   ALLOCATE(VARKN  (0:  nmfcell));               VARKN         = ZERO
   ALLOCATE(UNTB (0:  nmfcell,1:KBM1));          UNTB          = ZERO
   ALLOCATE(VNTB (0:  nmfcell,1:KBM1));          VNTB          = ZERO
   ALLOCATE(UNB  (0:  nmfcell,1:KBM1));          UNB           = ZERO
   ALLOCATE(VNB  (0:  nmfcell,1:KBM1));          VNB           = ZERO

   RETURN
   END SUBROUTINE ALLOC_OBC2_DATA
!==========================================================================|

!========================================================================
     SUBROUTINE FIND_OBSIDE
     USE ALL_VARS
     IMPLICIT NONE
     INTEGER :: I,ITMP,J,J1,I1,IERR
     INTEGER :: IA,IB
     INTEGER,ALLOCATABLE :: TEMP(:)
     INTEGER,ALLOCATABLE :: NODE_OB(:),CELL_OB(:)
     INTEGER              :: k,NCNT,itemp
     INTEGER, ALLOCATABLE :: temp1(:)
     REAL(SP),ALLOCATABLE :: RTEMP1(:,:),RTEMP2(:,:),RTEMP3(:,:,:),RTEMP4(:,:,:)
     INTEGER              :: i2,i3,ii,JN
     REAL(SP):: DELTX,DELTY,XTMP1,YTMP1,AA1,BB1,CC1,AA2,BB2,CC2
     REAL(SP)::TTIME
     ALLOCATE(NODE_OB(0:MT)); NODE_OB = 0
     ALLOCATE(CELL_OB(0:NT)); CELL_OB = 0
!-----------------------------------Jianzhong-------------------------!
!     DO I1=1,IOBCN
!       J=I_OBC_N(I1)
!       J1=NEXT_OBC(I1)
!       NODE_OB(J) = 1
!       DO I=1,NTVE(J)
!         CELL_OB(NBVE(J,I)) = 1
!       END DO
!     END DO  
     DO I1=1,IBCN(2)
       JN = OBC_LST(2,I1)
       J=I_OBC_N(JN)
       NODE_OB(J) = 1
       DO I=1,NTVE(J)
         CELL_OB(NBVE(J,I)) = 1
       END DO
     END DO 
!---------------------------------------------------------------------! 

     ALLOCATE(TEMP(NE));  TEMP = ZERO
     NOBE = 0
     
     DO I=1,NE
       IA = IEC(I,1)
       IB = IEC(I,2)
       IF(CELL_OB(IA) == 1 .OR. CELL_OB(IB) == 1)THEN
         NOBE = NOBE + 1
	 TEMP(NOBE) = I
       END IF
     END DO
     
     ALLOCATE(COBEDGE_LST(NOBE))
     COBEDGE_LST(1:NOBE) = TEMP(1:NOBE)
     DEALLOCATE(TEMP)
     
     ALLOCATE(TEMP(NCV));  TEMP = ZERO
     NOBCV = 0
     
     DO I=1,NCV
       IA = NIEC(I,1)
       IB = NIEC(I,2)
       IF(NODE_OB(IA) == 1 .OR. NODE_OB(IB) == 1)THEN
         NOBCV = NOBCV + 1
	 TEMP(NOBCV) = I
       END IF
     END DO
     
     ALLOCATE(NOBEDGE_LST(NOBCV))
     NOBEDGE_LST(1:NOBCV) = TEMP(1:NOBCV)
     DEALLOCATE(TEMP)
     
     DEALLOCATE(NODE_OB)
     DEALLOCATE(CELL_OB)

     INMF    =45
     INTCELL =46
     INTNODE =47
     INTELEL =48
     INTUV   =49
     CALL FOPEN(INMF,   TRIM(INPUT_DIR)//TRIM(CASENAME)//'_meanflow.dat'  ,"cfr")
     CALL FOPEN(INTCELL,TRIM(INPUT_DIR)//TRIM(CASENAME)//'_tide_cell.dat' ,"cfr")
     CALL FOPEN(INTNODE,TRIM(INPUT_DIR)//TRIM(CASENAME)//'_tide_node.dat' ,"cfr")
     CALL FOPEN(INTELEL,TRIM(INPUT_DIR)//TRIM(CASENAME)//'_tide_el.dat'   ,"cfr")
     CALL FOPEN(INTUV,  TRIM(INPUT_DIR)//TRIM(CASENAME)//'_tide_uv.dat'   ,"cfr")

!-------------------------------------------------------------------
!
!----Read in Tidal Current Time Series Data----------------------------------
!

     REWIND(INTCELL)
     READ(INTCELL,*) ntidecell_GL

     ntidecell   = 0
     ntidecell_i = 0
  IF (ntidecell_GL > 0) THEN

     ALLOCATE(I_TIDECELL_GL(ntidecell_GL))
     DO I=1,ntidecell_GL
        READ(INTCELL,*)I_TIDECELL_GL(I)
     ENDDO
     CLOSE(INTCELL)

!     IF(ntidecell_GL > 300) THEN
!       WRITE(*,*)'CHANGE FORMAT STATEMENT BELOW TO ACCOMODATE ntidecell_GL='
!       WRITE(*,*)ntidecell_GL,' NUMBER OF TIDAL OPEN BOUNDARY CELLS AND RECOMPILE'
!       CALL PSTOP
!     END IF

     REWIND(INTUV)

!------------------determine the julian forcing counting-------------------------
     IF(MSR)THEN
       CALL FOPEN(111,TRIM(INPUT_DIR)//TRIM(CASENAME)//'_elj_obc.dat',"cfr")
       NCNT = 0
       READ(111,*)
       READ(111,*)
       DO WHILE(.TRUE.)
         READ(111,*,IOSTAT=IOS)
         IF(IOS < 0)EXIT
         NCNT = NCNT + 1
       END DO
       IF(NCNT == 0)CALL FATAL_ERROR("JULIAN TIDE SELECTED BUT NO DATA IN FILE")
     END IF
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(NCNT,1,MPI_INTEGER,0,MPI_FVCOM_GROUP,IERR)
#    endif
     ELO_TM%NTIMES = NCNT
     ALLOCATE(ELO_TM%TIMES(NCNT))
     TTIME = 0.0_SP
     DO I=1,NCNT
       ELO_TM%TIMES(I) = TTIME
       TTIME = TTIME + 720 ! TTIME+DELTT
     END DO
!--------------------------------------------------------------------------------
     ALLOCATE(RTEMP1(ntidecell_GL,ELO_TM%NTIMES))   ! assumig we've known ELO_TM
     ALLOCATE(RTEMP2(ntidecell_GL,ELO_TM%NTIMES))
     ALLOCATE(RTEMP3(ntidecell_GL,KBM1,ELO_TM%NTIMES))
     ALLOCATE(RTEMP4(ntidecell_GL,KBM1,ELO_TM%NTIMES))
     RTEMP1 = 0.0_SP
     RTEMP2 = 0.0_SP
     RTEMP3 = 0.0_SP
     RTEMP4 = 0.0_SP

     IF(MSR)THEN
       DO I=1,ELO_TM%NTIMES
         READ(INTUV,'(I7,300f8.4)') itemp,(RTEMP1(J,I),J=1,ntidecell_GL)
         READ(INTUV,'(I7,300f8.4)') itemp,(RTEMP2(J,I),J=1,ntidecell_GL)
         DO k = 1,KBM1
            READ(INTUV,'(I7,300f8.4)') itemp,(RTEMP3(J,k,I),J=1,ntidecell_GL)
            READ(INTUV,'(I7,300f8.4)') itemp,(RTEMP4(J,k,I),J=1,ntidecell_GL)
         ENDDO
       END DO
     END IF

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(RTEMP1,ntidecell_GL*ELO_TM%NTIMES,MPI_F,0,MPI_FVCOM_GROUP,IERR)
     IF(PAR)CALL MPI_BCAST(RTEMP2,ntidecell_GL*ELO_TM%NTIMES,MPI_F,0,MPI_FVCOM_GROUP,IERR)
     IF(PAR)CALL MPI_BCAST(RTEMP3,ntidecell_GL*KBM1*ELO_TM%NTIMES,MPI_F,0,MPI_FVCOM_GROUP,IERR)
     IF(PAR)CALL MPI_BCAST(RTEMP4,ntidecell_GL*KBM1*ELO_TM%NTIMES,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#    endif
     close(INTUV)

!
!---Map to Local Domain----------------------------------------
!
     IF(SERIAL) THEN
        ntidecell   = ntidecell_GL
        ntidecell_i = ntidecell_GL
        ALLOCATE(I_TIDECELL_N(ntidecell))
        I_TIDECELL_N(:) = I_TIDECELL_GL(:)
        ALLOCATE(UATTS(ntidecell,ELO_TM%NTIMES))
        ALLOCATE(VATTS(ntidecell,ELO_TM%NTIMES))
        ALLOCATE(UTTS (ntidecell,KBM1,ELO_TM%NTIMES))
        ALLOCATE(VTTS (ntidecell,KBM1,ELO_TM%NTIMES))
        UATTS = RTEMP1
        VATTS = RTEMP2
        UTTS  = RTEMP3
        VTTS  = RTEMP4
     ENDIF

#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       ALLOCATE(TEMP1(ntidecell_GL))
       NCNT = 0
       DO I=1,ntidecell_GL
!         I1=ELID_X(I_TIDECELL_GL(I))
         I1=ELID(I_TIDECELL_GL(I))
	 IF(I1 /= 0)THEN
	   NCNT = NCNT + 1
	   TEMP1(NCNT) = I1
	 END IF
       END DO
       ntidecell_i = NCNT

       DO I=1,ntidecell_GL
         I1=ELID_X(I_TIDECELL_GL(I))
         I2=ELID(I_TIDECELL_GL(I))
	 IF(I1 /= 0 .and. I1 /= I2)THEN
	   NCNT = NCNT + 1
	   TEMP1(NCNT) = I1
	 END IF
       END DO
       ntidecell = NCNT
       IF(ntidecell > 0)THEN                          ! if ntidecell = 0
         ALLOCATE(I_TIDECELL_N(ntidecell))
	 I_TIDECELL_N(1:ntidecell) = TEMP1(1:ntidecell)
       END IF
       
       DEALLOCATE(TEMP1)

       if(ntidecell > 0)then
          ALLOCATE(UATTS(ntidecell,ELO_TM%NTIMES));         UATTS = ZERO
          ALLOCATE(VATTS(ntidecell,ELO_TM%NTIMES));         VATTS = ZERO
          ALLOCATE(UTTS (ntidecell,KBM1,ELO_TM%NTIMES));    UTTS  = ZERO
          ALLOCATE(VTTS (ntidecell,KBM1,ELO_TM%NTIMES));    VTTS  = ZERO

          NCNT = 0
          DO I=1,ntidecell_GL
!             I1=ELID_X(I_TIDECELL_GL(I))
             I1=ELID(I_TIDECELL_GL(I))
	     IF(I1 /= 0)THEN
	        NCNT = NCNT + 1
	        UATTS(NCNT,:)   = RTEMP1(I,:)
	        VATTS(NCNT,:)   = RTEMP2(I,:)
	        UTTS (NCNT,:,:) = RTEMP3(I,:,:)
	        VTTS (NCNT,:,:) = RTEMP4(I,:,:)
	     END IF
          END DO

          DO I=1,ntidecell_GL
             I1=ELID_X(I_TIDECELL_GL(I))
             I2=ELID(I_TIDECELL_GL(I))
             IF(I1 /= 0 .and. I1 /= I2)THEN
  	        NCNT = NCNT + 1
	        UATTS(NCNT,:)   = RTEMP1(I,:)
	        VATTS(NCNT,:)   = RTEMP2(I,:)
	        UTTS (NCNT,:,:) = RTEMP3(I,:,:)
	        VTTS (NCNT,:,:) = RTEMP4(I,:,:)
	     END IF
          END DO
       endif

     END IF
#   endif

     DEALLOCATE(RTEMP1,RTEMP2,RTEMP3,RTEMP4)

#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
#  endif
  
  ELSE ! if statement end for ntidecell_GL > 0
    close(INTCELL)
  END IF
!-------------------------------------------------------------------
!
!----Read in Tidal Elevation Time Series Data----------------------------------
!

     REWIND(INTNODE)
     READ(INTNODE,*) ntidenode_GL

     ntidenode   = 0
     ntidenode_i = 0
  IF (ntidenode_GL > 0) THEN

     ALLOCATE(I_TIDENODE_GL(ntidenode_GL))
     DO I=1,ntidenode_GL
        READ(INTNODE,*)I_TIDENODE_GL(I)
     ENDDO
     CLOSE(INTNODE)

!     IF(ntidenode_GL > 300 ) THEN
!       WRITE(*,*)'CHANGE FORMAT STATEMENT BELOW TO ACCOMODATE ntidenode_GL=',ntidenode_GL
!       WRITE(*,*)' NUMBER OF TIDAL OPEN BOUNDARY NODES AND RECOMPILE'
!       CALL PSTOP
!     END IF

     ALLOCATE(RTEMP1(ntidenode_GL, ELO_TM%NTIMES))
     RTEMP1 = 0.0_SP

     REWIND(INTELEL)
     IF(MSR)THEN
       DO I=1,ELO_TM%NTIMES
         READ(INTELEL,'(I7,300f8.4)') itemp,(RTEMP1(J,I),J=1,ntidenode_GL)
       END DO
     END IF
     CLOSE(INTELEL)

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL MPI_BCAST(RTEMP1,ntidenode_GL*ELO_TM%NTIMES,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#    endif

!
!---Map to Local Domain----------------------------------------
!
     IF(SERIAL) THEN
        ntidenode    = ntidenode_GL
        ntidenode_i  = ntidenode_GL
        ALLOCATE(I_TIDENODE_N(ntidenode))
        I_TIDENODE_N = I_TIDENODE_GL
        ALLOCATE(ELTTS(ntidenode,ELO_TM%NTIMES)); ELTTS = ZERO
        ELTTS = RTEMP1
     ENDIF

#    if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       ALLOCATE(TEMP1(ntidenode_GL))
       NCNT = 0
       DO I=1,ntidenode_GL
!         I1=NLID_X(I_TIDENODE_GL(I))
         I1=NLID(I_TIDENODE_GL(I))
	 IF(I1 /= 0)THEN
	   NCNT = NCNT + 1
	   TEMP1(NCNT) = I1
	 END IF
       END DO
       ntidenode_i = NCNT

       DO I=1,ntidenode_GL
         I1=NLID_X(I_TIDENODE_GL(I))
         I2=NLID(I_TIDENODE_GL(I))
	 IF(I1 /= 0 .and. I1 /= I2)THEN
	   NCNT = NCNT + 1
	   TEMP1(NCNT) = I1
	 END IF
       END DO
       ntidenode = NCNT
       IF(ntidenode > 0)THEN                          ! if ntidenode = 0
         ALLOCATE(I_TIDENODE_N(ntidenode))
	 I_TIDENODE_N(1:ntidenode) = TEMP1(1:ntidenode)
       END IF
       
       DEALLOCATE(TEMP1)

       if(ntidenode > 0)then
          ALLOCATE(ELTTS(ntidenode,ELO_TM%NTIMES)); ELTTS = ZERO

          NCNT = 0
          DO I=1,ntidenode_GL
!             I1=NLID_X(I_TIDENODE_GL(I))
             I1=NLID(I_TIDENODE_GL(I))
	     IF(I1 /= 0)THEN
	        NCNT = NCNT + 1
	        ELTTS(NCNT,:)   = RTEMP1(I,:)
	     END IF
          END DO

          DO I=1,ntidenode_GL
             I1=NLID_X(I_TIDENODE_GL(I))
             I2=NLID(I_TIDENODE_GL(I))
	     IF(I1 /= 0 .and. I1 /= I2)THEN
	        NCNT = NCNT + 1
	        ELTTS(NCNT,:)   = RTEMP1(I,:)
	     END IF
          END DO
       endif

     END IF
#   endif

     DEALLOCATE(RTEMP1)

#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
#  endif

  ELSE    ! if statement end for ntidenode_GL > 0
    close(INTNODE)
  END IF

   CALL READ_MEANFLOW
   CALL SET_BNDRY_MEANFLOW

!
!--- calculate mapping function
!
   ALLOCATE(I_OBC_CELL(0:NT),I_OBC_NODE(0:MT),I_OBC_CELL2(0:NT),I_OBC_NODE2(0:MT))

      I_OBC_NODE = 0   
! I_OBC_NODE is the mapping from local domain node index j(MT) to tidal open 
! bndy node index i. If not a tidal open bndy node, I_OBC_NODE(j)=0
    if(ntidenode > 0)then
        do j = 1, MT
        do i = 1, ntidenode
           I1 = I_TIDENODE_N(i)
           if (I1 == j) then
              I_OBC_NODE(J) = i
           endif
        enddo
        enddo
    endif

      I_OBC_NODE2 = 0   
! I_OBC_NODE2 is the mapping from local domain node index j(MT) to open 
! bndy node index i. If not a open bndy node, I_OBC_NODE2(j)=0
!-------------------------------Jianzhong----------------------------!
!    if(iobcn > 0)then
!        do j = 1, MT
!        do i = 1, iobcn
!           I1 = I_OBC_N(i)
!           if (I1 == j) then
!              I_OBC_NODE2(J) = i
!           endif
!        enddo
!        enddo
!    endif
    if(ibcn(2) > 0)then
        do j = 1, MT
        do i = 1, ibcn(2)
           jn=obc_lst(2,i)
           I1 = I_OBC_N(jn)
           if (I1 == j) then
              I_OBC_NODE2(J) = i
           endif
        enddo
        enddo
    endif
!--------------------------------------------------------------------!

     I_OBC_CELL = 0   
! I_OBC_CELL is the mapping from local domain cell index j(NT) to tidal open 
! bndy cell index i. if not a tidal open bndy cell, I_OBC_CELL(j)=0
     if(ntidecell > 0)then
        do j = 1, NT
        do i = 1, ntidecell            
           I1 = I_TIDECELL_N(i)        
           if (I1 == j) then
              I_OBC_CELL(J) = i
           endif
        enddo
        enddo
     endif

     I_OBC_CELL2 = 0   
! I_OBC_CELL2 is the mapping from local domain cell index j(NT) to mean flow open 
! bndy cell index i. if not a mean flow open bndy cell, I_OBC_CELL2(j)=0
     if(nmfcell > 0)then
        do j = 1, NT
        do i = 1, nmfcell            
           I1 = I_MFCELL_N(i)        
           if (I1 == j) then
              I_OBC_CELL2(J) = i
           endif
        enddo
        enddo
     endif


  IF(nmfcell > 0)THEN
     ALLOCATE (DLTN(nmfcell))
  DO I = 1,nmfcell
     II = I_MFCELL_N(I) 
     IF(NBE(II,1) == 0 .and. ISONB(nv(II,1)) /= 2) THEN
       DELTX = VX(NV(II,2))-VX(NV(II,3))
# if defined (SPHERICAL)
       IF(DELTX > 180.0_SP)THEN
         DELTX = -360.0_SP+DELTX
       ELSE IF(DELTX < -180.0_SP)THEN
         DELTX =  360.0_SP+DELTX	 
       END IF	 
# endif       
       DELTY = VY(NV(II,2))-VY(NV(II,3))
       AA1 = -DELTY
       BB1 = DELTX
       CC1 = -AA1*VX(NV(II,2))-BB1*VY(NV(II,2))
     ELSE IF(NBE(II,2) == 0 .and. ISONB(nv(II,2)) /= 2) THEN
       DELTX = VX(NV(II,3))-VX(NV(II,1))
# if defined (SPHERICAL)
       IF(DELTX > 180.0_SP)THEN
         DELTX = -360.0_SP+DELTX
       ELSE IF(DELTX < -180.0_SP)THEN
         DELTX =  360.0_SP+DELTX	 
       END IF	 
# endif       
       DELTY = VY(NV(II,3))-VY(NV(II,1))
       AA1 = -DELTY
       BB1 = DELTX
       CC1 = -AA1*VX(NV(II,3))-BB1*VY(NV(II,3))
     ELSE IF(NBE(II,3) == 0 .and. ISONB(nv(II,3)) /= 2) THEN
       DELTX = VX(NV(II,1))-VX(NV(II,2))
# if defined (SPHERICAL)
       IF(DELTX > 180.0_SP)THEN
         DELTX = -360.0_SP+DELTX
       ELSE IF(DELTX < -180.0_SP)THEN
         DELTX =  360.0_SP+DELTX	 
       END IF	 
# endif
       DELTY = VY(NV(II,1))-VY(NV(II,2))
       AA1 = -DELTY
       BB1 = DELTX
       CC1 = -AA1*VX(NV(II,1))-BB1*VY(NV(II,1))
     ELSE
       PRINT*,'something is wrong here 1'
       CALL PSTOP
     END IF

     AA2 = BB1
     BB2 = -AA1
     CC2 = -AA2*XC(II)-BB2*YC(II)

     XTMP1 = -(CC1*BB2-CC2*BB1)/(AA1*BB2-AA2*BB1)
     YTMP1 = -(CC1*AA2-CC2*AA1)/(BB1*AA2-BB2*AA1)
     
     DLTN(I) = SQRT((XC(II)-XTMP1)**2+(YC(II)-YTMP1)**2)
  END DO
  END IF

     RETURN
     END SUBROUTINE FIND_OBSIDE


!==============================================================================|
!  INTERPOLATION of VARIOUS TIME SERIES                                        |
!==============================================================================|

   SUBROUTINE BCOND_TIDE_2D
   USE CONTROL
   USE BCS
   USE MOD_OBCS
   INTEGER  L1,L2,IERR
   REAL(SP) :: FACT,UFACT,TIME1
   REAL(SP) :: TIMERK
   
   TIMERK = SECONDS(RKTIME-STARTTIME)
!   TIME1 = TIMERK * 86400.0_SP
   TIME1 = TIMERK

   CALL BRACKET(ELO_TM,TIME1,L1,L2,FACT,UFACT,IERR)
   IF(ntidecell > 0)THEN
     UATF(1:ntidecell) = UFACT*UATTS(1:ntidecell,L1) + FACT*UATTS(1:ntidecell,L2)
     VATF(1:ntidecell) = UFACT*VATTS(1:ntidecell,L1) + FACT*VATTS(1:ntidecell,L2)
     UATF(1:ntidecell) = UATF(1:ntidecell) * RAMP
     VATF(1:ntidecell) = VATF(1:ntidecell) * RAMP
   END IF

   IF(ntidenode > 0)THEN
     ELTF(1:ntidenode) = UFACT*ELTTS(1:ntidenode,L1) + FACT*ELTTS(1:ntidenode,L2)
     ELTF(1:ntidenode) = ELTF(1:ntidenode) * RAMP
   END IF

   RETURN
   END SUBROUTINE BCOND_TIDE_2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE BCOND_TIDE_3D
   USE BCS
   USE MOD_OBCS
   INTEGER  L1,L2,IERR
   REAL(SP) :: FACT,UFACT,TIME1

!   TIME1 = TIME * 86400.0_SP - dti     ! pay attention to this TIME (different between 2D and 3D)
   TIME1=SECONDS(intTime-StartTime)-DTI
   IF(ntidecell > 0)THEN
     CALL BRACKET(ELO_TM,TIME1,L1,L2,FACT,UFACT,IERR)
     UT(1:ntidecell,:) = UFACT*UTTS(1:ntidecell,:,L1) + FACT*UTTS(1:ntidecell,:,L2)
     VT(1:ntidecell,:) = UFACT*VTTS(1:ntidecell,:,L1) + FACT*VTTS(1:ntidecell,:,L2)
     UT(1:ntidecell,:) = UT(1:ntidecell,:) * RAMP
     VT(1:ntidecell,:) = VT(1:ntidecell,:) * RAMP
   END IF

   RETURN
   END SUBROUTINE BCOND_TIDE_3D


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==========================================================================|
   SUBROUTINE BCOND_NG_2D

!--------------------------------------------------------------------------|
!  NON-Gradient Open Boundary Condition (2-D)                              |
!--------------------------------------------------------------------------|

   USE ALL_VARS
   IMPLICIT NONE

   INTEGER :: I,J,J1

   IF(ntidecell > 0)THEN
   DO I = 1, ntidecell
      J = I_TIDECELL_N(I)
          UAP(I) = UA(J) - UAT(I)
          VAP(I) = VA(J) - VAT(I)
   END DO
   END IF

   IF(nmfcell > 0)THEN
   DO I = 1, nmfcell
      J = I_MFCELL_N(I)
      J1= I_OBC_CELL(J)
      UANT(I) = UAT(J1)
      VANT(I) = VAT(J1)
      UAN (I) = UA (J )
      VAN (I) = VA (J )
      UANP(I) = UAP(J1)
      VANP(I) = VAP(J1)
   END DO
   END IF

   RETURN
   END SUBROUTINE BCOND_NG_2D
!==========================================================================|

!==========================================================================|
   SUBROUTINE BCOND_NG_3D

!--------------------------------------------------------------------------|
!  NON-Gradient Open Boundary Condition (3-D)                              |
!--------------------------------------------------------------------------|

   USE ALL_VARS
   IMPLICIT NONE

   INTEGER :: I,J,J1,K

   IF(nmfcell > 0)THEN
   DO I = 1, nmfcell
      J = I_MFCELL_N(I)
      J1= I_OBC_CELL(J)
      DO K = 1,KBM1
         UNT(I,K) = UT(J1,K)
         VNT(I,K) = VT(J1,K)
         UN (I,K) = U (J ,K)
         VN (I,K) = V (J ,K)
     END DO
   END DO
   END IF

   RETURN
   END SUBROUTINE BCOND_NG_3D
!==========================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==========================================================================|
   SUBROUTINE BCOND_BKI_2D(KTT)

!--------------------------------------------------------------------------|
!  BLUMBERG AND KANTHA IMPLICIT OPEN BOUNDARY CONDITION (BKI)              |
!--------------------------------------------------------------------------|

   USE ALL_VARS
   IMPLICIT NONE

   INTEGER, INTENT(IN)  :: KTT
   INTEGER :: I,II,J1,J
   REAL(SP):: CC,CP
   REAL(SP):: coef


   coef = 10800.00_SP

   IF(ntidecell > 0)THEN
   DO I = 1, ntidecell
      J = I_TIDECELL_N(I)
          UAPF(I) = UAF(J) - UATF(I)
          VAPF(I) = VAF(J) - VATF(I)
          UAP(I)  = UAPF(I)
          VAP(I)  = VAPF(I)
   END DO
   END IF

  IF(nmfcell > 0)THEN
  DO I = 1,nmfcell
     II = I_MFCELL_N(I) 
     J1= I_OBC_CELL(II)

     CC = SQRT(GRAV_E(II)*D1(II))*ALPHA_RK(KTT)*DTE/DLTN(I)
     CP = CC + 1.0_SP
       UANTF(I) = (CC*UATF(J1) + UARKNT(I)*(1.0_SP-ALPHA_RK(KTT)*DTE/coef))/CP
       VANTF(I) = (CC*VATF(J1) + VARKNT(I)*(1.0_SP-ALPHA_RK(KTT)*DTE/coef))/CP
!       UANF (I) = (CC*UAF (II) + UARKN(I) *(1.0_SP-ALPHA_RK(KTT)*DTE/coef))/CP
!       VANF (I) = (CC*VAF (II) + VARKN(I) *(1.0_SP-ALPHA_RK(KTT)*DTE/coef))/CP
       UANF (I) = 0.0_SP
       VANF (I) = 0.0_SP

       UANPF(I) = UANF(I) - UANTF(I)
       VANPF(I) = VANF(I) - VANTF(I)
       UANT (I) = UANTF(I)
       VANT (I) = VANTF(I)
       UAN  (I) = UANF(I)
       VAN  (I) = VANF(I)
       UANP (I) = UANPF(I)
       VANP (I) = VANPF(I)
  END DO
  END IF

   RETURN
   END SUBROUTINE BCOND_BKI_2D
!==========================================================================|

!==========================================================================|
   SUBROUTINE BCOND_BKI_3D(KTT)

!--------------------------------------------------------------------------|
!  BLUMBERG AND KANTHA IMPLICIT OPEN BOUNDARY CONDITION (BKI)              |
!--------------------------------------------------------------------------|

   USE ALL_VARS
   IMPLICIT NONE

   INTEGER, INTENT(IN)  :: KTT
   INTEGER :: I,II,K,J1
   REAL(SP):: CC,CP
   REAL(SP):: coef

   coef = 10800.00_SP

  IF(nmfcell > 0)THEN
  DO I = 1,nmfcell
     II = I_MFCELL_N(I) 
     J1= I_OBC_CELL(II)

     CC = SQRT(GRAV_E(II)*D1(II))*DTI/DLTN(I)
     CP = CC + 1.0_SP
     DO K=1,KBM1
       UNT(I,K)  = (CC*UT(J1,K) + UNTB(I,K)*(1.0_SP-DTI/coef))/CP
       VNT(I,K)  = (CC*VT(J1,K) + VNTB(I,K)*(1.0_SP-DTI/coef))/CP
!       UN (I,K)  = (CC*U (II,K) + UNB (I,K)*(1.0_SP-DTI/coef))/CP
!       VN (I,K)  = (CC*V (II,K) + VNB (I,K)*(1.0_SP-DTI/coef))/CP
       UN (I,K)  = 0.0_SP
       VN (I,K)  = 0.0_SP
       IF (KTT == 2) THEN
          UNTB(I,K) = UNT(I,K)
          VNTB(I,K) = VNT(I,K)
          UNB (I,K) = UN (I,K)
          VNB (I,K) = VN (I,K)
       END IF
     END DO
  END DO
  END IF

   RETURN
   END SUBROUTINE BCOND_BKI_3D
!==========================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!--------------------------------------------------------------------------|
!               TEST WHETHER THE TWO BOUNDARY CONDITION WORK               |
!--------------------------------------------------------------------------|
   SUBROUTINE TEST_CELL(INDEX,MSG)
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif 

   IMPLICIT NONE

   CHARACTER(LEN=*)             :: MSG
   INTEGER                      :: INDEX
   INTEGER,PARAMETER            :: CELL1=59     ! IN TAIWAN STRAIT
   INTEGER,PARAMETER            :: CELL2=229    ! IN SECOND OPEN BOUNDARY



!   IF(CELL1==EGID(INDEX))THEN
!     WRITE(IPT_P,*)'IN '//MSG//' CELL OF TAIWAN STARIT  IN',MYID
!     write(ipt_p,*) "INDEX:",index,"egid(index)",egid(index),"; msg="//TRIM(MSG)//"; CELLS:",CELL1,CELL2
!   END IF
!   IF(CELL2==EGID(INDEX))THEN
!     WRITE(IPT_P,*)'IN '//MSG//' CELL OF PACIFIC BOUNDARY IN',MYID
!     write(ipt_p,*) "INDEX:",index,"egid(index)",egid(index),"; msg="//TRIM(MSG)//"; CELLS:",CELL1,CELL2
!   END IF

   RETURN
   END SUBROUTINE TEST_CELL
!==========================================================================|

   SUBROUTINE TEST_NODE(INDEX,MSG)
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif 
   IMPLICIT NONE

   CHARACTER(LEN=*)             :: MSG
   INTEGER                      :: INDEX
   INTEGER,PARAMETER            :: NODE1=30     ! IN TAIWAN STRAIT
   INTEGER,PARAMETER            :: NODE2=115    ! IN SECOND OPEN BOUNDARY


!   IF(NODE1==NGID(INDEX))THEN
!     WRITE(IPT_P,*)'IN '//MSG//' NODE OF TAIWAN STARIT IN',MYID
!     write(ipt_p,*) "INDEX:",index,"ngid(index)",NGID(index),"; msg="//TRIM(MSG)//"; NODES:",NODE1,NODE2
!   END IF
!   IF(NODE2==NGID(INDEX))THEN
!     WRITE(IPT_P,*)'IN '//MSG//' NODE OF PACIFIC BOUNDARY IN',MYID
!     write(ipt_p,*) "INDEX:",index,"ngid(index)",NGID(index),"; msg="//TRIM(MSG)//"; NODES:",NODE1,NODE2
!   END IF

   RETURN
   END SUBROUTINE TEST_NODE
!==========================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

# endif
END MODULE MOD_OBCS2
