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

MODULE MOD_RRKA 
# if defined (RRKF)
   USE RRKVAL
   USE CONTROL
   IMPLICIT NONE
   SAVE
      real(dp),PARAMETER :: one_db=dble(1.0),zero_db=dble(0.0),two_db=dble(2.0)
   INTEGER :: STDIM   
!   INTEGER NLOC
   INTEGER IDUM
   INTEGER NCYC
   INTEGER N1CYC
!   INTEGER ICYC
   INTEGER I_INITIAL
   INTEGER INORRK
   INTEGER OUTERR
   INTEGER,ALLOCATABLE :: STLOC(:)
   REAL(DP),ALLOCATABLE :: WKTMP(:)     !TEMP ARRAY FOR WK()
   REAL(DP),ALLOCATABLE :: STFCT(:)
   REAL(DP),ALLOCATABLE :: STANL(:)
   REAL(DP),ALLOCATABLE :: CRSTATE(:)
   REAL(DP),ALLOCATABLE :: KAL(:,:)
   REAL(DP),ALLOCATABLE :: MODDATA(:)
   REAL(DP),ALLOCATABLE :: OBSDATA(:)
  REAL(DP),ALLOCATABLE ::  INNOV(:)    !pengfei
   REAL(DP),ALLOCATABLE :: ERRVEC(:)
   REAL(DP),ALLOCATABLE :: R(:,:)
   REAL(DP),ALLOCATABLE :: OBSERR(:)
   REAL(DP),ALLOCATABLE :: STTRUE(:,:)
   REAL(DP),ALLOCATABLE :: STTRUE1(:)
   REAL(DP),ALLOCATABLE :: STTR1(:)
   REAL(DP),ALLOCATABLE :: STTR0(:)
!   REAL(DP),ALLOCATABLE :: SEOF(:,:)
!   REAL(DP),ALLOCATABLE :: STTEMP1(:)
!   REAL(DP),ALLOCATABLE :: STTEMP0(:)

!   REAL(DP),ALLOCATABLE :: STSD(:)
   REAL(SP),ALLOCATABLE :: RRKEL2(:)
   REAL(SP),ALLOCATABLE :: RRKU2(:,:)
   REAL(SP),ALLOCATABLE :: RRKV2(:,:)
   REAL(SP),ALLOCATABLE :: RRKT2(:,:)
   REAL(SP),ALLOCATABLE :: RRKS2(:,:)
!   REAL(DP)             :: ELSD, USD, TSD, SSD
   REAL(DP)             :: RRKSUM
   REAL(DP)             :: ERR2_INN_FCT
   REAL(DP)             :: ERR2_INN_AN1
   REAL(DP)             :: ERR2_TOT_FCT
   REAL(DP)             :: ERR2_TOT_AN1   
   REAL                 :: RNOBS
     REAL(SP), ALLOCATABLE :: UETMP(:,:)
   REAL(SP), ALLOCATABLE :: VETMP(:,:)
   REAL(SP), ALLOCATABLE :: SETMP(:,:)
   REAL(SP), ALLOCATABLE :: TETMP(:,:)
   REAL(SP), ALLOCATABLE :: ELETMP(:)
 
   CHARACTER(LEN=6)     :: FCYC
   CHARACTER(LEN=20)    :: TEXT
   CHARACTER(LEN=80)    :: ERRFILE
   CHARACTER(LEN=80)    :: AN1FILE
   CHARACTER(LEN=80)    :: FILENAME   

   INTEGER  TIMEN
   REAL(DP),ALLOCATABLE  :: EL_SRS(:,:),SRS_TMP(:)
   REAL(DP),ALLOCATABLE  :: TIME_SER(:)
   REAL(DP) AMP(6),PHAS(6),PHAI_IJ,FORCE

   CONTAINS

   SUBROUTINE RRK_RRKF
!-------------------------------------------------------------------------------|
! START OF UPDATING THE FORCAST BY THE RRKF                                     | 
!-------------------------------------------------------------------------------|

   USE LIMS
   USE CONTROL
   USE ALL_VARS
   USE RRKVAL
   USE MOD_RRK
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   
   INTEGER IDUMMY
   INTEGER I,J,K
   REAL(DP) :: ERR1
   REAL(DP) :: ERR2
   REAL(DP) :: ERR3
   
   INTEGER SS_DIM
!   INTEGER NLOC
   INTEGER IERR
   
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: UATMP,VATMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TTMP,STMP
   
   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + MGL
   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1


   CALL ALLOC_RRKA
   CALL OBS2VEC
!   OPEN(656,FILE='Y.dat')
!   write(656,'(f12.6)') obsdata
!   close(656)
   SS_DIM = RRK_NEOF
   
   IF (SERIAL) THEN
     IF(EL_ASSIM) THEN
       DO I=1, MGL
         RRKEL2(I) = EL(I)
       ENDDO
     ENDIF
     IF(UV_ASSIM) THEN  
       DO K=1, KBM1
         DO I=1, NGL
           RRKU2(I,K)  = U(I,K)
           RRKV2(I,K)  = V(I,K)
         ENDDO
       ENDDO 
     ENDIF
     IF(T_ASSIM) THEN
       DO K=1, KBM1
         DO I=1, MGL
           RRKT2(I,K)  = T1(I,K)
         ENDDO
       ENDDO
     ENDIF
     IF(S_ASSIM) THEN
       DO K=1, KBM1
         DO I=1, MGL
           RRKS2(I,K)  = S1(I,K)
         ENDDO
       ENDDO
     ENDIF
   ELSE
#    if defined(MULTIPROCESSOR)  
if (par) then
IF(EL_ASSIM) THEN
CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,EL, RRKEL2)
CALL MPI_BCAST(RRKEL2,MGL,MPI_F,0,MPI_FVCOM_GROUP,IERR)
END IF
IF(UV_ASSIM) THEN
CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,U, RRKU2)
CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,V, RRKV2)
CALL MPI_BCAST(RRKU2,NGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
CALL MPI_BCAST(RRKV2,NGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
END IF
IF(T_ASSIM) THEN
CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,T1, RRKT2)
CALL MPI_BCAST(RRKT2,MGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
END IF
IF(S_ASSIM) THEN
CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,S1, RRKS2)
CALL MPI_BCAST(RRKS2,MGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
END IF
end if



!       IF(PAR)THEN   
!         ALLOCATE(ELTMP(MGL))        ; ELTMP = 0.0_SP
!         ALLOCATE(UATMP(NGL))        ; UATMP = 0.0_SP
!         ALLOCATE(VATMP(NGL))        ; VATMP = 0.0_SP
!         ALLOCATE(UTMP(NGL,KB))      ; UTMP  = 0.0_SP
!         ALLOCATE(VTMP(NGL,KB))      ; VTMP  = 0.0_SP
!	 ALLOCATE(TTMP(MGL,KB))      ; TTMP  = 0.0_SP
!	 ALLOCATE(STMP(MGL,KB))      ; STMP  = 0.0_SP

!         CALL GATHER(LBOUND(EL,1), UBOUND(EL,1), M,MGL,1 ,MYID,NPROCS,NMAP,EL, ELTMP)
!         CALL GATHER(LBOUND(U,1),  UBOUND(U,1), N,NGL,KB,MYID,NPROCS,EMAP, U,  UTMP) 
!         CALL GATHER(LBOUND(V,1),  UBOUND(V,1), N,NGL,KB,MYID,NPROCS,EMAP, V,  VTMP) 
!         CALL GATHER(LBOUND(UA,1), UBOUND(UA,1), N,NGL,1 ,MYID,NPROCS,EMAP,UA, UATMP) 
!         CALL GATHER(LBOUND(VA,1), UBOUND(VA,1), N,NGL,1 ,MYID,NPROCS,EMAP,VA, VATMP)
!	 CALL GATHER(LBOUND(T1,1), UBOUND(T1,1), M,MGL,KB,MYID,NPROCS,EMAP,T1,  TTMP)
!	 CALL GATHER(LBOUND(S1,1), UBOUND(S1,1), M,MGL,KB,MYID,NPROCS,EMAP,S1,  STMP)
	 
!	 IF(EL_ASSIM) THEN
!	   DO I=1, MGL
!             RRKEL2(I) = ELTMP(I)
!           ENDDO
!           CALL MPI_BCAST(RRKEL2,MGL,MPI_F,0,MPI_FVCOM_GROUP,IERR)
!	 ENDIF
!         IF(UV_ASSIM) THEN
!	   DO K=1, KBM1
!             DO I=1, NGL
!               RRKU2(I,K)  = UTMP(I,K)
!               RRKV2(I,K)  = VTMP(I,K)
!             ENDDO
!           ENDDO
!           CALL MPI_BCAST(RRKU2,NGL*KB,MPI_F,0,MPI_FVCOM_GROUP,IERR)
!           CALL MPI_BCAST(RRKV2,NGL*KB,MPI_F,0,MPI_FVCOM_GROUP,IERR)
!	 ENDIF
!	 IF(T_ASSIM) THEN
!	   DO K=1, KBM1
!             DO I=1, MGL
!               RRKT2(I,K)  = TTMP(I,K)
!             ENDDO
!           ENDDO
!           CALL MPI_BCAST(RRKT2,MGL*KB,MPI_F,0,MPI_FVCOM_GROUP,IERR)
!	 ENDIF
!	 IF(S_ASSIM) THEN
!	   DO K=1, KBM1
!             DO I=1, MGL
!               RRKS2(I,K)  = STMP(I,K)
!             ENDDO
!           ENDDO
!	   CALL MPI_BCAST(RRKS2,MGL*KB,MPI_F,0,MPI_FVCOM_GROUP,IERR)
!	 ENDIF	 
       
!         DEALLOCATE(ELTMP,UTMP,VTMP,TTMP,STMP,UATMP,VATMP)
!       END IF
#    endif
   ENDIF
print *, 'save the forecast state'
!rrkf   Save the forecast state as 'stfct' 
    IDUMMY = 0
    IF(EL_ASSIM) THEN
      DO I=1, MGL
        IDUMMY = IDUMMY + 1
        STFCT(IDUMMY)= RRKEL2(I)
      ENDDO
    ENDIF
    IF(UV_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NGL 
           IDUMMY = IDUMMY + 1         
	   STFCT(IDUMMY) = RRKU2(J,I)
        ENDDO  
      ENDDO
      DO I=1, KBM1
        DO J=1, NGL
           IDUMMY = IDUMMY + 1
           STFCT(IDUMMY) = RRKV2(J,I)
        ENDDO
      ENDDO
    ENDIF
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL 
           IDUMMY = IDUMMY + 1         
	   STFCT(IDUMMY) = RRKT2(J,I)
        ENDDO  
      ENDDO
    ENDIF
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL 
           IDUMMY = IDUMMY + 1         
	   STFCT(IDUMMY) = RRKS2(J,I)
        ENDDO  
      ENDDO
    ENDIF
print *, 'get true ocean state'
!   Get true ocean state

    CALL READRESTART(NC_START,RRK_CYC)
    
    print *,'RRK_CYC', RRK_CYC

RRKEL2=0.0
RRKU2=0.0
RRKV2=0.0
RRKT2=0.0
RRKS2=0.0
   IF (SERIAL) THEN
     IF(EL_ASSIM) THEN
       DO I=1, MGL
         RRKEL2(I) = EL(I)
       ENDDO
     ENDIF
     IF(UV_ASSIM) THEN
       DO K=1, KBM1
         DO I=1, NGL
           RRKU2(I,K)  = U(I,K)
           RRKV2(I,K)  = V(I,K)
         ENDDO
       ENDDO
     ENDIF
     IF(T_ASSIM) THEN
       DO K=1, KBM1
         DO I=1, MGL
           RRKT2(I,K)  = T1(I,K)
         ENDDO
       ENDDO
     ENDIF
     IF(S_ASSIM) THEN
       DO K=1, KBM1
         DO I=1, MGL
           RRKS2(I,K)  = S1(I,K)
         ENDDO
       ENDDO
     ENDIF
   ELSE
#    if defined(MULTIPROCESSOR)
if (par) then
IF(EL_ASSIM) THEN
CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,EL, RRKEL2)
CALL MPI_BCAST(RRKEL2,MGL,MPI_F,0,MPI_FVCOM_GROUP,IERR)
END IF
IF(UV_ASSIM) THEN
CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,U, RRKU2)
CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,V, RRKV2)
CALL MPI_BCAST(RRKU2,NGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
CALL MPI_BCAST(RRKV2,NGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
END IF
IF(T_ASSIM) THEN
CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,T1, RRKT2)
CALL MPI_BCAST(RRKT2,MGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
END IF
IF(S_ASSIM) THEN
CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,S1, RRKS2)
CALL MPI_BCAST(RRKS2,MGL*KBM1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
END IF
end if
# endif
END IF


IDUMMY=0    
    IF(EL_ASSIM) THEN
      DO I=1, MGL
        IDUMMY = IDUMMY + 1
          STTRUE1(IDUMMY) = RRKEL2(I)
      ENDDO
    ENDIF

    IF(UV_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, NGL
          IDUMMY = IDUMMY + 1
	        STTRUE1(IDUMMY) = RRKU2(J,I)   
        ENDDO
      ENDDO
      DO I=1, KBM1
        DO J=1, NGL
          IDUMMY = IDUMMY + 1
	        STTRUE1(IDUMMY) = RRKV2(J,I)
        ENDDO
      ENDDO		
    ENDIF  
  
    IF(T_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL
          IDUMMY = IDUMMY + 1
    	    STTRUE1(IDUMMY) = RRKT2(J,I)
        ENDDO
      ENDDO
    ENDIF
  
    IF(S_ASSIM) THEN
      DO I=1, KBM1
        DO J=1, MGL
          IDUMMY = IDUMMY + 1
    	    STTRUE1(IDUMMY) = RRKS2(J,I)
        ENDDO
      ENDDO
    ENDIF
    
 print *, 'readobs', obsdata(1)
!------------------!only for idealized case to get obsdata    
!    CALL READOBS(STLOC,NLOC) 
       
!     DO I=1, NLOC
!      OBSDATA(I) = STTRUE1(STLOC(I))
!     ENDDO
!----------------------------------        

!   Perturb the observation with random errors of 
!   prescribed observatin error 'obserr'
!   Generate the random number of Gaussian deviation 
!   with mean of 0 and one standard deviation of 1
   
    DO I=1, NLOC
      CALL GASDEV(IDUM,RNOBS)
      OBSDATA(I) = OBSDATA(I) + wktmp(i)*RNOBS
    ENDDO

!   H(x^f)   -> moddata
!   Get the forecast values at the observation locations
    CALL H_MAPPING(STFCT,MODDATA)
    DO J=1,NLOC
!      MODDATA(J)=STFCT(STLOC(J))
!      write(myid+600,*) MODDATA(J), OBSDATA(j)
    ENDDO
    
!   Calculate innovation vector y' = y - H(x^f) ->moddata
    DO I = 1, NLOC
      INNOV(I) = OBSDATA(I) - MODDATA(I)
!      write(myid+800,*) obsdata(i),moddata(i),INNOV(I)
    ENDDO

!  Calculate the rms error at observational sites   
   TEXT='FCT(obs)'
   ERR1 = 0.0_DP
   ERR2 = 0.0_DP
   ERR3 = 0.0_DP
   DO K = 1, NLOC
     IF(ABS(INNOV(K)) > ERR1) THEN
        ERR1 = ABS(INNOV(K))
     ENDIF
     ERR2 = ERR2 + INNOV(K) * INNOV(K)
     ERR3 = ERR3 + ABS(INNOV(K))
   ENDDO
     ERR2 = DSQRT(ERR2/NLOC)
     ERR3 = ERR3 / NLOC
     ERR2_INN_FCT = ERR2   


!rrkf  Calculate the rms error of entire domain
   DO I=1,STDIM
     ERRVEC(I)=STFCT(I)-STTRUE1(I)
!     write(myid+700,*) STFCT(I), STTRUE1(I)
   ENDDO
   text='FCT(full)'
   ERR1 = 0.
   ERR2 = 0.
   ERR3 = 0.
   DO K = 1, STDIM
     IF(ABS(ERRVEC(K)) > ERR1) THEN
        ERR1 = ABS(ERRVEC(K))
     ENDIF
     ERR2 = ERR2 + ERRVEC(K) * ERRVEC(K)
     ERR3 = ERR3 + ABS(ERRVEC(K))
   ENDDO
     ERR2 = DSQRT(ERR2/STDIM)
     ERR3 = ERR3 / STDIM
     ERR2_TOT_FCT = ERR2

!  Read the Kalman gain matrix
   FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/rrK.dat'
   OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
   DO J=1, NLOC
     READ(INORRK) (KAL(I,J), I=1, SS_DIM)    
   ENDDO 
   CLOSE(INORRK)


!  Compute correction by applying gain matrix to innovation vector: Kr * (y - H(x^f))
   DO I=1,SS_DIM
     CRSTATE(I)=0.0_DP
     DO K=1,NLOC
        CRSTATE(I)= CRSTATE(I)+ KAL(I,K)*INNOV(K)
     ENDDO
   ENDDO  
   
!  Project this increment to the full state space(Er * Kr * (y-Hx^f)) and add it to x^f to get x^a

   DO I=1,STDIM
     RRKSUM=0.0_DP
     DO J=1,SS_DIM
       RRKSUM = RRKSUM + STSD(I)*SEOF(I,J)*CRSTATE(J)    
!      write(4500+myid,*) stsd(i),seof(i,j),crstate(j)   
     ENDDO
     STANL(I) = STFCT(I) + RRKSUM
   ENDDO


!  Transfer x^a to the model variable, El, U, and V
    IDUMMY=0
    IF(EL_ASSIM) THEN
      DO I=1, MGL
        IDUMMY = IDUMMY + 1
        RRKEL2(I)  = STANL(IDUMMY) 
      ENDDO
    ENDIF
    IF(UV_ASSIM) THEN 
      DO I=1, KBM1
        DO J=1, NGL
          IDUMMY = IDUMMY + 1
          RRKU2(J,I) = STANL(IDUMMY) 
        ENDDO
      ENDDO
      DO I=1, KBM1
        DO J=1, NGL
          IDUMMY = IDUMMY + 1
          RRKV2(J,I) = STANL(IDUMMY)
        ENDDO
      ENDDO
    ENDIF
    IF(T_ASSIM) THEN 
      DO I=1, KBM1
        DO J=1, MGL
          IDUMMY = IDUMMY + 1
          RRKT2(J,I) = STANL(IDUMMY) 
        ENDDO
      ENDDO
    ENDIF
    IF(S_ASSIM) THEN 
      DO I=1, KBM1
        DO J=1, MGL
          IDUMMY = IDUMMY + 1
          RRKS2(J,I) = STANL(IDUMMY) 
        ENDDO
      ENDDO
    ENDIF

   IF (SERIAL) THEN
     IF(EL_ASSIM) THEN
       DO I=1, MGL
         EL(I) = RRKEL2(I)
       ENDDO
       D  = H + EL
       ET = EL
       DT = D
       DO I=1, NGL
         EL1(I)=(EL(NVG(I,1)) + EL(NVG(I,2)) + EL(NVG(I,3)) )/3.0_DP
       ENDDO
       D1  = H1 + EL1
       ET1 = EL1
       DT1 = D1
     ENDIF
     IF(UV_ASSIM) THEN  
       DO K=1, KBM1
         DO I=1, NGL
           U(I,K)  =  RRKU2(I,K)
           V(I,K)  =  RRKV2(I,K) 
         ENDDO
       ENDDO
     ENDIF
     IF(T_ASSIM) THEN  
       DO K=1, KBM1
         DO I=1, MGL
           T1(I,K)  =  RRKT2(I,K)
         ENDDO
       ENDDO
     ENDIF
     IF(S_ASSIM) THEN  
       DO K=1, KBM1
         DO I=1, MGL
           S1(I,K)  =  RRKS2(I,K)
         ENDDO
       ENDDO
     ENDIF
   ELSE
#  if defined (MULTIPROCESSOR)
     IF(PAR) THEN
       IF(EL_ASSIM) THEN
         DO I=1, M
           EL(I)=RRKEL2(NGID(I))
         ENDDO
         D  = H + EL
         ET = EL
         DT = D
         DO I=1, N
           EL1(I)=(EL(NV(I,1)) + EL(NV(I,2)) + EL(NV(I,3)) )/3.0_DP
         ENDDO
         D1  = H1 + EL1
         ET1 = EL1
         DT1 = D1
       ENDIF
       IF(UV_ASSIM) THEN
         DO I=1, N
           DO J=1, KBM1
             U(I,J)=RRKU2(EGID(I),J)
             V(I,J)=RRKV2(EGID(I),J)
           ENDDO
         ENDDO
       ENDIF
       IF(T_ASSIM) THEN
         DO I=1, M
           DO J=1, KBM1
	     T1(I,J)=RRKT2(NGID(I),J)
	   ENDDO
	 ENDDO  
       ENDIF
       IF(S_ASSIM) THEN
         DO I=1, M
           DO J=1, KBM1
	     S1(I,J)=RRKS2(NGID(I),J)
	   ENDDO
	 ENDDO  
       ENDIF
     ENDIF
#  endif
   ENDIF 
  

!  Compute diagnostic on analysis 
!  Compute the analysis error at observational sites 
!   DO J=1,NLOC
!     MODDATA(J)=STANL(STLOC(J))
!   ENDDO
    CALL H_MAPPING(STANL,MODDATA)
   DO K = 1, NLOC
      ERRVEC(K) = OBSDATA(K) - MODDATA(K)
!      write(9000,*)  'obs,', OBSDATA(K) ,'mod', MODDATA(K)
   ENDDO
   TEXT='ANL(obs)'
   ERR1 = 0.0_DP
   ERR2 = 0.0_DP
   ERR3 = 0.0_DP
   DO K = 1, NLOC
     IF(ABS(ERRVEC(K)) > ERR1) THEN
        ERR1 = ABS(ERRVEC(K))
     ENDIF
     ERR2 = ERR2 + ERRVEC(K) * ERRVEC(K)
     ERR3 = ERR3 + ABS(ERRVEC(K))
   ENDDO
     ERR2 = DSQRT(ERR2/NLOC)
     ERR3 = ERR3 / NLOC
     ERR2_INN_AN1 = ERR2

!   write(1001,'(871f12.4)') (STFCT(i),i=1,871)  
!   write(1000,'(871f12.4)') (stanl(i),i=1,871)

!  Compute the analysis error of entire domain 
   DO I=1,STDIM
      ERRVEC(I)=STANL(I)-STTRUE1(I)
   ENDDO
   TEXT='ANL(full)'
   ERR1 = 0.
   ERR2 = 0.
   ERR3 = 0.
   DO K = 1, STDIM
     IF(ABS(ERRVEC(K)) > ERR1) THEN
        ERR1 = ABS(ERRVEC(K))
     ENDIF
     ERR2 = ERR2 + ERRVEC(K) * ERRVEC(K)
     ERR3 = ERR3 + ABS(ERRVEC(K))
   ENDDO
     ERR2 = DSQRT(ERR2/STDIM)
     ERR3 = ERR3 / STDIM
     ERR2_TOT_AN1 = ERR2
     
!  Write out error diagnostics for forecast and analysis

   IF(MSR) WRITE(75,'(I5,4E15.7)') ICYC,ERR2_INN_FCT, ERR2_TOT_FCT, &
                                       ERR2_INN_AN1, ERR2_TOT_AN1

   RETURN
   END SUBROUTINE RRK_RRKF



!   SUBROUTINE SETUP_RRKA 
!
!!------------------------------------------------------------------------------|
!!  SETUP RRK_ASSIMILATION RUN                                                  |
!!------------------------------------------------------------------------------|
!   
!   USE LIMS
!   USE CONTROL
!   USE MOD_RRK
!   USE RRKVAL, only: nloc
!   USE MOD_OBCS
!   IMPLICIT NONE
!
!   INTEGER I,J
!   INTEGER SS_DIM
!   INTEGER STDIM
!   
!   STDIM = 0
!   IF(EL_ASSIM) STDIM = STDIM + MGL
!   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
!   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
!   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
!
!   SS_DIM = RRK_NEOF
!
!!rrkf  Set the first seed integer for the random number generation in gasdev, 
!!rrkf  which will be used for the measurement error. 
!   IDUM = -31
!
!!READ THE OBSERVATION NUMBER AND LOCATION
!   CALL READOBS2
!   
!!rrkf  Read the one standard deviation and set observation error'
!   FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/avgstd.dat'
!   OPEN(INORRK,FILE=FILENAME,STATUS='OLD')
!   READ(INORRK,*)
!   READ(INORRK,*) ELSD, USD, TSD, SSD
!   CLOSE(INORRK)   
!   
!!rrkf  Set the magnitue of the measurement error as 1% of the spatially averaged one
!!rrkf  standard deviation.
!   CALL MAKEOBSERR   
!   
!!rrkf  Set a state vector consisting of the spatially averaged one standard deviation for later use.
!   CALL MAKESTSD   
!  
!!rrkf  Read the retained eofs from 'eof.cdf'
!   DO I=1, SS_DIM
!!      CALL READEOF2(I)
!      DO J=1,STDIM
!         SEOF(J,I)= STTEMP1(J)
!      ENDDO
!   ENDDO   
! 
!!rrkf  Read the Kalman gain matrix (Kr), which is creasted by 'fvcomrrK'. 
!   FILENAME = TRIM(OUTPUT_DIR)//'/rrktemp/rrK.dat'
!   OPEN(INORRK,FILE=FILENAME, FORM='UNFORMATTED')
!   DO J=1, NLOC
!     READ(INORRK) (KAL(I,J), I=1, SS_DIM)    
!   ENDDO 
!   CLOSE(INORRK)
!   
!!rrkf  Set the file name of the assimilation results, forecast and analysis error
!   ERRFILE = './ErrOut.dat'
!   OPEN(OUTERR,FILE=ERRFILE,STATUS='UNKNOWN')   
!
!!   I_INITIAL = IINT
!!   I_INITIAL = RRK_START_TIME
!!   NCYC      = (IEND-I_INITIAL)/DELTA_ASS
!!   N1CYC     = 1
!  
!   RETURN
!   END SUBROUTINE SETUP_RRKA



  SUBROUTINE ALLOC_RRKA

   USE LIMS
   USE CONTROL
   USE MOD_RRK
   IMPLICIT NONE
   
   INTEGER STDIM
   INTEGER SS_DIM

   STDIM = 0
   IF(EL_ASSIM) STDIM = STDIM + MGL
   IF(UV_ASSIM) STDIM = STDIM + 2*NGL*KBM1
   IF(T_ASSIM)  STDIM = STDIM + MGL*KBM1
   IF(S_ASSIM)  STDIM = STDIM + MGL*KBM1
   
   SS_DIM = RRK_NEOF

! ALLOCATE ARRYS

   IF(.not.ALLOCATED(STLOC)) ALLOCATE(STLOC(RRK_NOBSMAX))            ;STLOC     = 0 
   IF(.not.ALLOCATED(STFCT)) ALLOCATE(STFCT(STDIM))                  ;STFCT     = ZERO
   IF(.not.ALLOCATED(STANL)) ALLOCATE(STANL(STDIM))                  ;STANL     = ZERO
   IF(.not.ALLOCATED(CRSTATE)) ALLOCATE(CRSTATE(SS_DIM))               ;CRSTATE   = ZERO
   IF(.not.ALLOCATED(KAL)) ALLOCATE(KAL(SS_DIM,RRK_NOBSMAX))       ;KAL       = ZERO
   IF(.not.ALLOCATED(MODDATA)) ALLOCATE(MODDATA(RRK_NOBSMAX))          ;MODDATA   = ZERO
   IF(.not.ALLOCATED(OBSDATA)) ALLOCATE(OBSDATA(RRK_NOBSMAX))          ;OBSDATA   = ZERO
   IF(.not.ALLOCATED(INNOV)) ALLOCATE(INNOV(RRK_NOBSMAX))          ;INNOV   = ZERO
   IF(.not.ALLOCATED(ERRVEC)) ALLOCATE(ERRVEC(STDIM))                 ;ERRVEC    = ZERO
   IF(.not.ALLOCATED(R)) ALLOCATE(R(RRK_NOBSMAX,RRK_NOBSMAX))    ;R         = ZERO
   IF(.not.ALLOCATED(OBSERR)) ALLOCATE(OBSERR(STDIM))                 ;OBSERR    = ZERO
!   ALLOCATE(STTRUE(STDIM,(RRK_END-RRK_START)/DELTA_ASS+1))         ; STTRUE  = ZERO
   IF(.not.ALLOCATED(STTRUE1)) ALLOCATE(STTRUE1(STDIM))                ;STTRUE1   = ZERO
   IF(.not.ALLOCATED(STTR1)) ALLOCATE(STTR1(STDIM))                  ;STTR1     = ZERO
   IF(.not.ALLOCATED(STTR0)) ALLOCATE(STTR0(STDIM))                  ;STTR0     = ZERO
!   IF(.not.ALLOCATED(SEOF)) ALLOCATE(SEOF(STDIM,SS_DIM))            ;SEOF      = ZERO
   IF(.not.ALLOCATED(STTEMP1)) ALLOCATE(STTEMP1(STDIM))                ;STTEMP1   = ZERO
   IF(.not.ALLOCATED(STTEMP0)) ALLOCATE(STTEMP0(STDIM))                ;STTEMP0   = ZERO
!   IF(.not.ALLOCATED(STSD)) ALLOCATE(STSD(STDIM))                   ;STSD      = ZERO   
   IF(.not.ALLOCATED(RRKU2)) ALLOCATE(RRKU2(NGL,KB))                 ;RRKU2     = ZERO
   IF(.not.ALLOCATED(RRKV2)) ALLOCATE(RRKV2(NGL,KB))                 ;RRKV2     = ZERO
   IF(.not.ALLOCATED(RRKEL2)) ALLOCATE(RRKEL2(MGL))                   ;RRKEL2    = ZERO
   IF(.not.ALLOCATED(RRKT2)) ALLOCATE(RRKT2(MGL,KB))                 ;RRKT2     = ZERO
   IF(.not.ALLOCATED(RRKS2)) ALLOCATE(RRKS2(MGL,KB))                 ;RRKS2     = ZERO
   IF (EL_ASSIM) THEN
      IF(.NOT. ALLOCATED(ELETMP))  ALLOCATE(ELETMP(MGL))
   END IF
   IF(UV_ASSIM) THEN
     IF(.NOT. ALLOCATED(UETMP))   ALLOCATE(UETMP(NGL,KBM1))
     IF(.NOT. ALLOCATED(VETMP))   ALLOCATE(VETMP(NGL,KBM1))
   END IF
   IF(T_ASSIM)  THEN
   IF(.NOT. ALLOCATED(TETMP))  ALLOCATE(TETMP(MGL,KBM1))
   END IF
   IF(S_ASSIM)  THEN
      IF(.NOT. ALLOCATED(SETMP))   ALLOCATE(SETMP(MGL,KBM1))
   END IF
    IF(.NOT. ALLOCATED(wktmp)) allocate(wktmp(nloc)) 
! END ALLOCATION
   
  RETURN
  END SUBROUTINE ALLOC_RRKA

!=======================================================================|
!rrkf  Set the magnitue of the measurement error as 1% of the spatially |
!rrkf  averaged one standard deviation.                                 | 
!=======================================================================|

   SUBROUTINE MAKEOBSERR
      
    USE LIMS
    USE CONTROL
    USE MOD_RRK
    IMPLICIT NONE

    INTEGER IDUMMY
    INTEGER I,J
    
    IDUMMY = 0
    
    IF(EL_ASSIM) THEN
      DO I =1, MGL
        IDUMMY = IDUMMY + 1
        OBSERR(IDUMMY) = ELSD*DBLE(RRK_RSCALE)
      ENDDO
    ENDIF
    
    IF(UV_ASSIM) THEN
      DO J=1, 2*KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          OBSERR(IDUMMY) = USD*DBLE(RRK_RSCALE)
        ENDDO
      ENDDO
    ENDIF
    
    IF(T_ASSIM) THEN
      DO J=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          OBSERR(IDUMMY)= TSD*DBLE(RRK_RSCALE)
        ENDDO
      ENDDO
    ENDIF
    
    IF(S_ASSIM) THEN
      DO J=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          OBSERR(IDUMMY)= SSD*DBLE(RRK_RSCALE)
        ENDDO
      ENDDO
    ENDIF  
        
   RETURN
   END SUBROUTINE MAKEOBSERR

!=======================================================================!
!rrkf  Set a state vector consisting of the spatially averaged one      !
!rrkf  standard deviation for later use.                                !
!=======================================================================!
   SUBROUTINE MAKESTSD

    USE LIMS
    USE CONTROL
    USE MOD_RRK
    IMPLICIT NONE

    INTEGER IDUMMY
    INTEGER I,J    
    
    IDUMMY = 0
    
    IF(EL_ASSIM) THEN    
      DO I =1, MGL
        IDUMMY = IDUMMY + 1
        STSD(IDUMMY) = ELSD   
      ENDDO
    ENDIF

    IF(UV_ASSIM) THEN
      DO J=1, 2*KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          IF(RRK_OPTION == 0) THEN
	     STSD(IDUMMY) = USD
	  ELSE
	     STSD(IDUMMY) = USD*DSQRT(DBLE(KBM1)) 
	  ENDIF     
        ENDDO
      ENDDO
    ENDIF
    
    IF(T_ASSIM) THEN
      DO J=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          STSD(IDUMMY)= TSD
        ENDDO
      ENDDO
    ENDIF
    
    IF(S_ASSIM) THEN
      DO J=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          STSD(IDUMMY)= SSD
        ENDDO
      ENDDO
    ENDIF    

   RETURN
   END SUBROUTINE MAKESTSD

!======================================================================  
!rrkf   Get the observations  
!======================================================================
   SUBROUTINE GETOBSDATA

    USE LIMS
    USE CONTROL
    USE MOD_RRK
    IMPLICIT NONE
    
    INTEGER J
    
!    DO J=1,NLOC
!      OBSDATA(J)= STTRUE(STLOC(J),(IINT-1-RRK_START)/DELTA_ASS+1)
!    ENDDO

   RETURN
   END SUBROUTINE GETOBSDATA

!==========================================================================
!rrkf   Generate the random number of Gaussian deviation 
!rrkf   with mean of 0 and one standard deviation of 1
!==========================================================================
   SUBROUTINE GASDEV(IDUM,REV)
    
    USE LIMS
    USE CONTROL
    IMPLICIT NONE

    INTEGER IDUM
    REAL REV
    REAL REV2
!U    USES ran2
    INTEGER :: iset=0 
    REAL FAC
    REAL GSET
    REAL RSQ
    REAL v1
    REAL v2
    SAVE iset,gset

    IF(ISET==0) THEN
1     CALL RAN2(IDUM,REV2)
      V1=2.*REV2-1.
      CALL RAN2(IDUM,REV2)
      V2=2.*REV2-1.
      RSQ=V1**2+V2**2
      IF(RSQ >=1. .OR. RSQ == 0.) GOTO 1
      FAC=SQRT(-2.*LOG(RSQ)/RSQ)
      GSET=V1*FAC
      REV=V2*FAC
      ISET=1
    ELSE
      REV=GSET
      ISET=0
    ENDIF

   RETURN
   END SUBROUTINE GASDEV

!=====================================================================
  SUBROUTINE RAN2(IDUM,REV2)

    USE LIMS
    USE CONTROL
    IMPLICIT NONE
    
    INTEGER IDUM,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL AM
    REAL EPS
    REAL RNMX
    REAL REV2
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
               IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
               NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    
    IF(IDUM<=0) THEN
      IDUM=MAX(-IDUM,1)
      IDUM2=IDUM
      DO 11 J= NTAB+8,1,-1
        K=IDUM/IQ1
        IDUM=IA1*(IDUM-K*IQ1)-K*IR1
        IF (IDUM<0) IDUM=IDUM+IM1
        IF (J<=NTAB) IV(J)=IDUM
11    CONTINUE
      IY=IV(1)
    ENDIF
    K=IDUM/IQ1
    IDUM=IA1*(IDUM-K*IQ1)-K*IR1
    IF(IDUM<0) IDUM=IDUM+IM1
    K=IDUM2/IQ2
    IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
    IF(IDUM2<0) IDUM2=IDUM2+IM2
    J=1+IY/NDIV
    IY=IV(j)-IDUM2
    IV(J)=IDUM
    IF(IY<1) IY=IY+IMM1
    REV2 = MIN(AM*IY,RNMX)
    
  RETURN
  END SUBROUTINE RAN2


!======================================================================!
!                                                                      !
!======================================================================!
  SUBROUTINE READOBS2
    
   USE LIMS
   USE CONTROL
   USE MOD_RRK
   USE RRKVAL , ONLY : NLOC
   IMPLICIT NONE
   
     INTEGER ::  NUM  = 0 
     INTEGER ::  SWITCH = 0
     INTEGER ::  J,K
     INTEGER ::  IDUMMY = 0
     INTEGER ::  TMP
     CHARACTER(LEN=80) FILENAME
     CHARACTER(LEN=24) HEADINFO
     INTEGER LAY(RRK_NOBSMAX)

     FILENAME = TRIM(INPUT_DIR)//"/"//trim(casename)//"_assim_rrkf.dat"
    
     OPEN(INORRK,FILE=TRIM(FILENAME),FORM='FORMATTED')

     NLOC = 0

 100 READ(INORRK,'(A24)',END=200) HEADINFO
     IF(SWITCH/=1) THEN
       IF(HEADINFO=='!=== READ IN OBSERVATION') THEN
         SWITCH = 1
         GOTO 100
       ELSE
         GOTO 100
       ENDIF
     ENDIF 
     
     IF(TRIM(HEADINFO)=='!EL') THEN
       IF(EL_OBS) THEN
         READ(INORRK,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INORRK,*)  (STLOC(K), K=1,NLOC)	 
	 
       ENDIF

       IF(EL_ASSIM) THEN
	 IDUMMY = IDUMMY + MGL
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!UV') THEN
       IF(UV_OBS) THEN
         READ(INORRK,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC+NUM>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc+num, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INORRK,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)
	 READ(INORRK,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+NGL*(LAY(K)-1)
	 ENDDO   
	 
	 NLOC = NLOC + NUM
	 DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K-NUM)+NGL*KBM1+NGL*(LAY(K-NUM)-1)
	 ENDDO
	 
       ENDIF
       IF(UV_ASSIM) THEN  
          IDUMMY = IDUMMY + NGL*KBM1
       ENDIF
     ENDIF
     
     IF(TRIM(HEADINFO)=='!T') THEN
       IF(T_OBS) THEN
         READ(INORRK,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INORRK,*)  (STLOC(K), K=NLOC-NUM+1,NLOC)       
         READ(INORRK,*)  (LAY(K),   K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+MGL*(LAY(K)-1)
	 ENDDO   
         
       ENDIF   
       IF(T_ASSIM) THEN
         IDUMMY = IDUMMY + MGL*KBM1
       ENDIF
     ENDIF

     IF(TRIM(HEADINFO)=='!S') THEN
       IF(S_OBS) THEN
         READ(INORRK,*) NUM
	 NLOC = NLOC + NUM
         IF(NLOC>RRK_NOBSMAX) THEN
           WRITE(IPT,*) 'not enough storage for observations:', 'Nloc=', Nloc, 'Nobsmax=', RRK_NOBSMAX
           CALL PSTOP
         ENDIF
	 READ(INORRK,*)  (STLOC(K),K=NLOC-NUM+1,NLOC)        
         READ(INORRK,*)  (LAY(K),  K=NLOC-NUM+1,NLOC)
         DO K=NLOC-NUM+1, NLOC
	   STLOC(K)=STLOC(K)+IDUMMY+MGL*(LAY(K)-1)
	 ENDDO   
        
       ENDIF
       IF(S_ASSIM) THEN
         IDUMMY = IDUMMY + MGL*KBM1 
       ENDIF 
     ENDIF
     
     GOTO 100
 200 CONTINUE

     DO J=1, NLOC-1
       DO K=2, NLOC
	 IF(STLOC(K)<STLOC(J)) THEN
	   TMP = STLOC(J)
	   STLOC(J) = STLOC(K)
	   STLOC(K) = TMP 
         ENDIF
       ENDDO
     ENDDO

     CLOSE(INORRK)
    
  RETURN 
  END SUBROUTINE READOBS2  


      SUBROUTINE LSQFIT(EL,T,NUM,NCOMP,IDX,AMP2,PHA2)
!-----------------------------------------------------------------------
!     THIS SUBROUTINE IS USED FOR ANALYSIS AMPLITUDE AND PHASE BY LSQFIT
!
!     INPUT: 
!         EL(NUM)---ELEVATION (m)
!         T(NUM) ---TIME (s)
!         NUM    ---NUMBER OF TIMESERIES
!     OUTPUT:
!         AMP(6) ---AMPPLITUDE (m)
!         PHA(6) ---PHASE (deg)
!     from QXU, modified by LZG 10/25/05, 
!-----------------------------------------------------------------------
      IMPLICIT NONE
!      INTEGER, PARAMETER :: NCOMP=6, NCOMP2 = NCOMP*2+1
      REAL, PARAMETER, DIMENSION(6) :: &
                   !  S2       M2       N2       K1       P1       O1 
      PERIOD2 = (/43200.0, 44712.0, 45570.0, 86164.0, 86637.0, 92950.0/) !(sec)
              != 12.0000  12.4200  12.6583  23.9344  24.0658  25.8194 hours
      REAL, PARAMETER :: PI = 3.1415926
      
      INTEGER NCOMP, NCOMP2, IDX(6)
      INTEGER NUM,N,I,J,K,I1,I2,J1,J2
      REAL*8, DIMENSION(NUM) :: EL
      REAL*8, DIMENSION(NUM) :: T
      REAL*8  STEL,AEL
      REAL*8 A(NCOMP*2+1,NCOMP*2+1),B(NCOMP*2+1),F(NCOMP)
      REAL*8 AMP1(NCOMP),PHA1(NCOMP)
      REAL*8 AMP2(NCOMP),PHA2(NCOMP)
      REAL*8 PERIOD(NCOMP)
      
!      F = 2.0*PI/(PERIOD/3600.0)        !(1/HOUR)
!      F = 2.0*PI/PERIOD                 !(1/s)

      NCOMP2 = NCOMP*2+1
      
      N=0
      DO I=1, 6
         IF(IDX(I)==1) THEN       
           N = N + 1
	   PERIOD(N) = PERIOD2(I)
	 ENDIF
      ENDDO      
      
      F = 2.0*PI/PERIOD                  !(1/s)
     
      AEL = 0.0
      DO N=1,NUM
         AEL = AEL +EL(N)
      ENDDO
      AEL = AEL/FLOAT(NUM) 	 
      DO N=1,NUM
         EL(N)=EL(N)-AEL
      ENDDO
      STEL=0.0
      DO N=1,NUM
         STEL=STEL+EL(N)*EL(N)	 
      ENDDO
      STEL = SQRT(STEL/FLOAT(NUM))
      DO N=1,NUM
         EL(N)=EL(N)/STEL
      ENDDO

    
      DO J = 1, NCOMP2
         DO K = 1, NCOMP2
            A(J,K) = 0.0
         ENDDO
      ENDDO
      DO J = 1, NCOMP2
         B(J) = 0.0
      ENDDO
           
      DO N = 1,NUM
         A(1,1)    = A(1,1)    + 1
	 DO I=1,NCOMP
	    I1 = I*2
	    I2 = I1+1
	    A(1,I1)= A(1,I1)   + COS(F(I)*T(N))
	    A(1,I2)= A(1,I2)   + SIN(F(i)*T(N))
	 ENDDO 
	 DO I=1,NCOMP  
	    I1 = I*2
	    I2 = I1+1
	    DO J=I,NCOMP
	       J1 = J*2
	       J2 = J1+1
	       A(I1,J1) = A(I1,J1) + COS(F(I)*T(N))* COS(F(J)*T(N))
	       A(I1,J2) = A(I1,J2) + COS(F(I)*T(N))* SIN(F(J)*T(N))
	       A(I2,J1) = A(I2,J1) + SIN(F(I)*T(N))* COS(F(J)*T(N))
	       A(I2,J2) = A(I2,J2) + SIN(F(I)*T(N))* SIN(F(J)*T(N))
	    ENDDO   
	 ENDDO      
         
         B(1) = B(1) + EL(N)
	 DO I=1,NCOMP
	    I1 = I*2
	    I2 = I1+1
	    B(I1) = B(I1) + EL(N)*COS(F(I)*T(N))
	    B(I2) = B(I2) + EL(N)*SIN(F(I)*T(N))
	 ENDDO   
      ENDDO
      DO I=2,NCOMP2
         DO J=1,I
            A(I,J)=A(J,I)
         ENDDO
      ENDDO
         
      CALL GAUSSJ_2(A,NCOMP2,NCOMP2,B,1,1)
      
      DO I=1,NCOMP
         I1 = I*2
	 I2 = I1+1
         AMP1(I) = SQRT(B(I1)*B(I1)+B(I2)*B(I2))*STEL
	 PHA1(I) = ATAN2(B(I2),B(I1))*180/PI
         IF(PHA1(I).LT.0) PHA1(I)=PHA1(I)+360.0
      ENDDO  
      AMP2 = AMP1
      PHA2 = PHA1 	   
          
      RETURN
      END SUBROUTINE LSQFIT
     

!--------------------------------------------------      
      SUBROUTINE GAUSSJ_2(A,N,NP,B,M,MP)
      IMPLICIT NONE
      INTEGER M,MP,N,NP,NMAX
      REAL*8 A(NP,NP),B(NP,MP)
      PARAMETER (NMAX=50)
      INTEGER I,ICOL,IROW,J,K,L,LL,INDXC(NMAX),INDXR(NMAX),&
             IPIV(NMAX)
      REAL*8 BIG,DUM,PIVINV
      DO J=1,N
         IPIV(J)=0
      ENDDO
      DO 22 I=1,N
         BIG=0.
         DO 13 J=1,N
            IF(IPIV(J).NE.1)THEN
              DO 12 K=1,N
                 IF (IPIV(K).EQ.0) THEN
                   IF (ABS(A(J,K)).GE.BIG)THEN
                     BIG=ABS(A(J,K))
                     IROW=J
                     ICOL=K
                   ENDIF

                 ELSE IF (IPIV(K).GT.1) THEN
                   PAUSE 'SINGULAR MATRIX IN GAUSSJ'
                 ENDIF
12            CONTINUE
            ENDIF
13       CONTINUE
         IPIV(ICOL)=IPIV(ICOL)+1
         IF (IROW.NE.ICOL) THEN
           DO 14 L=1,N
              DUM=A(IROW,L)
              A(IROW,L)=A(ICOL,L)
              A(ICOL,L)=DUM
14         CONTINUE
           DO 15 L=1,M
              DUM=B(IROW,L)
              B(IROW,L)=B(ICOL,L)
              B(ICOL,L)=DUM
15         CONTINUE
         ENDIF
         INDXR(I)=IROW
         INDXC(I)=ICOL
         IF (A(ICOL,ICOL).EQ.0.) PAUSE 'SINGULAR MATRIX IN GAUSSJ'

         PIVINV=1./A(ICOL,ICOL)
         A(ICOL,ICOL)=1.
         DO 16 L=1,N
            A(ICOL,L)=A(ICOL,L)*PIVINV
16       CONTINUE
         DO 17 L=1,M
            B(ICOL,L)=B(ICOL,L)*PIVINV
17       CONTINUE
         DO 21 LL=1,N
            IF(LL.NE.ICOL)THEN
              DUM=A(LL,ICOL)
              A(LL,ICOL)=0.
              DO 18 L=1,N
                 A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18            CONTINUE
              DO 19 L=1,M
                 B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19            CONTINUE
            ENDIF
21       CONTINUE
22    CONTINUE

      DO 24 L=N,1,-1
         IF(INDXR(L).NE.INDXC(L))THEN

           DO 23 K=1,N
              DUM=A(K,INDXR(L))
              A(K,INDXR(L))=A(K,INDXC(L))
              A(K,INDXC(L))=DUM
23         CONTINUE
         ENDIF
24    CONTINUE

      RETURN
      END SUBROUTINE GAUSSJ_2
      
   
   SUBROUTINE READMEDM2
   
   USE ALL_VARS
   USE MOD_RRK
   IMPLICIT NONE

   INTEGER NUMS
   integer i,k,idump,IINTT,NGL9,MGL9
   REAL(SP) THOUR9
   CHARACTER DIR*120,FILE_NO*4,FILENAME*120
   REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::UGL,VGL,WGL,KMGL
   REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::S1GL,T1GL,RHO1GL
   REAL(SP), ALLOCATABLE, DIMENSION(:)     ::ELGL,UAGL,VAGL

!   if(MOD(IEND,DELTA_ASS) .ne. 0) then
!      print*,'-------Error in read ture obs data'
!      print*,'iint,int(iint/DELTA_ASS)*DELTA_ASS=',iint,int(iint/DELTA_ASS)*DELTA_ASS
!      stop
!   endif
   
!=== read true state of current time step to sttr0
   
!   WRITE(FILE_NO,'(I4.4)') ISTART/24
!   DIR = '../medm_bck/'
!   FILENAME='chn_sim'//FILE_NO//'.dat'
!   OPEN(1,file=TRIM(DIR)//TRIM(FILENAME),STATUS='OLD',FORM='UNFORMATTED')
!
!   READ(1) IINTT,NGL9,MGL9,THOUR9  
!
!     DO I=1,NGL
!        READ(1) (UGL(I,K),VGL(I,K),WGL(I,K),KMGL(I,K), K=1,KBM1)
!     ENDDO
!     DO I=1,MGL
!        READ(1) ELGL(I),(T1GL(I,K),S1GL(I,K),RHO1GL(I,K),K=1,KBM1)
!     ENDDO     

    CALL READRESTART(NC_START,REF_START_TIME)
!    CALL GR2ST(0)
!    STTR0 = STTEMP0


!=== read true state of next time step to sttr1
     
!   WRITE(FILE_NO,'(I4.4)') ISTART/24+1
!   DIR = '../medm_bck/'
!   FILENAME='chn_sim'//FILE_NO//'.dat'
!   OPEN(1,file=TRIM(DIR)//TRIM(FILENAME),STATUS='OLD',FORM='UNFORMATTED')

     
     DEALLOCATE(UGL,VGL,WGL,KMGL)
     DEALLOCATE(ELGL,UAGL,VAGL)
     DEALLOCATE(S1GL,T1GL,RHO1GL)
     CLOSE(1)
     RETURN
     
     END SUBROUTINE READMEDM2
!------------------------
   SUBROUTINE OBS2VEC
   use mod_rrkf_obs
   IMPLICIT NONE
   REAL(DP) ::H_VEC(NLOC)
   REAL(DP):: VEC(STDIM)
   real    :: var1,var2,w1,w2
   INTEGER :: I,K,NLAY,IDUMMY
   call VEC2VAR(VEC)
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1,N_ASSIM_EL
       NLAY      = EL0_OBS(I)%N_LAYERS
       if (nlay .ne. 1) then
           print *, 'impossible nlay for elevation is ',nlay
           call pstop
       end if
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= EL0_OBS(I)%EL(I,K)
         WKTMP(IDUMMY)= OBSERR_EL
       END DO
     END DO
    END IF

   IF(UV_ASSIM) THEN
     DO I=1,N_ASSIM_CUR
       NLAY      = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= CUR_OBS(I)%UO(1,K) ! only one time , i.e. current time
         WKTMP(IDUMMY)= OBSERR_UV
       END DO
     END DO
     DO I=1,N_ASSIM_CUR
       NLAY      = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= CUR_OBS(I)%VO(1,K) ! only one time , i.e. current time
         WKTMP(IDUMMY)= OBSERR_UV
       END DO
     END DO
   END IF



   IF(T_ASSIM) THEN

     DO I=1,N_ASSIM_T
       NLAY      = T0_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= T0_OBS(I)%TEMP(i,K) !only one time , i.e. current time
         WKTMP(IDUMMY)= OBSERR_T
       END DO
     END DO
   END IF


   IF(S_ASSIM) THEN
     DO I=1,N_ASSIM_S
       NLAY      = S0_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         OBSDATA(IDUMMY)= S0_OBS(I)%SAL(1,K) !only one time , i.e. current time
         WKTMP(IDUMMY)= OBSERR_S
       END DO
     END DO
   END IF
   END SUBROUTINE OBS2VEC
!====================================================
!------------pengfei-----------------------------------




   SUBROUTINE H_MAPPING(VEC,H_VEC)
   use rrkval
   use mod_rrkf_obs
   IMPLICIT NONE
   REAL(DP) ::H_VEC(NLOC)
   REAL(DP):: VEC(STDIM)
   real    :: var1,var2,w1,w2
   INTEGER :: I,K,NLAY,IDUMMY
   H_VEC=ZERO_DB
   call VEC2VAR(VEC)
   IDUMMY = 0
   IF(EL_ASSIM) THEN
     DO I=1,N_ASSIM_EL
       NLAY      = EL0_OBS(I)%N_LAYERS
       if (nlay .ne. 1) then
           print *, 'impossible nlay for elevation is ',nlay
           call pstop
       end if
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = ELETMP(EL0_OBS(I)%INNODE)
!         var2 = ELETMP(EL0_OBS(I)%INNODE)
!         W1 = EL0_OBS(I)%S_WEIGHT(K,1)
!         W2 = EL0_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1 !*W1 + var2*W2
       END DO
     END DO
    END IF

   IF(UV_ASSIM) THEN
     DO I=1,N_ASSIM_CUR
       NLAY      = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = UETMP(CUR_OBS(I)%N_CELL,CUR_OBS(I)%S_INT(K,1))
         var2 = UETMP(CUR_OBS(I)%N_CELL,CUR_OBS(I)%S_INT(K,2))
         W1 =   CUR_OBS(I)%S_WEIGHT(K,1)
         W2 =   CUR_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1*W1 + var2*W2
       END DO
     END DO
     DO I=1,N_ASSIM_CUR
       NLAY      = CUR_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = VETMP(CUR_OBS(I)%N_CELL,CUR_OBS(I)%S_INT(K,1))
         var2 = VETMP(CUR_OBS(I)%N_CELL,CUR_OBS(I)%S_INT(K,2))
         W1 =   CUR_OBS(I)%S_WEIGHT(K,1)
         W2 =   CUR_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1*W1 + var2*W2
       END DO
     END DO
   END IF



   IF(T_ASSIM) THEN

     DO I=1,N_ASSIM_T
       NLAY      = T0_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = TETMP(T0_OBS(I)%INNODE,T0_OBS(I)%S_INT(K,1))
         var2 = TETMP(T0_OBS(I)%INNODE,T0_OBS(I)%S_INT(K,2))
         W1 =   T0_OBS(I)%S_WEIGHT(K,1)
         W2 =   T0_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1*W1 + var2*W2
       END DO
     END DO
   END IF


   IF(S_ASSIM) THEN
     DO I=1,N_ASSIM_S
       NLAY      = S0_OBS(I)%N_LAYERS
       DO K=1,NLAY
         IDUMMY = IDUMMY + 1
         var1 = SETMP(S0_OBS(I)%INNODE,S0_OBS(I)%S_INT(K,1))
         var2 = SETMP(S0_OBS(I)%INNODE,S0_OBS(I)%S_INT(K,2))
         W1 =   S0_OBS(I)%S_WEIGHT(K,1)
         W2 =   S0_OBS(I)%S_WEIGHT(K,2)
         H_VEC(IDUMMY)= var1*W1 + var2*W2
       END DO
     END DO
   END IF
   END SUBROUTINE H_MAPPING
!------------pengfei-----------------------------------
   SUBROUTINE VAR2VEC(VEC)
   USE LIMS
   USE CONTROL
   IMPLICIT NONE
   INTEGER :: IDUMMY,I,K
   REAL(DP):: VEC(STDIM)
!   ALLOCATE(VEC(STDIM))
   IDUMMY = 0
   IF(EL_ASSIM) THEN
      DO I=1, MGL
         IDUMMY = IDUMMY + 1
         VEC(IDUMMY) = ELETMP(I)
      ENDDO

   ENDIF

   IF(UV_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          VEC(IDUMMY) = UETMP(I,K) 
        ENDDO
      ENDDO


      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          VEC(IDUMMY) = VETMP(I,K)
        ENDDO
      ENDDO

   ENDIF

   IF(T_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          VEC(IDUMMY) = TETMP(I,K)
        ENDDO
      ENDDO
   ENDIF

   IF(S_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          VEC(IDUMMY) = SETMP(I,K)
        ENDDO
      ENDDO
   ENDIF
!   DEALLOCATE(VEC)
   END SUBROUTINE VAR2VEC

!------------------------------------------------------

!------------------------------------------------------
   SUBROUTINE VEC2VAR(VEC)
   USE LIMS
   USE CONTROL
   IMPLICIT NONE
   INTEGER :: IDUMMY,I,K
   REAL(DP):: VEC(STDIM)
!   ALLOCATE(VEC(STDIM))
   IDUMMY = 0
   IF(EL_ASSIM) THEN
      DO I=1, MGL
         IDUMMY = IDUMMY + 1
         ELETMP(I) = VEC(IDUMMY)
      ENDDO

   ENDIF

   IF(UV_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          UETMP(I,K) = VEC(IDUMMY)
        ENDDO
      ENDDO


      DO K=1, KBM1
        DO I=1, NGL
          IDUMMY = IDUMMY + 1
          VETMP(I,K) = VEC(IDUMMY)
        ENDDO
      ENDDO

   ENDIF

   IF(T_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          TETMP(I,K) = VEC(IDUMMY)
        ENDDO
      ENDDO
   ENDIF

   IF(S_ASSIM) THEN
      DO K=1, KBM1
        DO I=1, MGL
          IDUMMY = IDUMMY + 1
          SETMP(I,K) = VEC(IDUMMY)
        ENDDO
      ENDDO
   ENDIF
!   DEALLOCATE(VEC)
   END SUBROUTINE VEC2VAR


!====================================================
# endif     
END MODULE MOD_RRKA
