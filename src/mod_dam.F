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

MODULE MOD_DAM
# if defined (THIN_DAM)
   USE ALL_VARS, ONLY : nvg,EL,ELF,H,ART1,MT,NT,ZERO,D,D1,DT,H,H1,ELF1,DZ,KB,NV
   USE ALL_VARS, ONLY : xc,yc,vx,vy,alpha,nbe,ne,IENODE,ISONB,IEC
   USE MOD_PREC
   USE MOD_PAR
   USE control, only : par,serial,msr,casename,input_dir
   USE lims,    only : m,myid
   USE Mod_Utils, only: pstop
# if defined (SPHERICAL)
   USE MOD_SPHERICAL
# endif

   IMPLICIT NONE
   SAVE
   REAL(SP), PARAMETER  :: EL_EPS = 1.0E-04_SP
   character(len=120)   :: cellfile
   character(len=120)   :: nodefile

   INTEGER :: NODE_DAM1_GL,NODE_DAM1_N
   INTEGER :: NODE_DAM2_GL,NODE_DAM2_N
   INTEGER :: NODE_DAM3_GL,NODE_DAM3_N

   INTEGER :: CELL_DAM_GL,CELL_DAM_N

   
   INTEGER ,ALLOCATABLE :: I_NODE_DAM1_GL(:,:),I_NODE_DAM1_N(:,:)
   INTEGER ,ALLOCATABLE :: I_NODE_DAM2_GL(:,:),I_NODE_DAM2_N(:,:)
   INTEGER ,ALLOCATABLE :: I_NODE_DAM3_GL(:,:),I_NODE_DAM3_N(:,:)

   REAL(SP) ,ALLOCATABLE :: DAM1_SPONGE_GL(:,:),DAM1_SPONGE_N(:,:)
   REAL(SP) ,ALLOCATABLE :: DAM2_SPONGE_GL(:,:),DAM2_SPONGE_N(:,:)
   REAL(SP) ,ALLOCATABLE :: DAM3_SPONGE_GL(:,:),DAM3_SPONGE_N(:,:)
   

   INTEGER ,ALLOCATABLE :: I_CELL_DAM_GL(:,:),I_CELL_DAM_N(:,:)
   INTEGER ,ALLOCATABLE :: CLP_CELL(:)  ! Clamping cell which shares
                                        ! two egdes with two dam cell
   INTEGER ,ALLOCATABLE :: CLP_KDAM(:,:)
   REAL(SP),ALLOCATABLE :: CLP_ALPHA(:) ! Angle between Clamping cell
                                        ! and dam boundary 

   REAL(SP),ALLOCATABLE :: D_DAM1_GL(:,:),D_DAM2_GL(:,:),D_DAM3_GL(:,:)
   REAL(SP),ALLOCATABLE :: D_DAM1_N(:,:),D_DAM2_N(:,:),D_DAM3_N(:,:)

   INTEGER ,ALLOCATABLE :: NODE_DAM1(:,:),NODE_DAM2(:,:),NODE_DAM3(:,:)
   INTEGER ,ALLOCATABLE :: CELL_DAM(:,:)
   REAL(SP),ALLOCATABLE :: D_DAM1(:,:),D_DAM2(:,:),D_DAM3(:,:)
   INTEGER              :: NN_DAM1,NN_DAM2,NN_DAM3,NC_DAM

   REAL(SP), ALLOCATABLE :: DAM1_R_SPG(:), DAM1_C_SPG(:)
   REAL(SP), ALLOCATABLE :: DAM2_R_SPG(:), DAM2_C_SPG(:)
   REAL(SP), ALLOCATABLE :: DAM3_R_SPG(:), DAM3_C_SPG(:)
   
   REAL(SP), ALLOCATABLE :: A1U_DAM(:,:)      
   REAL(SP), ALLOCATABLE :: A2U_DAM(:,:)     
   REAL(SP), ALLOCATABLE :: AW0_DAM(:,:)
   REAL(SP), ALLOCATABLE :: AWX_DAM(:,:)
   REAL(SP), ALLOCATABLE :: AWY_DAM(:,:) 
   
   REAL(SP), ALLOCATABLE :: KDAM(:),KDAM1(:)
   REAL(SP), ALLOCATABLE :: DAM_SPONGE(:)  

   INTEGER,  ALLOCATABLE :: NBE_DAM(:)
   INTEGER,  ALLOCATABLE :: IS_DAM(:)
   INTEGER,  ALLOCATABLE :: N_DAM_MATCH(:,:)  ! Nodes connected at
   ! one point 
   INTEGER,  ALLOCATABLE :: E_DAM_MATCH(:)  ! Matching cell against
   ! one dam cell.
   INTEGER,  ALLOCATABLE :: EDGE_MATCH(:)  ! Matching edges against
   INTEGER,  ALLOCATABLE :: EDGE_DAM(:,:)  ! one dam cell.
#  if defined (SPHERICAL)
   REAL(SP), ALLOCATABLE :: DLTXNE_DAM_MATCH(:)
   REAL(SP), ALLOCATABLE :: DLTYNE_DAM_MATCH(:)
#  endif
   logical               :: fexist

   integer               :: nn,corner_proc_id
   CONTAINS
!==============================================================================|
!
!==============================================================================|
   SUBROUTINE READ_DAM
   IMPLICIT NONE
   INTEGER :: I
   INTEGER :: NCNT,I1,I2,I3,I4,J1,J2,JJ,IA,IB
   INTEGER :: CLP_COUNT,COUNT,ITMP,ITMP1,ITMP2,ITMP3,ITMP4
  
   INTEGER, ALLOCATABLE :: TEMP1(:,:)
   REAL,    ALLOCATABLE :: TEMP2(:,:),TEMP3(:),TEMP4(:)

   REAL(SP)  TEMP,DTMP,C_SPONGE,ALPHA1
   REAL(DP) X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: read_dam" 

  !ensure sediment parameter file exists
  cellfile = "./"//trim(input_dir)//"/"//trim(casename)//'_dam_cell.dat'
  nodefile = "./"//trim(input_dir)//"/"//trim(casename)//'_dam_node.dat'
  inquire(file=trim(cellfile),exist=fexist)
  if(.not.fexist)then
    write(*,*)'dam cell file: ',trim(cellfile),' does not exist'
    write(*,*)'stopping'
    call pstop
  end if

  inquire(file=trim(nodefile),exist=fexist)
  if(.not.fexist)then
    write(*,*)'dam node file: ',trim(nodefile),' does not exist'
    write(*,*)'stopping'
    call pstop
  end if

   OPEN(1,FILE=trim(nodefile))

!---read in type 1 dam. ---------------------
   READ(1,*)
   READ(1,*) NN_DAM1
   
   ALLOCATE(NODE_DAM1(NN_DAM1,2));     NODE_DAM1 = 0
   ALLOCATE(D_DAM1(NN_DAM1,2));        D_DAM1    = 0
   ALLOCATE(DAM1_SPONGE_GL(NN_DAM1,2));DAM1_SPONGE_GL = 0.0

   DO I=1,NN_DAM1
     READ(1,*) NODE_DAM1(I,1),NODE_DAM1(I,2),D_DAM1(I,1),D_DAM1(I,2)
   END DO

!---read in type 2 dam. ---------------------   
   READ(1,*)
   READ(1,*) NN_DAM2
   
   ALLOCATE(NODE_DAM2(NN_DAM2,3));    NODE_DAM2 = 0
   ALLOCATE(D_DAM2(NN_DAM2,3));       D_DAM2    = 0
   ALLOCATE(DAM2_SPONGE_GL(NN_DAM1,2));DAM2_SPONGE_GL = 0.0
   
   DO I=1,NN_DAM2
     READ(1,*) NODE_DAM2(I,1),NODE_DAM2(I,2),NODE_DAM2(I,3),    &
               D_DAM2(I,1),D_DAM2(I,2),D_DAM2(I,3)
   END DO
   
!---read in type 3 dam. ---------------------   
   READ(1,*)
   READ(1,*) NN_DAM3
   
   ALLOCATE(NODE_DAM3(NN_DAM2,4));    NODE_DAM3 = 0
   ALLOCATE(D_DAM3(NN_DAM3,4));       D_DAM3    = 0
   ALLOCATE(DAM3_SPONGE_GL(NN_DAM1,2));DAM3_SPONGE_GL = 0.0
   
   DO I=1,NN_DAM3
     READ(1,*) NODE_DAM3(I,1),NODE_DAM3(I,2),NODE_DAM3(I,3),NODE_DAM3(I,4), &
               D_DAM3(I,1),D_DAM3(I,2),D_DAM3(I,3),D_DAM3(I,4)
   END DO
   
   CLOSE(1)
!---read in cells list. ---------------------     
   OPEN(1,FILE=trim(cellfile))
   READ(1,*) NC_DAM
   
   ALLOCATE(CELL_DAM(NC_DAM,2));        CELL_DAM = ZERO
   DO I=1,NC_DAM
     READ(1,*) CELL_DAM(I,1),CELL_DAM(I,2)
   END DO
   CLOSE(1)
!
!---Map to Local Domain----------------------

   IF(SERIAL)THEN
     NODE_DAM1_GL = NN_DAM1
     NODE_DAM2_GL = NN_DAM2
     NODE_DAM3_GL = NN_DAM3
     NODE_DAM1_N  = NN_DAM1
     NODE_DAM2_N  = NN_DAM2
     NODE_DAM3_N  = NN_DAM3
     CELL_DAM_GL  = NC_DAM
     CELL_DAM_N   = NC_DAM

     ALLOCATE(I_NODE_DAM1_GL(NODE_DAM1_GL,2))
     ALLOCATE(I_NODE_DAM2_GL(NODE_DAM2_GL,3))
     ALLOCATE(I_NODE_DAM3_GL(NODE_DAM3_GL,4))

     ALLOCATE(I_NODE_DAM1_N(NODE_DAM1_N,2))
     ALLOCATE(I_NODE_DAM2_N(NODE_DAM2_N,3))
     ALLOCATE(I_NODE_DAM3_N(NODE_DAM3_N,4))

     ALLOCATE(D_DAM1_GL(NODE_DAM1_GL,2))
     ALLOCATE(D_DAM2_GL(NODE_DAM2_GL,3))
     ALLOCATE(D_DAM3_GL(NODE_DAM3_GL,4))

     ALLOCATE(D_DAM1_N(NODE_DAM1_N,2))
     ALLOCATE(D_DAM2_N(NODE_DAM2_N,3))
     ALLOCATE(D_DAM3_N(NODE_DAM3_N,4))

     ALLOCATE(I_CELL_DAM_GL(CELL_DAM_GL,2))
     ALLOCATE(I_CELL_DAM_N(CELL_DAM_N,2))

     ALLOCATE(DAM1_R_SPG(NODE_DAM1_N))
     ALLOCATE(DAM2_R_SPG(NODE_DAM2_N))
     ALLOCATE(DAM3_R_SPG(NODE_DAM3_N))
     ALLOCATE(DAM1_C_SPG(NODE_DAM1_N))
     ALLOCATE(DAM2_C_SPG(NODE_DAM2_N))
     ALLOCATE(DAM3_C_SPG(NODE_DAM3_N))



     I_NODE_DAM1_GL = NODE_DAM1
     I_NODE_DAM2_GL = NODE_DAM2
     I_NODE_DAM3_GL = NODE_DAM3

     I_NODE_DAM1_N = NODE_DAM1
     I_NODE_DAM2_N = NODE_DAM2
     I_NODE_DAM3_N = NODE_DAM3
    
     D_DAM1_GL = D_DAM1
     D_DAM2_GL = D_DAM2
     D_DAM3_GL = D_DAM3

     D_DAM1_N  = D_DAM1
     D_DAM2_N  = D_DAM2
     D_DAM3_N  = D_DAM3

     I_CELL_DAM_GL = CELL_DAM
     I_CELL_DAM_N  = CELL_DAM

     DAM1_R_SPG = DAM1_SPONGE_GL(:,1)
     DAM2_R_SPG = DAM1_SPONGE_GL(:,1)
     DAM3_R_SPG = DAM1_SPONGE_GL(:,1)
     DAM1_C_SPG = DAM1_SPONGE_GL(:,2)
     DAM2_C_SPG = DAM1_SPONGE_GL(:,2)
     DAM3_C_SPG = DAM1_SPONGE_GL(:,2)

   END IF


#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       NODE_DAM1_GL = NN_DAM1
       NODE_DAM2_GL = NN_DAM2
       NODE_DAM3_GL = NN_DAM3

       CELL_DAM_GL  = NC_DAM      
       ALLOCATE(I_NODE_DAM1_GL(NODE_DAM1_GL,2))
       ALLOCATE(I_NODE_DAM2_GL(NODE_DAM2_GL,3))
       ALLOCATE(I_NODE_DAM3_GL(NODE_DAM3_GL,4))

       ALLOCATE(D_DAM1_GL(NODE_DAM1_GL,2))
       ALLOCATE(D_DAM2_GL(NODE_DAM2_GL,3))
       ALLOCATE(D_DAM3_GL(NODE_DAM3_GL,4))
       ALLOCATE(I_CELL_DAM_GL(CELL_DAM_GL,2))

       I_NODE_DAM1_GL = NODE_DAM1
       I_NODE_DAM2_GL = NODE_DAM2
       I_NODE_DAM3_GL = NODE_DAM3

       D_DAM1_GL = D_DAM1
       D_DAM2_GL = D_DAM2
       D_DAM3_GL = D_DAM3

       I_CELL_DAM_GL = CELL_DAM
 
!--------------------Type 1 dam node-----------------------
       ALLOCATE(TEMP1(NODE_DAM1_GL,2))
       ALLOCATE(TEMP2(NODE_DAM1_GL,2))
       ALLOCATE(TEMP3(NODE_DAM1_GL))
       ALLOCATE(TEMP4(NODE_DAM1_GL))
       NCNT = 0
       DO I=1,NODE_DAM1_GL
         I1=NLID(I_NODE_DAM1_GL(I,1))
         I2=NLID(I_NODE_DAM1_GL(I,2))
	 IF(I1 /= 0.AND.I2 /= 0)THEN
           print*,"dam1 node myid:",myid,I_NODE_DAM1_GL(I,1),I_NODE_DAM1_GL(I,2)
	   NCNT = NCNT + 1
	   TEMP1(NCNT,1) = I1
	   TEMP1(NCNT,2) = I2
           TEMP2(NCNT,1) = D_DAM1_GL(I,1)
           TEMP2(NCNT,2) = D_DAM1_GL(I,2)
           TEMP3(NCNT) = DAM1_SPONGE_GL(I,1)
           TEMP4(NCNT) = DAM1_SPONGE_GL(I,2)
	 END IF
       END DO
       NODE_DAM1_N = NCNT
  
       IF(NODE_DAM1_N > 0)THEN      
         ALLOCATE(I_NODE_DAM1_N(NODE_DAM1_N,2))
         ALLOCATE(D_DAM1_N(NODE_DAM1_N,2))
         ALLOCATE(DAM1_R_SPG(NODE_DAM1_N))
         ALLOCATE(DAM1_C_SPG(NODE_DAM1_N))
	 I_NODE_DAM1_N(1:NODE_DAM1_N,:) = TEMP1(1:NODE_DAM1_N,:)
         D_DAM1_N(1:NODE_DAM1_N,:)      = TEMP2(1:NODE_DAM1_N,:)
         DAM1_R_SPG(1:NODE_DAM1_N)      = TEMP3(1:NODE_DAM1_N)
         DAM1_C_SPG(1:NODE_DAM1_N)      = TEMP4(1:NODE_DAM1_N)
       END IF
       
       DEALLOCATE(TEMP1)
       DEALLOCATE(TEMP2)
       DEALLOCATE(TEMP3)
       DEALLOCATE(TEMP4)
  
!---------------------Type 2 dam node-----------------------
       ALLOCATE(TEMP1(NODE_DAM2_GL,3))
       ALLOCATE(TEMP2(NODE_DAM2_GL,3))
       ALLOCATE(TEMP3(NODE_DAM2_GL))
       ALLOCATE(TEMP4(NODE_DAM2_GL))
       NCNT = 0
       DO I=1,NODE_DAM2_GL
         I1=NLID(I_NODE_DAM2_GL(I,1))
         I2=NLID(I_NODE_DAM2_GL(I,2))
         I3=NLID(I_NODE_DAM2_GL(I,3))
	 IF(I1 /= 0.AND.I2 /= 0.AND.I3 /= 0)THEN
           print*,"dam2 node myid:",myid,I_NODE_DAM2_GL(I,1)&
                &,I_NODE_DAM2_GL(I,2),I_NODE_DAM2_GL(I,3)
	   NCNT = NCNT + 1
	   TEMP1(NCNT,1) = I1
	   TEMP1(NCNT,2) = I2
	   TEMP1(NCNT,3) = I3
           TEMP2(NCNT,1) = D_DAM2_GL(I,1)
           TEMP2(NCNT,2) = D_DAM2_GL(I,2)
           TEMP2(NCNT,3) = D_DAM2_GL(I,3)
           TEMP3(NCNT) = DAM2_SPONGE_GL(I,1)
           TEMP4(NCNT) = DAM2_SPONGE_GL(I,2)
	 END IF
       END DO
       NODE_DAM2_N = NCNT

       IF(NODE_DAM2_N > 0)THEN      
         ALLOCATE(I_NODE_DAM2_N(NODE_DAM2_N,3))
         ALLOCATE(D_DAM2_N(NODE_DAM2_N,3))
         ALLOCATE(DAM2_R_SPG(NODE_DAM2_N))
         ALLOCATE(DAM2_C_SPG(NODE_DAM2_N))
	 I_NODE_DAM2_N(1:NODE_DAM2_N,:) = TEMP1(1:NODE_DAM2_N,:)
         D_DAM2_N(1:NODE_DAM2_N,:)     = TEMP2(1:NODE_DAM2_N,:)
         DAM2_R_SPG(1:NODE_DAM2_N)      = TEMP3(1:NODE_DAM2_N)
         DAM2_C_SPG(1:NODE_DAM2_N)      = TEMP4(1:NODE_DAM2_N)
       END IF
       
       DEALLOCATE(TEMP1)
       DEALLOCATE(TEMP2)
       DEALLOCATE(TEMP3)
       DEALLOCATE(TEMP4)

!---------------------Type 3 dam node-----------------------
       ALLOCATE(TEMP1(NODE_DAM3_GL,4))
       ALLOCATE(TEMP2(NODE_DAM3_GL,4))
       ALLOCATE(TEMP3(NODE_DAM3_GL))
       ALLOCATE(TEMP4(NODE_DAM3_GL))
       NCNT = 0
       DO I=1,NODE_DAM3_GL
         I1=NLID(I_NODE_DAM3_GL(I,1))
         I2=NLID(I_NODE_DAM3_GL(I,2))
         I3=NLID(I_NODE_DAM3_GL(I,3))
         I4=NLID(I_NODE_DAM3_GL(I,4))
	 IF(I1 /= 0.AND.I2 /= 0.AND.I3 /= 0.AND.I4 /= 0)THEN
	   NCNT = NCNT + 1
	   TEMP1(NCNT,1) = I1
	   TEMP1(NCNT,2) = I2
	   TEMP1(NCNT,3) = I3
	   TEMP1(NCNT,4) = I4
           TEMP2(NCNT,1) = D_DAM3_GL(I,1)
           TEMP2(NCNT,2) = D_DAM3_GL(I,2)
           TEMP2(NCNT,3) = D_DAM3_GL(I,3)
           TEMP2(NCNT,4) = D_DAM3_GL(I,4)
           TEMP3(NCNT) = DAM3_SPONGE_GL(I,1)
           TEMP4(NCNT) = DAM3_SPONGE_GL(I,2)
	 END IF
       END DO
       NODE_DAM3_N = NCNT

       IF(NODE_DAM3_N > 0)THEN      
         ALLOCATE(I_NODE_DAM3_N(NODE_DAM3_N,4))
         ALLOCATE(D_DAM3_N(NODE_DAM3_N,4))
         ALLOCATE(DAM3_R_SPG(NODE_DAM3_N))
         ALLOCATE(DAM3_C_SPG(NODE_DAM3_N))
	 I_NODE_DAM3_N(1:NODE_DAM3_N,:) = TEMP1(1:NODE_DAM3_N,:)
         D_DAM3_N(1:NODE_DAM3_N,:)     = TEMP2(1:NODE_DAM3_N,:)
         DAM3_R_SPG(1:NODE_DAM3_N)      = TEMP3(1:NODE_DAM3_N)
         DAM3_C_SPG(1:NODE_DAM3_N)      = TEMP4(1:NODE_DAM3_N)
       END IF
       
       DEALLOCATE(TEMP1)
       DEALLOCATE(TEMP2)
       DEALLOCATE(TEMP3)
       DEALLOCATE(TEMP4)
!--------------------- dam cell -----------------------
       ALLOCATE(TEMP1(CELL_DAM_GL,2))
       NCNT = 0
       DO I=1,CELL_DAM_GL
         I1=ELID(I_CELL_DAM_GL(I,1))
         I2=ELID(I_CELL_DAM_GL(I,2))
	 IF(I1 /= 0.AND.I2 /= 0)THEN
	   NCNT = NCNT + 1
	   TEMP1(NCNT,1) = I1
	   TEMP1(NCNT,2) = I2
	 END IF
       END DO
       CELL_DAM_N = NCNT

       IF(CELL_DAM_N > 0)THEN      
         ALLOCATE(I_CELL_DAM_N(CELL_DAM_N,2))
	 I_CELL_DAM_N(1:CELL_DAM_N,:) = TEMP1(1:CELL_DAM_N,:)
       END IF
       
       DEALLOCATE(TEMP1)

     END IF
#   endif

   ALLOCATE(A1U_DAM(0:NT,4))         ;A1U_DAM = ZERO
   ALLOCATE(A2U_DAM(0:NT,4))         ;A2U_DAM = ZERO 
   ALLOCATE(AWX_DAM(0:NT,3))         ;AWX_DAM = ZERO 
   ALLOCATE(AWY_DAM(0:NT,3))         ;AWY_DAM = ZERO 
   ALLOCATE(AW0_DAM(0:NT,3))         ;AW0_DAM = ZERO 

   ALLOCATE(KDAM(0:MT))              ;KDAM    = 0 
   ALLOCATE(KDAM1(0:NT))             ;KDAM1   = 0 
   ALLOCATE(NBE_DAM(0:NT))           ;NBE_DAM = 0

   ALLOCATE(DAM_SPONGE(0:NT))        ;DAM_SPONGE = 0

   ALLOCATE(CLP_CELL(0:NT))          ;CLP_CELL = 0
   ALLOCATE(CLP_KDAM(0:NT,3))        ;CLP_KDAM = 0
   ALLOCATE(CLP_ALPHA(0:NT))         ;CLP_ALPHA = ZERO

   ALLOCATE(IS_DAM(0:MT))            ;IS_DAM = 0
   ALLOCATE(N_DAM_MATCH(0:MT,4))     ;N_DAM_MATCH = 0
   ALLOCATE(E_DAM_MATCH(0:NT))       ;E_DAM_MATCH = 0

#  if defined (SPHERICAL)
   ALLOCATE(DLTXNE_DAM_MATCH(0:NE))       ;DLTXNE_DAM_MATCH = 0.0
   ALLOCATE(DLTYNE_DAM_MATCH(0:NE))       ;DLTYNE_DAM_MATCH = 0.0
#  endif

   ALLOCATE(EDGE_DAM(CELL_DAM_N,2))  ;EDGE_DAM = 0
   DO I=1,NE
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     DO NN=1,CELL_DAM_N
        COUNT = 0
        IF(J1==NV(I_CELL_DAM_N(NN,1),1))COUNT = COUNT + 1
        IF(J1==NV(I_CELL_DAM_N(NN,1),2))COUNT = COUNT + 1
        IF(J1==NV(I_CELL_DAM_N(NN,1),3))COUNT = COUNT + 1
        IF(J2==NV(I_CELL_DAM_N(NN,1),1))COUNT = COUNT + 1
        IF(J2==NV(I_CELL_DAM_N(NN,1),2))COUNT = COUNT + 1
        IF(J2==NV(I_CELL_DAM_N(NN,1),3))COUNT = COUNT + 1
        IF(COUNT==2.AND.ISONB(J1)==1.AND.ISONB(J2)==1)EDGE_DAM(NN,1)=I
        COUNT = 0
        IF(J1==NV(I_CELL_DAM_N(NN,2),1))COUNT = COUNT + 1
        IF(J1==NV(I_CELL_DAM_N(NN,2),2))COUNT = COUNT + 1
        IF(J1==NV(I_CELL_DAM_N(NN,2),3))COUNT = COUNT + 1
        IF(J2==NV(I_CELL_DAM_N(NN,2),1))COUNT = COUNT + 1
        IF(J2==NV(I_CELL_DAM_N(NN,2),2))COUNT = COUNT + 1
        IF(J2==NV(I_CELL_DAM_N(NN,2),3))COUNT = COUNT + 1
        IF(COUNT==2.AND.ISONB(J1)==1.AND.ISONB(J2)==1)EDGE_DAM(NN,2)=I
     END DO
   END DO


   DO I=1,NT
     DO NN=1,CELL_DAM_N
       IF(I==I_CELL_DAM_N(NN,1))E_DAM_MATCH(I)=I_CELL_DAM_N(NN,2)
       IF(I==I_CELL_DAM_N(NN,2))E_DAM_MATCH(I)=I_CELL_DAM_N(NN,1)
     END DO
   END DO

#  if defined (SPHERICAL)
   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     IF(IB==0.AND.E_DAM_MATCH(IA)/=0)THEN
       DO NN=1,NE
          IF(IEC(NN,1)==E_DAM_MATCH(IA).AND.IEC(NN,2)==0)THEN
            DLTXNE_DAM_MATCH(I)=DLTXNE(NN,1)
            DLTYNE_DAM_MATCH(I)=DLTYNE(NN,1)
          END IF
       END DO 
     END IF
   END DO
#  endif

   DO I=1,NODE_DAM1_N
     IS_DAM(I_NODE_DAM1_N(I,1))=1
     IS_DAM(I_NODE_DAM1_N(I,2))=1
     N_DAM_MATCH(I_NODE_DAM1_N(I,1),1)=1
     N_DAM_MATCH(I_NODE_DAM1_N(I,2),1)=1

     N_DAM_MATCH(I_NODE_DAM1_N(I,1),2)=I_NODE_DAM1_N(I,2)
     N_DAM_MATCH(I_NODE_DAM1_N(I,2),2)=I_NODE_DAM1_N(I,1)
   END DO
   DO I=1,NODE_DAM2_N
     IS_DAM(I_NODE_DAM2_N(I,1))=1
     IS_DAM(I_NODE_DAM2_N(I,2))=1
     IS_DAM(I_NODE_DAM2_N(I,3))=1
     N_DAM_MATCH(I_NODE_DAM2_N(I,1),1)=2
     N_DAM_MATCH(I_NODE_DAM2_N(I,2),1)=2
     N_DAM_MATCH(I_NODE_DAM2_N(I,3),1)=2

     N_DAM_MATCH(I_NODE_DAM2_N(I,1),2)=I_NODE_DAM2_N(I,2)
     N_DAM_MATCH(I_NODE_DAM2_N(I,1),3)=I_NODE_DAM2_N(I,3)

     N_DAM_MATCH(I_NODE_DAM2_N(I,2),2)=I_NODE_DAM2_N(I,3)
     N_DAM_MATCH(I_NODE_DAM2_N(I,2),3)=I_NODE_DAM2_N(I,1)

     N_DAM_MATCH(I_NODE_DAM2_N(I,3),2)=I_NODE_DAM2_N(I,1)
     N_DAM_MATCH(I_NODE_DAM2_N(I,3),3)=I_NODE_DAM2_N(I,2)
   END DO
   DO I=1,NODE_DAM3_N
     IS_DAM(I_NODE_DAM3_N(I,1))=1
     IS_DAM(I_NODE_DAM3_N(I,2))=1
     IS_DAM(I_NODE_DAM3_N(I,3))=1
     IS_DAM(I_NODE_DAM3_N(I,4))=1
     N_DAM_MATCH(I_NODE_DAM3_N(I,1),1)=3
     N_DAM_MATCH(I_NODE_DAM3_N(I,2),1)=3
     N_DAM_MATCH(I_NODE_DAM3_N(I,3),1)=3
     N_DAM_MATCH(I_NODE_DAM3_N(I,4),1)=3

     N_DAM_MATCH(I_NODE_DAM3_N(I,1),2)=I_NODE_DAM3_N(I,2)
     N_DAM_MATCH(I_NODE_DAM3_N(I,1),3)=I_NODE_DAM3_N(I,3)
     N_DAM_MATCH(I_NODE_DAM3_N(I,1),4)=I_NODE_DAM3_N(I,4)

     N_DAM_MATCH(I_NODE_DAM3_N(I,2),2)=I_NODE_DAM3_N(I,1)
     N_DAM_MATCH(I_NODE_DAM3_N(I,2),3)=I_NODE_DAM3_N(I,3)
     N_DAM_MATCH(I_NODE_DAM3_N(I,2),4)=I_NODE_DAM3_N(I,4)

     N_DAM_MATCH(I_NODE_DAM3_N(I,3),2)=I_NODE_DAM3_N(I,1)
     N_DAM_MATCH(I_NODE_DAM3_N(I,3),3)=I_NODE_DAM3_N(I,2)
     N_DAM_MATCH(I_NODE_DAM3_N(I,3),4)=I_NODE_DAM3_N(I,4)

     N_DAM_MATCH(I_NODE_DAM3_N(I,4),2)=I_NODE_DAM3_N(I,1)
     N_DAM_MATCH(I_NODE_DAM3_N(I,4),3)=I_NODE_DAM3_N(I,2)
     N_DAM_MATCH(I_NODE_DAM3_N(I,4),4)=I_NODE_DAM3_N(I,3)

   END DO

   IF(CELL_DAM_N>0)THEN
     DO I=1,NT
       CLP_COUNT=0
       ALPHA1=0.0
       DO NN=1,CELL_DAM_N
          ITMP1=I_CELL_DAM_N(NN,1)
          ITMP2=I_CELL_DAM_N(NN,2)
          IF(ITMP1==I.OR.ITMP2==I)GOTO 1111
       END DO
       COUNT = 0
       DO NN=1,NODE_DAM1_N
          ITMP1=I_NODE_DAM1_N(NN,1)
          ITMP2=I_NODE_DAM1_N(NN,2)
          IF(ITMP1==NV(I,1).OR.ITMP1==NV(I,2).OR.ITMP1==NV(I,3).OR. &
          &  ITMP2==NV(I,1).OR.ITMP2==NV(I,2).OR.ITMP2==NV(I,3) )THEN
             COUNT = COUNT + 1
          END IF
       END DO
       DO NN=1,NODE_DAM2_N
          ITMP1=I_NODE_DAM2_N(NN,1)
          ITMP2=I_NODE_DAM2_N(NN,2)
          ITMP3=I_NODE_DAM2_N(NN,3)
          IF(ITMP1==NV(I,1).OR.ITMP1==NV(I,2).OR.ITMP1==NV(I,3).OR. &
          &  ITMP2==NV(I,1).OR.ITMP2==NV(I,2).OR.ITMP2==NV(I,3).OR. &
          &  ITMP3==NV(I,1).OR.ITMP3==NV(I,2).OR.ITMP3==NV(I,3) )THEN
             COUNT = COUNT + 1
          END IF
       END DO
       DO NN=1,NODE_DAM3_N
          ITMP1=I_NODE_DAM3_N(NN,1)
          ITMP2=I_NODE_DAM3_N(NN,2)
          ITMP3=I_NODE_DAM3_N(NN,3)
          ITMP4=I_NODE_DAM3_N(NN,4)
          IF(ITMP1==NV(I,1).OR.ITMP1==NV(I,2).OR.ITMP1==NV(I,3).OR. &
          &  ITMP2==NV(I,1).OR.ITMP2==NV(I,2).OR.ITMP2==NV(I,3).OR. &
          &  ITMP3==NV(I,1).OR.ITMP3==NV(I,2).OR.ITMP3==NV(I,3).OR. &
          &  ITMP4==NV(I,1).OR.ITMP4==NV(I,2).OR.ITMP4==NV(I,3) )THEN
             COUNT = COUNT + 1
          END IF
       END DO
       IF(COUNT == 0) GOTO 1111

       DO NN=1,CELL_DAM_N
          ITMP=I_CELL_DAM_N(NN,1)
          IF(ITMP/=I)THEN
            IF(NBE(I,1)==ITMP.OR.NBE(I,2)==ITMP.OR.NBE(I,3)=&
                 &=ITMP)THEN
              CLP_COUNT=CLP_COUNT+1
              ALPHA1=ALPHA1+ALPHA(ITMP)
!              CLP_KDAM(I,CLP_COUNT+1)=ITMP
              CLP_KDAM(I,CLP_COUNT)=ITMP
            END IF
          END IF
          ITMP=I_CELL_DAM_N(NN,2)
          IF(ITMP/=I)THEN
            IF(NBE(I,1)==ITMP.OR.NBE(I,2)==ITMP.OR.NBE(I,3)=&
                 &=ITMP)THEN
              CLP_COUNT=CLP_COUNT+1
              ALPHA1=ALPHA1+ALPHA(ITMP)
!              CLP_KDAM(I,CLP_COUNT+1)=ITMP
              CLP_KDAM(I,CLP_COUNT)=ITMP
            END IF
          END IF
       END DO
1111   IF(CLP_COUNT>=1)THEN
         CLP_CELL(I)=1
         CLP_ALPHA(I) = ALPHA1/CLP_COUNT
         CLP_KDAM(I,1)  = CLP_COUNT 
       END IF
     END DO
   END IF   

# if !defined (SPHERICAL)
    DO I=1,NT
       DO I1=1,NODE_DAM1_N
          DTMP=(XC(I)-VX(I_NODE_DAM1_N(I1,1)))**2+(YC(I)&
               &-VY(I_NODE_DAM1_N(I1,1)))**2
          DTMP=SQRT(DTMP)/DAM1_R_SPG(I1)

          IF(DTMP <= 1.0_SP) THEN
             C_SPONGE=DAM1_C_SPG(I1)*(1.0_SP-DTMP)
             DAM_SPONGE(I)=MAX(C_SPONGE,DAM_SPONGE(I))
          END IF
       END DO

       DO I1=1,NODE_DAM2_N
          DTMP=(XC(I)-VX(I_NODE_DAM2_N(I1,1)))**2+(YC(I)&
               &-VY(I_NODE_DAM2_N(I1,1)))**2
          DTMP=SQRT(DTMP)/DAM2_R_SPG(I1)

          IF(DTMP <= 1.0_SP) THEN
             C_SPONGE=DAM2_C_SPG(I1)*(1.0_SP-DTMP)
             DAM_SPONGE(I)=MAX(C_SPONGE,DAM_SPONGE(I))
          END IF
       END DO

       DO I1=1,NODE_DAM3_N
          DTMP=(XC(I)-VX(I_NODE_DAM3_N(I1,1)))**2+(YC(I)&
               &-VY(I_NODE_DAM3_N(I1,1)))**2
          DTMP=SQRT(DTMP)/DAM3_R_SPG(I1)

          IF(DTMP <= 1.0_SP) THEN
             C_SPONGE=DAM3_C_SPG(I1)*(1.0_SP-DTMP)
             DAM_SPONGE(I)=MAX(C_SPONGE,DAM_SPONGE(I))
          END IF
       END DO
    END DO

# else
    ! SPHERICAL

    DO I=1,NT
       DO I1=1,NODE_DAM1_N
          X1_DP=XC(I)
          Y1_DP=YC(I)
          X2_DP=VX(I_NODE_DAM1_N(I1,1))
          Y2_DP=VY(I_NODE_DAM1_N(I1,1))
          CALL ARC(X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP)
          DTMP=DTMP_DP/DAM1_R_SPG(I1)

          IF(DTMP <= 1.0_SP) THEN
             C_SPONGE=DAM1_C_SPG(I1)*(1.0_SP-DTMP)
             DAM_SPONGE(I)=MAX(C_SPONGE,DAM_SPONGE(I))
          END IF
       END DO

       DO I1=1,NODE_DAM2_N
          X1_DP=XC(I)
          Y1_DP=YC(I)
          X2_DP=VX(I_NODE_DAM2_N(I1,1))
          Y2_DP=VY(I_NODE_DAM2_N(I1,1))
          CALL ARC(X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP)
          DTMP=DTMP_DP/DAM2_R_SPG(I1)

          IF(DTMP <= 1.0_SP) THEN
             C_SPONGE=DAM2_C_SPG(I1)*(1.0_SP-DTMP)
             DAM_SPONGE(I)=MAX(C_SPONGE,DAM_SPONGE(I))
          END IF
       END DO

       DO I1=1,NODE_DAM3_N
          X1_DP=XC(I)
          Y1_DP=YC(I)
          X2_DP=VX(I_NODE_DAM3_N(I1,1))
          Y2_DP=VY(I_NODE_DAM3_N(I1,1))
          CALL ARC(X1_DP,Y1_DP,X2_DP,Y2_DP,DTMP_DP)
          DTMP=DTMP_DP/DAM3_R_SPG(I1)

          IF(DTMP <= 1.0_SP) THEN
             C_SPONGE=DAM3_C_SPG(I1)*(1.0_SP-DTMP)
             DAM_SPONGE(I)=MAX(C_SPONGE,DAM_SPONGE(I))
          END IF
       END DO
    END DO

# endif    

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: read_dam" 

   RETURN
   END SUBROUTINE READ_DAM
!==============================================================================|     
!
!==============================================================================|
   SUBROUTINE ELE_DAM1
   
   IMPLICIT NONE
   REAL(SP) :: D1_TMP,D2_TMP,D_ELF,D_ELF1,D_ELF2,DH1,DH2
   INTEGER  :: I,I1,I2

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: ele_dam1 "  

   IF(SERIAL)THEN
     DO I=1,NN_DAM1
       I1 = NODE_DAM1(I,1)
       I2 = NODE_DAM1(I,2)

       D1_TMP = H(I1)+ELF(I1)
       D2_TMP = H(I2)+ELF(I2)

       IF(D1_TMP > D_DAM1(I,1) .AND. D2_TMP > D_DAM1(I,2))THEN
         IF(ABS(ELF(I1)-ELF(I2)) > EL_EPS)THEN
           D_ELF1 = D1_TMP-D_DAM1(I,1)
           D_ELF2 = D2_TMP-D_DAM1(I,2)
           D_ELF  = (D_ELF1*ART1(I1)+D_ELF2*ART1(I2))/(ART1(I1)+ART1(I2))
           ELF(I1) = D_DAM1(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM1(I,2) + D_ELF - H(I2)
         END IF
       ELSE IF(D1_TMP > D_DAM1(I,1) .AND. D2_TMP <= D_DAM1(I,2))THEN  
         D_ELF = D1_TMP-D_DAM1(I,1)
         ELF(I1) = ELF(I1)-D_ELF
         ELF(I2) = ELF(I2)+D_ELF*ART1(I1)/ART1(I2)
         D2_TMP = H(I2)+ELF(I2)
         IF(D2_TMP > D_DAM1(I,2))THEN
           D_ELF2 = D2_TMP-D_DAM1(I,2)
           D_ELF  = D_ELF2*ART1(I2)/(ART1(I1)+ART1(I2))
           ELF(I1) = D_DAM1(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM1(I,2) + D_ELF - H(I2)
         END IF  
       ELSE IF(D2_TMP > D_DAM1(I,2) .AND. D1_TMP <= D_DAM1(I,1))THEN 
         D_ELF = D2_TMP-D_DAM1(I,2)
         ELF(I1) = ELF(I1)+D_ELF*ART1(I2)/ART1(I1)
         ELF(I2) = ELF(I2)-D_ELF
         D1_TMP = H(I1)+ELF(I1)
         IF(D1_TMP > D_DAM1(I,1))THEN
           D_ELF1 = D1_TMP-D_DAM1(I,1)
           D_ELF  = D_ELF1*ART1(I1)/(ART1(I1)+ART1(I2))
           ELF(I1) = D_DAM1(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM1(I,2) + D_ELF - H(I2)
         END IF  
       END IF
     END DO
   END IF

#   if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I=1,NODE_DAM1_N
       I1  = I_NODE_DAM1_N(I,1)
       I2  = I_NODE_DAM1_N(I,2)
       DH1 = D_DAM1_N(I,1)
       DH2 = D_DAM1_N(I,2)

       D1_TMP = H(I1)+ELF(I1)
       D2_TMP = H(I2)+ELF(I2)
       
       IF(D1_TMP > DH1 .AND. D2_TMP > DH2)THEN
         IF(ABS(ELF(I1)-ELF(I2)) > EL_EPS)THEN
           D_ELF1 = D1_TMP-DH1
           D_ELF2 = D2_TMP-DH2
           D_ELF  = (D_ELF1*ART1(I1)+D_ELF2*ART1(I2))/(ART1(I1)+ART1(I2))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
         END IF
       ELSE IF(D1_TMP > DH1 .AND. D2_TMP <= DH2)THEN  
         D_ELF = D1_TMP-DH1
         ELF(I1) = ELF(I1)-D_ELF
         ELF(I2) = ELF(I2)+D_ELF*ART1(I1)/ART1(I2)
         D2_TMP = H(I2)+ELF(I2)
         IF(D2_TMP > DH2)THEN
           D_ELF2 = D2_TMP-DH2
           D_ELF  = D_ELF2*ART1(I2)/(ART1(I1)+ART1(I2))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
         END IF  
       ELSE IF(D2_TMP > DH2 .AND. D1_TMP <= DH1)THEN 
         D_ELF = D2_TMP-DH2
         ELF(I1) = ELF(I1)+D_ELF*ART1(I2)/ART1(I1)
         ELF(I2) = ELF(I2)-D_ELF
         D1_TMP = H(I1)+ELF(I1)
         IF(D1_TMP > DH1)THEN
           D_ELF1 = D1_TMP-DH1
           D_ELF  = D_ELF1*ART1(I1)/(ART1(I1)+ART1(I2))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
         END IF  
       END IF
     END DO
   END IF
#   endif

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: ele_dam1 "  

   RETURN
   END SUBROUTINE ELE_DAM1 
!==============================================================================|     
!
!==============================================================================|
   SUBROUTINE ELE_DAM2
   
   IMPLICIT NONE
   REAL(SP) :: D1_TMP,D2_TMP,D3_TMP,D_ELF,D_ELF1,D_ELF2,D_ELF3
   INTEGER  :: I,I1,I2,I3
   REAL(SP) :: DH1,DH2,DH3

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: ele_dam2 "  
   IF(SERIAL)THEN
     DO I=1,NN_DAM2
       I1 = NODE_DAM2(I,1)
       I2 = NODE_DAM2(I,2)
       I3 = NODE_DAM2(I,3)
       
       D1_TMP = H(I1)+ELF(I1)
       D2_TMP = H(I2)+ELF(I2)
       D3_TMP = H(I3)+ELF(I3)
     
       IF(D1_TMP > D_DAM2(I,1) .AND. D2_TMP > D_DAM2(I,2) .AND.   &
          D3_TMP > D_DAM2(I,3))THEN
         IF(ABS(ELF(I1)-ELF(I2)) > EL_EPS .OR.    &
            ABS(ELF(I2)-ELF(I3)) > EL_EPS .OR.    &
            ABS(ELF(I3)-ELF(I1)) > EL_EPS)THEN
           D_ELF1 = D1_TMP-D_DAM2(I,1)
           D_ELF2 = D2_TMP-D_DAM2(I,2)
           D_ELF3 = D3_TMP-D_DAM2(I,3)
           D_ELF  = (D_ELF1*ART1(I1)+D_ELF2*ART1(I2)+D_ELF3*ART1(I3))/   &
	          (ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = D_DAM2(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM2(I,2) + D_ELF - H(I2)
	   ELF(I3) = D_DAM2(I,3) + D_ELF - H(I3)
         END IF
       ELSE IF(D1_TMP > D_DAM2(I,1) .AND.                       &
             (D2_TMP <= D_DAM2(I,2) .AND. D3_TMP <= D_DAM2(I,3)))THEN  
         D_ELF1 = (ELF(I2)*ART1(I2)+ELF(I3)*ART1(I3))/(ART1(I2)+ART1(I3))
         ELF(I2) = D_ELF1
         ELF(I3) = D_ELF1
         
         D_ELF = D1_TMP-D_DAM2(I,1)
         ELF(I1) = ELF(I1)-D_ELF
       
         ELF(I2) = ELF(I2)+D_ELF*ART1(I1)/(ART1(I2)+ART1(I3))
         ELF(I3) = ELF(I2)
         D2_TMP = H(I2)+ELF(I2)
         IF(D2_TMP > D_DAM2(I,2))THEN
           D_ELF2 = D2_TMP-D_DAM2(I,2)
           D_ELF  = D_ELF2*(ART1(I2)+ART1(I3))/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = D_DAM2(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM2(I,2) + D_ELF - H(I2)
           ELF(I3) = D_DAM2(I,3) + D_ELF - H(I3)
         END IF  
       ELSE IF(D2_TMP > D_DAM2(I,2) .AND.                       &
              (D3_TMP <= D_DAM2(I,3) .AND. D1_TMP <= D_DAM2(I,1)))THEN  
         D_ELF1 = (ELF(I3)*ART1(I3)+ELF(I1)*ART1(I1))/(ART1(I3)+ART1(I1))
         ELF(I3) = D_ELF1
         ELF(I1) = D_ELF1
         
         D_ELF = D2_TMP-D_DAM2(I,2)
         ELF(I2) = ELF(I2)-D_ELF
       
         ELF(I3) = ELF(I3)+D_ELF*ART1(I2)/(ART1(I3)+ART1(I1))
         ELF(I1) = ELF(I3)
         D3_TMP = H(I3)+ELF(I3)
         IF(D3_TMP > D_DAM2(I,3))THEN
           D_ELF3 = D3_TMP-D_DAM2(I,3)
           D_ELF  = D_ELF3*(ART1(I3)+ART1(I1))/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = D_DAM2(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM2(I,2) + D_ELF - H(I2)
           ELF(I3) = D_DAM2(I,3) + D_ELF - H(I3)
         END IF  
       ELSE IF(D3_TMP > D_DAM2(I,3) .AND.                       &
              (D1_TMP <= D_DAM2(I,1) .AND. D2_TMP <= D_DAM2(I,2)))THEN  
         D_ELF1 = (ELF(I1)*ART1(I1)+ELF(I2)*ART1(I2))/(ART1(I1)+ART1(I2))
         ELF(I1) = D_ELF1
         ELF(I2) = D_ELF1
         
         D_ELF = D3_TMP-D_DAM2(I,3)
         ELF(I3) = ELF(I3)-D_ELF
       
         ELF(I1) = ELF(I1)+D_ELF*ART1(I3)/(ART1(I1)+ART1(I2))
         ELF(I2) = ELF(I1)
         D1_TMP = H(I1)+ELF(I1)
         IF(D1_TMP > D_DAM2(I,1))THEN
           D_ELF1 = D1_TMP-D_DAM2(I,1)
           D_ELF  = D_ELF1*(ART1(I1)+ART1(I2))/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = D_DAM2(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM2(I,2) + D_ELF - H(I2)
           ELF(I3) = D_DAM2(I,3) + D_ELF - H(I3)
         END IF  
       ELSE IF(D1_TMP > D_DAM2(I,1) .AND. D2_TMP > D_DAM2(I,2) .AND.  &
               D3_TMP <= D_DAM2(I,3))THEN
         D_ELF1 = D1_TMP - D_DAM2(I,1)
         D_ELF2 = D2_TMP - D_DAM2(I,2)
         D_ELF = (D_ELF1*ART1(I1)+D_ELF2*ART1(I2))/(ART1(I1)+ART1(I2))
         ELF(I1) = D_DAM2(I,1) + D_ELF -H(I1)      
         ELF(I2) = D_DAM2(I,2) + D_ELF -H(I2)
       
         D1_TMP = H(I1)+ELF(I1)
         D_ELF = D1_TMP - D_DAM2(I,1)
         ELF(I3) = ELF(I3)+D_ELF*(ART1(I1)+ART1(I2))/ART1(I3)
       
         D3_TMP = H(I3)+ELF(I3)
         IF(D3_TMP > D_DAM2(I,3))THEN
           D_ELF3 = D3_TMP-D_DAM2(I,3)
	   D_ELF = D_ELF3*ART1(I3)/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = D_DAM2(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM2(I,2) + D_ELF - H(I2)
           ELF(I3) = D_DAM2(I,3) + D_ELF - H(I3)
         END IF
       ELSE IF(D2_TMP > D_DAM2(I,2) .AND. D3_TMP > D_DAM2(I,3) .AND.  &
               D1_TMP <= D_DAM2(I,1))THEN
         D_ELF2 = D2_TMP - D_DAM2(I,2)
         D_ELF3 = D3_TMP - D_DAM2(I,3)
         D_ELF = (D_ELF2*ART1(I2)+D_ELF3*ART1(I3))/(ART1(I2)+ART1(I3))
         ELF(I2) = D_DAM2(I,2) + D_ELF -H(I2)      
         ELF(I3) = D_DAM2(I,3) + D_ELF -H(I3)
       
         D2_TMP = H(I2)+ELF(I2)
         D_ELF = D2_TMP - D_DAM2(I,2)
         ELF(I1) = ELF(I1)+D_ELF*(ART1(I2)+ART1(I3))/ART1(I1)
       
         D1_TMP = H(I1)+ELF(I1)
         IF(D1_TMP > D_DAM2(I,1))THEN
           D_ELF1 = D1_TMP-D_DAM2(I,1)
	   D_ELF = D_ELF1*ART1(I1)/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = D_DAM2(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM2(I,2) + D_ELF - H(I2)
           ELF(I3) = D_DAM2(I,3) + D_ELF - H(I3)
         END IF
       ELSE IF(D3_TMP > D_DAM2(I,3) .AND. D1_TMP > D_DAM2(I,1) .AND.  &
               D2_TMP <= D_DAM2(I,2))THEN
         D_ELF3 = D3_TMP - D_DAM2(I,3)
         D_ELF1 = D1_TMP - D_DAM2(I,1)
         D_ELF = (D_ELF3*ART1(I3)+D_ELF1*ART1(I1))/(ART1(I3)+ART1(I1))
         ELF(I3) = D_DAM2(I,3) + D_ELF -H(I3)      
         ELF(I1) = D_DAM2(I,1) + D_ELF -H(I1)
       
         D3_TMP = H(I3)+ELF(I3)
         D_ELF = D3_TMP - D_DAM2(I,3)
         ELF(I2) = ELF(I2)+D_ELF*(ART1(I3)+ART1(I1))/ART1(I2)
       
         D2_TMP = H(I2)+ELF(I2)
         IF(D2_TMP > D_DAM2(I,2))THEN
           D_ELF2 = D2_TMP-D_DAM2(I,2)
	   D_ELF = D_ELF2*ART1(I2)/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = D_DAM2(I,1) + D_ELF - H(I1)
           ELF(I2) = D_DAM2(I,2) + D_ELF - H(I2)
           ELF(I3) = D_DAM2(I,3) + D_ELF - H(I3)
         END IF
       END IF 
     END DO
   END IF

#  if defined (MULTIPROCESSOR)
     DO I=1,NODE_DAM2_N
       I1 = I_NODE_DAM2_N(I,1)
       I2 = I_NODE_DAM2_N(I,2)
       I3 = I_NODE_DAM2_N(I,3)
       DH1 = D_DAM2_N(I,1)       
       DH2 = D_DAM2_N(I,2)
       DH3 = D_DAM2_N(I,3) 
 
       D1_TMP = H(I1)+ELF(I1)
       D2_TMP = H(I2)+ELF(I2)
       D3_TMP = H(I3)+ELF(I3)
     
       IF(D1_TMP > DH1 .AND. D2_TMP > DH2 .AND.   &
          D3_TMP > DH3)THEN
         IF(ABS(ELF(I1)-ELF(I2)) > EL_EPS .OR.    &
            ABS(ELF(I2)-ELF(I3)) > EL_EPS .OR.    &
            ABS(ELF(I3)-ELF(I1)) > EL_EPS)THEN
           D_ELF1 = D1_TMP-DH1
           D_ELF2 = D2_TMP-DH2
           D_ELF3 = D3_TMP-DH3
           D_ELF  = (D_ELF1*ART1(I1)+D_ELF2*ART1(I2)+D_ELF3*ART1(I3))/   &
	          (ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
	   ELF(I3) = DH3 + D_ELF - H(I3)
         END IF
       ELSE IF(D1_TMP > DH1 .AND.                       &
             (D2_TMP <= DH2 .AND. D3_TMP <= DH3))THEN  
         D_ELF1 = (ELF(I2)*ART1(I2)+ELF(I3)*ART1(I3))/(ART1(I2)+ART1(I3))
         ELF(I2) = D_ELF1
         ELF(I3) = D_ELF1
         
         D_ELF = D1_TMP-DH1
         ELF(I1) = ELF(I1)-D_ELF
       
         ELF(I2) = ELF(I2)+D_ELF*ART1(I1)/(ART1(I2)+ART1(I3))
         ELF(I3) = ELF(I2)
         D2_TMP = H(I2)+ELF(I2)
         IF(D2_TMP > DH2)THEN
           D_ELF2 = D2_TMP-DH2
           D_ELF  = D_ELF2*(ART1(I2)+ART1(I3))/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
         END IF  
       ELSE IF(D2_TMP > DH2 .AND.                       &
              (D3_TMP <= DH3 .AND. D1_TMP <= DH1))THEN  
         D_ELF1 = (ELF(I3)*ART1(I3)+ELF(I1)*ART1(I1))/(ART1(I3)+ART1(I1))
         ELF(I3) = D_ELF1
         ELF(I1) = D_ELF1
         
         D_ELF = D2_TMP-DH2
         ELF(I2) = ELF(I2)-D_ELF
       
         ELF(I3) = ELF(I3)+D_ELF*ART1(I2)/(ART1(I3)+ART1(I1))
         ELF(I1) = ELF(I3)
         D3_TMP = H(I3)+ELF(I3)
         IF(D3_TMP > DH3)THEN
           D_ELF3 = D3_TMP-DH3
           D_ELF  = D_ELF3*(ART1(I3)+ART1(I1))/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
         END IF  
       ELSE IF(D3_TMP > DH3 .AND.                       &
              (D1_TMP <= DH1 .AND. D2_TMP <= DH2))THEN  
         D_ELF1 = (ELF(I1)*ART1(I1)+ELF(I2)*ART1(I2))/(ART1(I1)+ART1(I2))
         ELF(I1) = D_ELF1
         ELF(I2) = D_ELF1
         
         D_ELF = D3_TMP-DH3
         ELF(I3) = ELF(I3)-D_ELF
       
         ELF(I1) = ELF(I1)+D_ELF*ART1(I3)/(ART1(I1)+ART1(I2))
         ELF(I2) = ELF(I1)
         D1_TMP = H(I1)+ELF(I1)
         IF(D1_TMP > DH1)THEN
           D_ELF1 = D1_TMP-DH1
           D_ELF  = D_ELF1*(ART1(I1)+ART1(I2))/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
         END IF  
       ELSE IF(D1_TMP > DH1 .AND. D2_TMP > DH2 .AND.  &
               D3_TMP <= DH3)THEN
         D_ELF1 = D1_TMP - DH1
         D_ELF2 = D2_TMP - DH2
         D_ELF = (D_ELF1*ART1(I1)+D_ELF2*ART1(I2))/(ART1(I1)+ART1(I2))
         ELF(I1) = DH1 + D_ELF -H(I1)      
         ELF(I2) = DH2 + D_ELF -H(I2)
       
         D1_TMP = H(I1)+ELF(I1)
         D_ELF = D1_TMP - DH1
         ELF(I3) = ELF(I3)+D_ELF*(ART1(I1)+ART1(I2))/ART1(I3)
       
         D3_TMP = H(I3)+ELF(I3)
         IF(D3_TMP > DH3)THEN
           D_ELF3 = D3_TMP-DH3
	   D_ELF = D_ELF3*ART1(I3)/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
         END IF
       ELSE IF(D2_TMP > DH2 .AND. D3_TMP > DH3 .AND.  &
               D1_TMP <= DH1)THEN
         D_ELF2 = D2_TMP - DH2
         D_ELF3 = D3_TMP - DH3
         D_ELF = (D_ELF2*ART1(I2)+D_ELF3*ART1(I3))/(ART1(I2)+ART1(I3))
         ELF(I2) = DH2 + D_ELF -H(I2)      
         ELF(I3) = DH3 + D_ELF -H(I3)
       
         D2_TMP = H(I2)+ELF(I2)
         D_ELF = D2_TMP - DH2
         ELF(I1) = ELF(I1)+D_ELF*(ART1(I2)+ART1(I3))/ART1(I1)
       
         D1_TMP = H(I1)+ELF(I1)
         IF(D1_TMP > DH1)THEN
           D_ELF1 = D1_TMP-DH1
	   D_ELF = D_ELF1*ART1(I1)/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
         END IF
       ELSE IF(D3_TMP > DH3 .AND. D1_TMP > DH1 .AND.  &
               D2_TMP <= DH2)THEN
         D_ELF3 = D3_TMP - DH3
         D_ELF1 = D1_TMP - DH1
         D_ELF = (D_ELF3*ART1(I3)+D_ELF1*ART1(I1))/(ART1(I3)+ART1(I1))
         ELF(I3) = DH3 + D_ELF -H(I3)      
         ELF(I1) = DH1 + D_ELF -H(I1)
       
         D3_TMP = H(I3)+ELF(I3)
         D_ELF = D3_TMP - DH3
         ELF(I2) = ELF(I2)+D_ELF*(ART1(I3)+ART1(I1))/ART1(I2)
       
         D2_TMP = H(I2)+ELF(I2)
         IF(D2_TMP > DH2)THEN
           D_ELF2 = D2_TMP-DH2
	   D_ELF = D_ELF2*ART1(I2)/(ART1(I1)+ART1(I2)+ART1(I3))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
         END IF
       END IF 
     END DO
#  endif

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: ele_dam2 "  
  
   RETURN
   END SUBROUTINE ELE_DAM2 
!==============================================================================|     
!
!==============================================================================|
   SUBROUTINE ELE_DAM3
   
   IMPLICIT NONE
   REAL(SP) :: D1_TMP,D2_TMP,D3_TMP,D4_TMP
   REAL(SP) :: D_ELF,D_ELF1,D_ELF2,D_ELF3,D_ELF4
   REAL(SP) :: DH1,DH2,DH3,DH4
   INTEGER  :: I,I1,I2,I3,I4
 
   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: ele_dam3 " 
     DO I=1,NODE_DAM3_N
       I1 = I_NODE_DAM3_N(I,1)
       I2 = I_NODE_DAM3_N(I,2)
       I3 = I_NODE_DAM3_N(I,3)
       I4 = I_NODE_DAM3_N(I,4)
   
       DH1 = D_DAM3_N(I,1)     
       DH2 = D_DAM3_N(I,2) 
       DH3 = D_DAM3_N(I,3) 
       DH4 = D_DAM3_N(I,4) 

       D1_TMP = H(I1)+ELF(I1)
       D2_TMP = H(I2)+ELF(I2)
       D3_TMP = H(I3)+ELF(I3)
       D4_TMP = H(I4)+ELF(I4)
     
       IF(D1_TMP > DH1 .AND. D2_TMP > DH2 .AND.   &
          D3_TMP > DH3 .AND. D4_TMP > DH4)THEN
         IF(ABS(ELF(I1)-ELF(I2)) > EL_EPS .OR.    &
            ABS(ELF(I2)-ELF(I3)) > EL_EPS .OR.    &
	    ABS(ELF(I3)-ELF(I4)) > EL_EPS .OR.    &
	    ABS(ELF(I4)-ELF(I1)) > EL_EPS)THEN
           D_ELF1 = D1_TMP-DH1
           D_ELF2 = D2_TMP-DH2
	   D_ELF3 = D3_TMP-DH3
	   D_ELF4 = D4_TMP-DH4
           D_ELF  = (D_ELF1*ART1(I1)+D_ELF2*ART1(I2)+D_ELF3*ART1(I3)+   &
	             D_ELF4*ART1(I4))/                                  &
	            (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
	   ELF(I3) = DH3 + D_ELF - H(I3)
	   ELF(I4) = DH4 + D_ELF - H(I4)
         END IF
       ELSE IF(D1_TMP >  DH1 .AND. D2_TMP <= DH2 .AND.   &
               D3_TMP <= DH3 .AND. D4_TMP <= DH4)THEN  
         D_ELF1 = (ELF(I2)*ART1(I2)+ELF(I3)*ART1(I3)+ELF(I4)*ART1(I4))/  &
                    (ART1(I2)+ART1(I3)+ART1(I4))
         ELF(I2) = D_ELF1
         ELF(I3) = D_ELF1
         ELF(I4) = D_ELF1
         
         D_ELF = D1_TMP-DH1
         ELF(I1) = ELF(I1)-D_ELF
       
         ELF(I2) = ELF(I2)+D_ELF*ART1(I1)/(ART1(I2)+ART1(I3)+ART1(I4))
         ELF(I3) = ELF(I2)
         ELF(I4) = ELF(I2)
         D2_TMP = H(I2)+ELF(I2)
         IF(D2_TMP > DH2)THEN
           D_ELF2 = D2_TMP-DH2
           D_ELF  = D_ELF2*(ART1(I2)+ART1(I3)+ART1(I4))/               &
	            (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF  
       ELSE IF(D2_TMP >  DH2 .AND. D3_TMP <= DH3 .AND.   &
               D4_TMP <= DH4 .AND. D1_TMP <= DH1)THEN  
         D_ELF1 = (ELF(I3)*ART1(I3)+ELF(I4)*ART1(I4)+ELF(I1)*ART1(I1))/  &
                  (ART1(I3)+ART1(I4)+ART1(I1))
         ELF(I3) = D_ELF1
         ELF(I4) = D_ELF1
         ELF(I1) = D_ELF1
         
         D_ELF = D2_TMP-DH2
         ELF(I2) = ELF(I2)-D_ELF
       
         ELF(I3) = ELF(I3)+D_ELF*ART1(I2)/(ART1(I3)+ART1(I4)+ART1(I1))
         ELF(I4) = ELF(I3)
         ELF(I1) = ELF(I3)
         D3_TMP = H(I3)+ELF(I3)
         IF(D3_TMP > DH3)THEN
           D_ELF3 = D3_TMP-DH3
           D_ELF  = D_ELF3*(ART1(I3)+ART1(I4)+ART1(I1))/               &
	            (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF  
       ELSE IF(D3_TMP >  DH3 .AND. D4_TMP <= DH4 .AND.   &
               D1_TMP <= DH1 .AND. D2_TMP <= DH2)THEN  
         D_ELF1 = (ELF(I4)*ART1(I4)+ELF(I1)*ART1(I1)+ELF(I2)*ART1(I2))/  &
                  (ART1(I4)+ART1(I1)+ART1(I2))
         ELF(I4) = D_ELF1
         ELF(I1) = D_ELF1
         ELF(I2) = D_ELF1
         
         D_ELF = D3_TMP-DH3
         ELF(I3) = ELF(I3)-D_ELF
       
         ELF(I4) = ELF(I4)+D_ELF*ART1(I3)/(ART1(I4)+ART1(I1)+ART1(I2))
         ELF(I1) = ELF(I4)
         ELF(I2) = ELF(I4)
         D4_TMP = H(I4)+ELF(I4)
         IF(D4_TMP > DH4)THEN
           D_ELF4 = D4_TMP-DH4
           D_ELF  = D_ELF4*(ART1(I4)+ART1(I1)+ART1(I2))/               &
	   (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF  
       ELSE IF(D4_TMP >  DH4 .AND. D1_TMP <= DH1 .AND.   &
               D2_TMP <= DH2 .AND. D3_TMP <= DH3)THEN  
         D_ELF1 = (ELF(I1)*ART1(I1)+ELF(I2)*ART1(I2)+ELF(I3)*ART1(I3))/  &
                  (ART1(I1)+ART1(I2)+ART1(I3))
         ELF(I1) = D_ELF1
         ELF(I2) = D_ELF1
         ELF(I3) = D_ELF1
         
         D_ELF = D4_TMP-DH4
         ELF(I4) = ELF(I4)-D_ELF
       
         ELF(I1) = ELF(I1)+D_ELF*ART1(I4)/(ART1(I1)+ART1(I2)+ART1(I3))
         ELF(I2) = ELF(I1)
         ELF(I3) = ELF(I1)
         D1_TMP = H(I1)+ELF(I1)
         IF(D1_TMP > DH1)THEN
           D_ELF1 = D1_TMP-DH1
           D_ELF  = D_ELF1*(ART1(I1)+ART1(I2)+ART1(I3))/               &
	   (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF  
       ELSE IF(D1_TMP >  DH1 .AND. D2_TMP > DH2 .AND.  &
               D3_TMP <= DH3 .AND. D4_TMP <= DH4)THEN
         D_ELF1 = D1_TMP - DH1
         D_ELF2 = D2_TMP - DH2
         D_ELF = (D_ELF1*ART1(I1)+D_ELF2*ART1(I2))/(ART1(I1)+ART1(I2))
         ELF(I1) = DH1 + D_ELF -H(I1)      
         ELF(I2) = DH2 + D_ELF -H(I2)
       
         D_ELF = (ELF(I3)*ART1(I3)+ELF(I4)*ART1(I4))/(ART1(I3)+ART1(I4))
         ELF(I3) = D_ELF
         ELF(I4) = D_ELF
       
         D1_TMP = H(I1)+ELF(I1)
         D_ELF = D1_TMP - DH1
         ELF(I3) = ELF(I3)+D_ELF*(ART1(I1)+ART1(I2))/(ART1(I3)+ART1(I4))
         ELF(I4) = ELF(I4)
       
         D3_TMP = H(I3)+ELF(I3)
         IF(D3_TMP > DH3)THEN
           D_ELF3 = D3_TMP-DH3
	   D_ELF = D_ELF3*(ART1(I3)+ART1(I4))/                   &
	           (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF
       ELSE IF(D2_TMP >  DH2 .AND. D3_TMP >  DH3 .AND.  &
               D4_TMP <= DH4 .AND. D1_TMP <= DH1)THEN
         D_ELF2 = D2_TMP - DH2
         D_ELF3 = D3_TMP - DH3
         D_ELF = (D_ELF2*ART1(I2)+D_ELF3*ART1(I3))/(ART1(I2)+ART1(I3))
         ELF(I2) = DH2 + D_ELF -H(I2)      
         ELF(I3) = DH3 + D_ELF -H(I3)
       
         D_ELF = (ELF(I4)*ART1(I4)+ELF(I1)*ART1(I1))/(ART1(I4)+ART1(I1))
         ELF(I4) = D_ELF
         ELF(I1) = D_ELF
       
         D2_TMP = H(I2)+ELF(I2)
         D_ELF = D2_TMP - DH2
         ELF(I4) = ELF(I4)+D_ELF*(ART1(I2)+ART1(I3))/(ART1(I4)+ART1(I1))
         ELF(I1) = ELF(I4)
       
         D4_TMP = H(I4)+ELF(I4)
         IF(D4_TMP > DH4)THEN
           D_ELF4 = D4_TMP-DH4
	   D_ELF = D_ELF4*(ART1(I4)+ART1(I1))/                 &
	           (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF
       ELSE IF(D3_TMP >  DH3 .AND. D4_TMP >  DH4 .AND.  &
               D1_TMP <= DH1 .AND. D2_TMP <= DH2)THEN
         D_ELF3 = D3_TMP - DH3
         D_ELF4 = D4_TMP - DH4
         D_ELF = (D_ELF3*ART1(I3)+D_ELF4*ART1(I4))/(ART1(I3)+ART1(I4))
         ELF(I3) = DH3 + D_ELF -H(I3)      
         ELF(I4) = DH4 + D_ELF -H(I4)
       
         D_ELF = (ELF(I1)*ART1(I1)+ELF(I2)*ART1(I2))/(ART1(I1)+ART1(I2))
         ELF(I1) = D_ELF
         ELF(I2) = D_ELF
       
         D3_TMP = H(I3)+ELF(I3)
         D_ELF = D3_TMP - DH3
         ELF(I1) = ELF(I1)+D_ELF*(ART1(I3)+ART1(I4))/(ART1(I1)+ART1(I2))
         ELF(I2) = ELF(I1)
       
         D1_TMP = H(I1)+ELF(I1)
         IF(D1_TMP > DH1)THEN
           D_ELF1 = D1_TMP-DH1
	   D_ELF = D_ELF1*(ART1(I1)+ART1(I2))/                  &
	           (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF
       ELSE IF(D4_TMP >  DH4 .AND. D1_TMP >  DH1 .AND.  &
               D2_TMP <= DH2 .AND. D3_TMP <= DH3)THEN
         D_ELF4 = D4_TMP - DH4
         D_ELF1 = D1_TMP - DH1
         D_ELF = (D_ELF4*ART1(I4)+D_ELF1*ART1(I1))/(ART1(I4)+ART1(I1))
         ELF(I4) = DH4 + D_ELF -H(I4)      
         ELF(I1) = DH1 + D_ELF -H(I1)
       
         D_ELF = (ELF(I2)*ART1(I2)+ELF(I3)*ART1(I3))/(ART1(I2)+ART1(I3))
         ELF(I2) = D_ELF
         ELF(I3) = D_ELF
       
         D4_TMP = H(I4)+ELF(I4)
         D_ELF = D4_TMP - DH4
         ELF(I2) = ELF(I2)+D_ELF*(ART1(I4)+ART1(I1))/(ART1(I2)+ART1(I3))
         ELF(I3) = ELF(I2)
       
         D2_TMP = H(I2)+ELF(I2)
         IF(D2_TMP > DH2)THEN
           D_ELF2 = D2_TMP-DH2
           D_ELF = D_ELF2*(ART1(I2)+ART1(I3))/                  &
    	         (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF
       ELSE IF(D1_TMP > DH1 .AND. D2_TMP >  DH2 .AND.  &
               D3_TMP > DH3 .AND. D4_TMP <= DH4)THEN
         D_ELF1 = D1_TMP - DH1
         D_ELF2 = D2_TMP - DH2
         D_ELF3 = D3_TMP - DH3
         D_ELF = (D_ELF1*ART1(I1)+D_ELF2*ART1(I2)+D_ELF3*ART1(I3))/    &
                 (ART1(I1)+ART1(I2)+ART1(I3))
         ELF(I1) = DH1 + D_ELF -H(I1)      
         ELF(I2) = DH2 + D_ELF -H(I2)
         ELF(I3) = DH3 + D_ELF -H(I3)
       
         D1_TMP = H(I1)+ELF(I1)
         D_ELF = D1_TMP - DH1
         ELF(I4) = ELF(I4)+D_ELF*(ART1(I1)+ART1(I2)+ART1(I3))/ART1(I4)
       
         D4_TMP = H(I4)+ELF(I4)
         IF(D4_TMP > DH4)THEN
            D_ELF4 = D4_TMP-DH4
            D_ELF = D_ELF4*ART1(I4)/(ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF
       ELSE IF(D2_TMP > DH2 .AND. D3_TMP >  DH3 .AND.  &
               D4_TMP > DH4 .AND. D1_TMP <= DH1)THEN
         D_ELF2 = D2_TMP - DH2
         D_ELF3 = D3_TMP - DH3
         D_ELF4 = D4_TMP - DH4
         D_ELF = (D_ELF2*ART1(I2)+D_ELF3*ART1(I3)+D_ELF4*ART1(I4))/    &
                 (ART1(I2)+ART1(I3)+ART1(I4))
         ELF(I2) = DH2 + D_ELF -H(I2)      
         ELF(I3) = DH3 + D_ELF -H(I3)
         ELF(I4) = DH4 + D_ELF -H(I4)
       
         D2_TMP = H(I2)+ELF(I2)
         D_ELF = D2_TMP - DH2
         ELF(I1) = ELF(I1)+D_ELF*(ART1(I2)+ART1(I3)+ART1(I4))/ART1(I1)
       
         D1_TMP = H(I1)+ELF(I1)
         IF(D1_TMP > DH1)THEN
           D_ELF1 = D1_TMP-DH1
           D_ELF = D_ELF1*ART1(I1)/(ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF
       ELSE IF(D3_TMP > DH3 .AND. D4_TMP >  DH4 .AND.  &
               D1_TMP > DH1 .AND. D2_TMP <= DH2)THEN
         D_ELF3 = D3_TMP - DH3
         D_ELF4 = D4_TMP - DH4
         D_ELF1 = D1_TMP - DH1
         D_ELF = (D_ELF3*ART1(I3)+D_ELF4*ART1(I4)+D_ELF1*ART1(I1))/    &
                 (ART1(I3)+ART1(I4)+ART1(I1))
         ELF(I3) = DH3 + D_ELF -H(I3)      
         ELF(I4) = DH4 + D_ELF -H(I4)
         ELF(I1) = DH1 + D_ELF -H(I1)
       
         D3_TMP = H(I3)+ELF(I3)
         D_ELF = D3_TMP - DH3
         ELF(I2) = ELF(I2)+D_ELF*(ART1(I3)+ART1(I4)+ART1(I1))/ART1(I2)
       
         D2_TMP = H(I2)+ELF(I2)
         IF(D2_TMP > DH2)THEN
           D_ELF2 = D2_TMP-DH2
    	   D_ELF = D_ELF2*ART1(I2)/(ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF
       ELSE IF(D4_TMP > DH4 .AND. D1_TMP >  DH1 .AND.  &
               D2_TMP > DH2 .AND. D3_TMP <= DH3)THEN
         D_ELF4 = D4_TMP - DH4
         D_ELF1 = D1_TMP - DH1
         D_ELF2 = D2_TMP - DH2
         D_ELF = (D_ELF4*ART1(I4)+D_ELF1*ART1(I1)+D_ELF2*ART1(I2))/    &
                 (ART1(I4)+ART1(I1)+ART1(I2))
         ELF(I4) = DH4 + D_ELF -H(I4)      
         ELF(I1) = DH1 + D_ELF -H(I1)
         ELF(I2) = DH2 + D_ELF -H(I2)
       
         D4_TMP = H(I4)+ELF(I4)
         D_ELF = D4_TMP - DH4
         ELF(I3) = ELF(I3)+D_ELF*(ART1(I4)+ART1(I1)+ART1(I2))/ART1(I3)
       
         D3_TMP = H(I3)+ELF(I3)
         IF(D3_TMP > DH3)THEN
           D_ELF3 = D3_TMP-DH3
           D_ELF = D_ELF3*ART1(I3)/(ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))
           ELF(I1) = DH1 + D_ELF - H(I1)
           ELF(I2) = DH2 + D_ELF - H(I2)
           ELF(I3) = DH3 + D_ELF - H(I3)
           ELF(I4) = DH4 + D_ELF - H(I4)
         END IF
       END IF 
     END DO


   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: ele_dam3 "    

   RETURN
   END SUBROUTINE ELE_DAM3 
!==============================================================================|     
!
!==============================================================================|
   SUBROUTINE SHAPE_COEF_DAM

!----------------------------------------------------------------------!
!  This subroutine is used to calculate the coefficient for a linear   !
!  function on the x-y plane, i.e.:                                    !
!                                                                      !
!                     r(x,y;phai)=phai_c+cofa1*x+cofa2*y               !
!                                                                      !
!  This subroutine is used for dam cell boundary condition cases       !
!----------------------------------------------------------------------!

   USE ALL_VARS
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif
   IMPLICIT NONE
   REAL(DP) X1,X2,X3,Y1,Y2,Y3,DELT,AI1,AI2,AI3,BI1,BI2,BI3,CI1,CI2,CI3
   REAL(DP) DELTX,DELTY
   INTEGER  I,II,J,JJ,J1,J2,I1,JJ1,JJ2
# if defined (SPHERICAL)
   REAL(DP) XXC1,YYC1,XXC2,YYC2,XXC3,YYC3,SIDE,X1_DP,Y1_DP,X2_DP,Y2_DP
# endif

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: shape_ceof_dam " 

!
!---------------dam boundary cells------------------------------------!
!

   DO II = 1, CELL_DAM_N

!-------one side of dam ----------------------------------
   JJ1 = I_CELL_DAM_N(II,1)
   JJ2 = I_CELL_DAM_N(II,2)
   IF(JJ1 /= 0 .AND.JJ2 /= 0 )THEN
     I  = I_CELL_DAM_N(II,1)
     I1 = I_CELL_DAM_N(II,2)

     IF(ISBCE(I) == 1) THEN
       DO J = 1, 3
         IF(NBE(I,J) == 0) JJ = J
       END DO
       J1 = JJ+1-INT((JJ+1)/4)*3
       J2 = JJ+2-INT((JJ+2)/4)*3

       DELTX = VX(NV(I,J1))-VX(NV(I,J2))
# if defined (SPHERICAL)
       IF(DELTX > 180.0_SP)THEN
         DELTX = -360.0_SP+DELTX
       ELSE IF(DELTX < -180.0_SP)THEN
         DELTX =  360.0_SP+DELTX	 
       END IF	 
# endif       
       DELTY = VY(NV(I,J1))-VY(NV(I,J2))

       IF(JJ == 1)THEN
#        if defined (SPHERICAL)
         Y1 = YC(I1)-YC(I)
         Y2 = YC(NBE(I,J1))-YC(I)
         Y3 = YC(NBE(I,J2))-YC(I)
         Y1 = TPI*Y1
         Y2 = TPI*Y2
         Y3 = TPI*Y3

         X1_DP = XC(I)
         Y1_DP = YC(I)

         X2_DP = XC(I1)
         Y2_DP = YC(I1)
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X1 = SIDE

         X2_DP = XC(NBE(I,J1))
         Y2_DP = YC(NBE(I,J1))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X2 = SIDE

         X2_DP = XC(NBE(I,J2))
         Y2_DP = YC(NBE(I,J2))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X3 = SIDE
#        else
         X1 = XC(I1)-XC(I)
         Y1 = YC(I1)-YC(I)
         X2 = XC(NBE(I,J1))-XC(I)
         Y2 = YC(NBE(I,J1))-YC(I)
         X3 = XC(NBE(I,J2))-XC(I)
         Y3 = YC(NBE(I,J2))-YC(I)
#        endif
       ELSE IF(JJ == 2)THEN
#        if defined (SPHERICAL)
         Y1 = YC(NBE(I,J2))-YC(I)
         Y2 = YC(I1)-YC(I)
         Y3 = YC(NBE(I,J1))-YC(I)
         Y1 = TPI*Y1
         Y2 = TPI*Y2
         Y3 = TPI*Y3

         X1_DP = XC(I)
         Y1_DP = YC(I)

         X2_DP = XC(NBE(I,J2))
         Y2_DP = YC(NBE(I,J2))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X1 = SIDE

         X2_DP = XC(I1)
         Y2_DP = YC(I1)
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X2 = SIDE

         X2_DP = XC(NBE(I,J1))
         Y2_DP = YC(NBE(I,J1))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X3 = SIDE
#        else
         X1 = XC(NBE(I,J2))-XC(I)
         Y1 = YC(NBE(I,J2))-YC(I)
         X2 = XC(I1)-XC(I)
         Y2 = YC(I1)-YC(I)
         X3 = XC(NBE(I,J1))-XC(I)
         Y3 = YC(NBE(I,J1))-YC(I)
#        endif
       ELSE
#        if defined (SPHERICAL)
         Y1 = YC(NBE(I,J1))-YC(I)
         Y2 = YC(NBE(I,J2))-YC(I)
         Y3 = YC(I1)-YC(I)
         Y1 = TPI*Y1
         Y2 = TPI*Y2
         Y3 = TPI*Y3

         X1_DP = XC(I)
         Y1_DP = YC(I)
         X2_DP = XC(NBE(I,J1))
         Y2_DP = YC(NBE(I,J1))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X1 = SIDE

         X2_DP = XC(NBE(I,J2))
         Y2_DP = YC(NBE(I,J2))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X2 = SIDE

         X2_DP = XC(I1)
         Y2_DP = YC(I1)
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X3 = SIDE
#        else
         X1 = XC(NBE(I,J1))-XC(I)
         y1 = YC(NBE(I,J1))-YC(I)
         X2 = XC(NBE(I,J2))-XC(I)
         Y2 = YC(NBE(I,J2))-YC(I)
         X3 = XC(I1)-XC(I)
         Y3 = YC(I1)-YC(I)
#        endif
       END IF
     END IF

     X1 = X1/1000.0_SP
     X2 = X2/1000.0_SP
     X3 = X3/1000.0_SP
     Y1 = Y1/1000.0_SP
     Y2 = Y2/1000.0_SP
     Y3 = Y3/1000.0_SP

     DELT = (X1*Y2-X2*Y1)**2+(X1*Y3-X3*Y1)**2+(X2*Y3-X3*Y2)**2
     DELT = DELT*1000._SP

     A1U_DAM(I,1) = (Y1+Y2+Y3)*(X1*Y1+X2*Y2+X3*Y3)- &
                    (X1+X2+X3)*(Y1**2+Y2**2+Y3**2)
     A1U_DAM(I,1) = A1U_DAM(I,1)/DELT
     A1U_DAM(I,2) = (Y1**2+Y2**2+Y3**2)*X1-(X1*Y1+X2*Y2+X3*Y3)*Y1
     A1U_DAM(I,2) = A1U_DAM(I,2)/DELT
     A1U_DAM(I,3) = (Y1**2+Y2**2+Y3**2)*X2-(X1*Y1+X2*Y2+X3*Y3)*Y2
     A1U_DAM(I,3) = A1U_DAM(I,3)/DELT
     A1U_DAM(I,4) = (Y1**2+Y2**2+Y3**2)*X3-(X1*Y1+X2*Y2+X3*Y3)*Y3
     A1U_DAM(I,4) = A1U_DAM(I,4)/DELT

!     A2U_DAM(I,1) = (X1+X2+X3)*(X1*X1+X2*X2+X3*X3)- &
!                    (Y1+Y2+Y3)*(X1**2+X2**2+X3**2)
     A2U_DAM(I,1) = (X1+X2+X3)*(X1*Y1+X2*Y2+X3*Y3)- &
                    (Y1+Y2+Y3)*(X1**2+X2**2+X3**2)
     A2U_DAM(I,1) = A2U_DAM(I,1)/DELT
     A2U_DAM(I,2) = (X1**2+X2**2+X3**2)*Y1-(X1*Y1+X2*Y2+X3*Y3)*X1
     A2U_DAM(I,2) = A2U_DAM(I,2)/DELT
     A2U_DAM(I,3) = (X1**2+X2**2+X3**2)*Y2-(X1*Y1+X2*Y2+X3*Y3)*X2
     A2U_DAM(I,3) = A2U_DAM(I,3)/DELT
     A2U_DAM(I,4) = (X1**2+X2**2+X3**2)*Y3-(X1*Y1+X2*Y2+X3*Y3)*X3
     A2U_DAM(I,4) = A2U_DAM(I,4)/DELT
!     if((egid(i)==6794.and.egid(i1)==6793).or.(egid(i)=&
!          &=6793.and.egid(i1)==6794))then
!       print*,'ia_a1u:',egid(i),A1U_DAM(I,1),A1U_DAM(I,2),A1U_DAM(I,3),A1U_DAM(I,4)
!       print*,'ia_a2u:',egid(i1),A2U_DAM(I,1),A2U_DAM(I,2),A2U_DAM(I,3),A2U_DAM(I,4)
!     endif

#    if defined (SPHERICAL)
     X1 = VX(NV(I,1))
     X2 = VX(NV(I,2))
     X3 = VX(NV(I,3))
     Y1 = VY(NV(I,1))
     Y2 = VY(NV(I,2))
     Y3 = VY(NV(I,3))

     AI1 = TPI*(Y2-Y3)
     AI2 = TPI*(Y3-Y1)
     AI3 = TPI*(Y1-Y2)
     CALL ARCX(X2,Y2,X3,Y3,SIDE)
     BI1 = SIDE
     CALL ARCX(X3,Y3,X1,Y1,SIDE)
     BI2 = SIDE
     CALL ARCX(X1,Y1,X2,Y2,SIDE)
     BI3 = SIDE

     X2_DP = XC(I)
     Y2_DP = YC(I)
     CALL ARCC(X1,Y1,X2_DP,Y2_DP,XXC1,YYC1)
     CALL ARCC(X2,Y2,X2_DP,Y2_DP,XXC2,YYC2)
     CALL ARCC(X3,Y3,X2_DP,Y2_DP,XXC3,YYC3)

     CI1 = TPI*(X2-XC(I))*TPI*(Y3-YC(I))*COS(DEG2RAD*YYC2)-&
           TPI*(X3-XC(I))*TPI*(Y2-YC(I))*COS(DEG2RAD*YYC3)

     CI2 = TPI*(X3-XC(I))*TPI*(Y1-YC(I))*COS(DEG2RAD*YYC3)-&
           TPI*(X1-XC(I))*TPI*(Y3-YC(I))*COS(DEG2RAD*YYC1)

     CI3 = TPI*(X1-XC(I))*TPI*(Y2-YC(I))*COS(DEG2RAD*YYC1)-&
           TPI*(X2-XC(I))*TPI*(Y1-YC(I))*COS(DEG2RAD*YYC2)
# else
     X1 = VX(NV(I,1))-XC(I)
     X2 = VX(NV(I,2))-XC(I)
     X3 = VX(NV(I,3))-XC(I)
     Y1 = VY(NV(I,1))-YC(I)
     Y2 = VY(NV(I,2))-YC(I)
     Y3 = VY(NV(I,3))-YC(I)


     AI1 = Y2-Y3
     AI2 = Y3-Y1
     AI3 = Y1-Y2
     BI1 = X3-X2
     BI2 = X1-X3
     BI3 = X2-X1
     CI1 = X2*Y3-X3*Y2
     CI2 = X3*Y1-X1*Y3
     CI3 = X1*Y2-X2*Y1
# endif

     AW0_DAM(I,1) = -CI1/2./ART(I)
     AW0_DAM(I,2) = -CI2/2./ART(I)
     AW0_DAM(I,3) = -CI3/2./ART(I)
     AWX_DAM(I,1) = -AI1/2./ART(I)
     AWX_DAM(I,2) = -AI2/2./ART(I)
     AWX_DAM(I,3) = -AI3/2./ART(I)
     AWY_DAM(I,1) = -BI1/2./ART(I)
     AWY_DAM(I,2) = -BI2/2./ART(I)
     AWY_DAM(I,3) = -BI3/2./ART(I)

!-----other side of dam -------------------------------
     I  = I_CELL_DAM_N(II,2)
     I1 = I_CELL_DAM_N(II,1)

     IF(ISBCE(I) == 1) THEN
       DO J = 1, 3
         IF(NBE(I,J) == 0) JJ = J
       END DO
       J1 = JJ+1-INT((JJ+1)/4)*3
       J2 = JJ+2-INT((JJ+2)/4)*3

       DELTX = VX(NV(I,J1))-VX(NV(I,J2))
# if defined (SPHERICAL)
       IF(DELTX > 180.0_SP)THEN
         DELTX = -360.0_SP+DELTX
       ELSE IF(DELTX < -180.0_SP)THEN
         DELTX =  360.0_SP+DELTX	 
       END IF	 
# endif       
       DELTY = VY(NV(I,J1))-VY(NV(I,J2))

       IF(JJ == 1)THEN
#        if defined (SPHERICAL)
         Y1 = YC(I1)-YC(I)
         Y2 = YC(NBE(I,J1))-YC(I)
         Y3 = YC(NBE(I,J2))-YC(I)
         Y1 = TPI*Y1
         Y2 = TPI*Y2
         Y3 = TPI*Y3

         X1_DP = XC(I)
         Y1_DP = YC(I)
         X2_DP = XC(I1)
         Y2_DP = YC(I1)
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X1 = SIDE

         X2_DP = XC(NBE(I,J1))
         Y2_DP = YC(NBE(I,J1))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X2 = SIDE

         X2_DP = XC(NBE(I,J2))
         Y2_DP = YC(NBE(I,J2))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X3 = SIDE
#        else
         X1 = XC(I1)-XC(I)
         Y1 = YC(I1)-YC(I)
         X2 = XC(NBE(I,J1))-XC(I)
         Y2 = YC(NBE(I,J1))-YC(I)
         X3 = XC(NBE(I,J2))-XC(I)
         Y3 = YC(NBE(I,J2))-YC(I)
#        endif
       ELSE IF(JJ == 2)THEN
#        if defined (SPHERICAL)
         Y1 = YC(NBE(I,J2))-YC(I)
         Y2 = YC(I1)-YC(I)
         Y3 = YC(NBE(I,J1))-YC(I)
         Y1 = TPI*Y1
         Y2 = TPI*Y2
         Y3 = TPI*Y3

         X1_DP = XC(I)
         Y1_DP = YC(I)
         X2_DP = XC(NBE(I,J2))
         Y2_DP = YC(NBE(I,J2))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X1 = SIDE

         X2_DP = XC(I1)
         Y2_DP = YC(I1)
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X2 = SIDE

         X2_DP = XC(NBE(I,J1))
         Y2_DP = YC(NBE(I,J1))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X3 = SIDE
#        else
         X1 = XC(NBE(I,J2))-XC(I)
         Y1 = YC(NBE(I,J2))-YC(I)
         X2 = XC(I1)-XC(I)
         Y2 = YC(I1)-YC(I)
         X3 = XC(NBE(I,J1))-XC(I)
         Y3 = YC(NBE(I,J1))-YC(I)
#        endif
       ELSE
#        if defined (SPHERICAL)
         Y1 = YC(NBE(I,J1))-YC(I)
         Y2 = YC(NBE(I,J2))-YC(I)
         Y3 = YC(I1)-YC(I)
         Y1 = TPI*Y1
         Y2 = TPI*Y2
         Y3 = TPI*Y3

         X1_DP = XC(I)
         Y1_DP = YC(I)
         X2_DP = XC(NBE(I,J1))
         Y2_DP = YC(NBE(I,J1))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X1 = SIDE

         X2_DP = XC(NBE(I,J2))
         Y2_DP = YC(NBE(I,J2))
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X2 = SIDE

         X2_DP = XC(I1)
         Y2_DP = YC(I1)
         CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
         X3 = SIDE
#        else
         X1 = XC(NBE(I,J1))-XC(I)
         y1 = YC(NBE(I,J1))-YC(I)
         X2 = XC(NBE(I,J2))-XC(I)
         Y2 = YC(NBE(I,J2))-YC(I)
         X3 = XC(I1)-XC(I)
         Y3 = YC(I1)-YC(I)
#        endif
       END IF
     END IF

     X1 = X1/1000.0_SP
     X2 = X2/1000.0_SP
     X3 = X3/1000.0_SP
     Y1 = Y1/1000.0_SP
     Y2 = Y2/1000.0_SP
     Y3 = Y3/1000.0_SP

     DELT = (X1*Y2-X2*Y1)**2+(X1*Y3-X3*Y1)**2+(X2*Y3-X3*Y2)**2
     DELT = DELT*1000._SP

     A1U_DAM(I,1) = (Y1+Y2+Y3)*(X1*Y1+X2*Y2+X3*Y3)- &
                    (X1+X2+X3)*(Y1**2+Y2**2+Y3**2)
     A1U_DAM(I,1) = A1U_DAM(I,1)/DELT
     A1U_DAM(I,2) = (Y1**2+Y2**2+Y3**2)*X1-(X1*Y1+X2*Y2+X3*Y3)*Y1
     A1U_DAM(I,2) = A1U_DAM(I,2)/DELT
     A1U_DAM(I,3) = (Y1**2+Y2**2+Y3**2)*X2-(X1*Y1+X2*Y2+X3*Y3)*Y2
     A1U_DAM(I,3) = A1U_DAM(I,3)/DELT
     A1U_DAM(I,4) = (Y1**2+Y2**2+Y3**2)*X3-(X1*Y1+X2*Y2+X3*Y3)*Y3
     A1U_DAM(I,4) = A1U_DAM(I,4)/DELT

!     A2U_DAM(I,1) = (X1+X2+X3)*(X1*X1+X2*X2+X3*X3)- &
!                    (Y1+Y2+Y3)*(X1**2+X2**2+X3**2)
     A2U_DAM(I,1) = (X1+X2+X3)*(X1*Y1+X2*Y2+X3*Y3)- &
                    (Y1+Y2+Y3)*(X1**2+X2**2+X3**2)
     A2U_DAM(I,1) = A2U_DAM(I,1)/DELT
     A2U_DAM(I,2) = (X1**2+X2**2+X3**2)*Y1-(X1*Y1+X2*Y2+X3*Y3)*X1
     A2U_DAM(I,2) = A2U_DAM(I,2)/DELT
     A2U_DAM(I,3) = (X1**2+X2**2+X3**2)*Y2-(X1*Y1+X2*Y2+X3*Y3)*X2
     A2U_DAM(I,3) = A2U_DAM(I,3)/DELT
     A2U_DAM(I,4) = (X1**2+X2**2+X3**2)*Y3-(X1*Y1+X2*Y2+X3*Y3)*X3
     A2U_DAM(I,4) = A2U_DAM(I,4)/DELT
!     if((egid(i)==6794.and.egid(i1)==6793).or.(egid(i)=&
!          &=6793.and.egid(i1)==6794))then
!       print*,'ia_a1u:',egid(i),A1U_DAM(I,1),A1U_DAM(I,2),A1U_DAM(I,3),A1U_DAM(I,4)
!       print*,'ia_a2u:',egid(i1),A2U_DAM(I,1),A2U_DAM(I,2),A2U_DAM(I,3),A2U_DAM(I,4)
!     endif
#    if defined (SPHERICAL)
     X1 = VX(NV(I,1))
     X2 = VX(NV(I,2))
     X3 = VX(NV(I,3))
     Y1 = VY(NV(I,1))
     Y2 = VY(NV(I,2))
     Y3 = VY(NV(I,3))

     AI1 = TPI*(Y2-Y3)
     AI2 = TPI*(Y3-Y1)
     AI3 = TPI*(Y1-Y2)
     CALL ARCX(X2,Y2,X3,Y3,SIDE)
     BI1 = SIDE
     CALL ARCX(X3,Y3,X1,Y1,SIDE)
     BI2 = SIDE
     CALL ARCX(X1,Y1,X2,Y2,SIDE)
     BI3 = SIDE

     X2_DP = XC(I)
     Y2_DP = YC(I)
     CALL ARCC(X1,Y1,X2_DP,Y2_DP,XXC1,YYC1)
     CALL ARCC(X2,Y2,X2_DP,Y2_DP,XXC2,YYC2)
     CALL ARCC(X3,Y3,X2_DP,Y2_DP,XXC3,YYC3)

     CI1 = TPI*(X2-XC(I))*TPI*(Y3-YC(I))*COS(DEG2RAD*YYC2)-&
           TPI*(X3-XC(I))*TPI*(Y2-YC(I))*COS(DEG2RAD*YYC3)

     CI2 = TPI*(X3-XC(I))*TPI*(Y1-YC(I))*COS(DEG2RAD*YYC3)-&
           TPI*(X1-XC(I))*TPI*(Y3-YC(I))*COS(DEG2RAD*YYC1)

     CI3 = TPI*(X1-XC(I))*TPI*(Y2-YC(I))*COS(DEG2RAD*YYC1)-&
           TPI*(X2-XC(I))*TPI*(Y1-YC(I))*COS(DEG2RAD*YYC2)
# else
     X1 = VX(NV(I,1))-XC(I)
     X2 = VX(NV(I,2))-XC(I)
     X3 = VX(NV(I,3))-XC(I)
     Y1 = VY(NV(I,1))-YC(I)
     Y2 = VY(NV(I,2))-YC(I)
     Y3 = VY(NV(I,3))-YC(I)


     AI1 = Y2-Y3
     AI2 = Y3-Y1
     AI3 = Y1-Y2
     BI1 = X3-X2
     BI2 = X1-X3
     BI3 = X2-X1
     CI1 = X2*Y3-X3*Y2
     CI2 = X3*Y1-X1*Y3
     CI3 = X1*Y2-X2*Y1
# endif

     AW0_DAM(I,1) = -CI1/2./ART(I)
     AW0_DAM(I,2) = -CI2/2./ART(I)
     AW0_DAM(I,3) = -CI3/2./ART(I)
     AWX_DAM(I,1) = -AI1/2./ART(I)
     AWX_DAM(I,2) = -AI2/2./ART(I)
     AWX_DAM(I,3) = -AI3/2./ART(I)
     AWY_DAM(I,1) = -BI1/2./ART(I)
     AWY_DAM(I,2) = -BI2/2./ART(I)
     AWY_DAM(I,3) = -BI3/2./ART(I)
    END IF
   END DO

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: shape_ceof_dam " 

   RETURN
   END SUBROUTINE SHAPE_COEF_DAM
!=============================================================================|
!
!=============================================================================|
   SUBROUTINE GET_KDAM
   
   IMPLICIT NONE
   INTEGER  :: I,I1,I2,I3,I4,J,K
   REAL(SP) :: DTMP

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: get_kdam "    

   DO I=1,NODE_DAM1_N
     I1 = I_NODE_DAM1_N(I,1)
     I2 = I_NODE_DAM1_N(I,2)
     D(I1) = H(I1)+ELF(I1)
     D(I2) = H(I2)+ELF(I2)
     
     KDAM(I1) = 0
     DTMP = 0.0_SP
     DO K=1,KB
       DTMP = DTMP + DZ(I1,K)*D(I1)
       IF(DTMP <= D(I1)-D_DAM1_N(I,1))THEN
         KDAM(I1) = KDAM(I1) + 1
       END IF
     END DO  	 
     
     KDAM(I2) = 0
     DTMP = 0.0_SP
     DO K=1,KB
       DTMP = DTMP + DZ(I2,K)*D(I2)
       IF(DTMP <= D(I2)-D_DAM1_N(I,2))THEN
         KDAM(I2) = KDAM(I2) + 1
       END IF
     END DO  	 
   END DO  
   
   DO I=1,NODE_DAM2_N
     I1 = I_NODE_DAM2_N(I,1)
     I2 = I_NODE_DAM2_N(I,2)
     I3 = I_NODE_DAM2_N(I,3)
     D(I1) = H(I1)+ELF(I1)
     D(I2) = H(I2)+ELF(I2)
     D(I3) = H(I3)+ELF(I3)
     
     KDAM(I1) = 0
     DTMP = 0.0_SP
     DO K=1,KB
       DTMP = DTMP + DZ(I1,K)*D(I1)
       IF(DTMP <= D(I1)-D_DAM2_N(I,1))THEN
         KDAM(I1) = KDAM(I1) + 1
       END IF
     END DO  	 
     
     KDAM(I2) = 0
     DTMP = 0.0_SP
     DO K=1,KB
       DTMP = DTMP + DZ(I2,K)*D(I2)
       IF(DTMP <= D(I2)-D_DAM2_N(I,2))THEN
         KDAM(I2) = KDAM(I2) + 1
       END IF
     END DO  	 
     
     KDAM(I3) = 0
     DTMP = 0.0_SP
     DO K=1,KB
       DTMP = DTMP + DZ(I3,K)*D(I3)
       IF(DTMP <= D(I3)-D_DAM2_N(I,3))THEN
         KDAM(I3) = KDAM(I3) + 1
       END IF
     END DO  	 
   END DO  
   
   DO I=1,NODE_DAM3_N
     I1 = I_NODE_DAM3_N(I,1)
     I2 = I_NODE_DAM3_N(I,2)
     I3 = I_NODE_DAM3_N(I,3)
     I4 = I_NODE_DAM3_N(I,4)
     D(I1) = H(I1)+ELF(I1)
     D(I2) = H(I2)+ELF(I2)
     D(I3) = H(I3)+ELF(I3)
     D(I4) = H(I4)+ELF(I4)
     
     KDAM(I1) = 0
     DTMP = 0.0_SP
     DO K=1,KB
       DTMP = DTMP + DZ(I1,K)*D(I1)
       IF(DTMP <= D(I1)-D_DAM3_N(I,1))THEN
         KDAM(I1) = KDAM(I1) + 1
       END IF
     END DO  	 
     
     KDAM(I2) = 0
     DTMP = 0.0_SP
     DO K=1,KB
       DTMP = DTMP + DZ(I2,K)*D(I2)
       IF(DTMP <= D(I2)-D_DAM3_N(I,2))THEN
         KDAM(I2) = KDAM(I2) + 1
       END IF
     END DO  	 
     
     KDAM(I3) = 0
     DTMP = 0.0_SP
     DO K=1,KB
       DTMP = DTMP + DZ(I3,K)*D(I3)
       IF(DTMP <= D(I3)-D_DAM3_N(I,3))THEN
         KDAM(I3) = KDAM(I3) + 1
       END IF
     END DO  	 
     
     KDAM(I4) = 0
     DTMP = 0.0_SP
     DO K=1,KB
       DTMP = DTMP + DZ(I4,K)*D(I4)
       IF(DTMP <= D(I4)-D_DAM3_N(I,4))THEN
         KDAM(I4) = KDAM(I4) + 1
       END IF
     END DO  	 
   END DO  
   
   DO I=1,CELL_DAM_N
     I1 = I_CELL_DAM_N(I,1)
     I2 = I_CELL_DAM_N(I,2)
     D1(I1) = H1(I1)+ELF1(I1)
     D1(I2) = H1(I2)+ELF1(I2)
     
     NBE_DAM(I1) = I2
     NBE_DAM(I2) = I1
     
     KDAM1(I1)=1000
     KDAM1(I2)=1000
     DO J=1,3
       IF(KDAM(NV(I1,J)) >= 1) KDAM1(I1) = MIN(KDAM(NV(I1,J)),KDAM1(I1))
       IF(KDAM(NV(I2,J)) >= 1) KDAM1(I2) = MIN(KDAM(NV(I2,J)),KDAM1(I2))
     END DO 
     IF(KDAM1(I1) == 1000) KDAM1(I1) = 0 
     IF(KDAM1(I2) == 1000) KDAM1(I2) = 0 
   END DO  

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: get_kdam " 

   RETURN
   END SUBROUTINE GET_KDAM   
!============================================================================|
!
!============================================================================|
   SUBROUTINE GET_KDAM_INTERNAL
   IMPLICIT NONE
   
   INTEGER  :: I,I1,I2,I3,I4,J,K
   REAL(SP) :: DTMP

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: get_kdam_internal" 
   
   DO I=1,NODE_DAM1_N
     I1 = I_NODE_DAM1_N(I,1)
     I2 = I_NODE_DAM1_N(I,2)

     IF(KDAM(I1) > 0 .AND. (1-DZ(I1,1))*DT(I1) < D_DAM1_N(I,1))THEN
       KDAM(I1) = 0
     END IF

     IF(KDAM(I2) > 0 .AND. (1-DZ(I2,1))*DT(I2) < D_DAM1_N(I,2))THEN
       KDAM(I2) = 0
     END IF
   END DO  	 
     
   DO I=1,NODE_DAM2_N
     I1 = I_NODE_DAM2_N(I,1)
     I2 = I_NODE_DAM2_N(I,2)
     I3 = I_NODE_DAM2_N(I,3)
     
     IF(KDAM(I1) > 0 .AND. (1-DZ(I1,1))*DT(I1) < D_DAM2_N(I,1))THEN
       KDAM(I1) = 0
     END IF

     IF(KDAM(I2) > 0 .AND. (1-DZ(I2,1))*DT(I2) < D_DAM2_N(I,2))THEN
       KDAM(I2) = 0
     END IF

     IF(KDAM(I3) > 0 .AND. (1-DZ(I3,1))*DT(I3) < D_DAM2_N(I,3))THEN
       KDAM(I3) = 0
     END IF
   END DO  
   
   DO I=1,NODE_DAM3_N
     I1 = I_NODE_DAM3_N(I,1)
     I2 = I_NODE_DAM3_N(I,2)
     I3 = I_NODE_DAM3_N(I,3)
     I4 = I_NODE_DAM3_N(I,4)

     IF(KDAM(I1) > 0 .AND. (1-DZ(I1,1))*DT(I1) < D_DAM3_N(I,1))THEN
       KDAM(I1) = 0
     END IF

     IF(KDAM(I2) > 0 .AND. (1-DZ(I2,1))*DT(I2) < D_DAM3_N(I,2))THEN
       KDAM(I2) = 0
     END IF

     IF(KDAM(I3) > 0 .AND. (1-DZ(I3,1))*DT(I3) < D_DAM3_N(I,3))THEN
       KDAM(I3) = 0
     END IF

     IF(KDAM(I4) > 0 .AND. (1-DZ(I4,1))*DT(I4) < D_DAM3_N(I,4))THEN
       KDAM(I4) = 0
     END IF
   END DO  
   
   DO I=1,CELL_DAM_N
     I1 = I_CELL_DAM_N(I,1)
     I2 = I_CELL_DAM_N(I,2)
     
     KDAM1(I1)=1000
     KDAM1(I2)=1000
     DO J=1,3
       IF(KDAM(NV(I1,J)) >= 1) KDAM1(I1) = MIN(KDAM(NV(I1,J)),KDAM1(I1))
       IF(KDAM(NV(I2,J)) >= 1) KDAM1(I2) = MIN(KDAM(NV(I2,J)),KDAM1(I2))
     END DO 
     IF(KDAM1(I1) == 1000) KDAM1(I1) = 0 
     IF(KDAM1(I2) == 1000) KDAM1(I2) = 0 
   END DO  

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: get_kdam_internal" 

   RETURN
   END SUBROUTINE GET_KDAM_INTERNAL

# endif
END MODULE MOD_DAM       
