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

# if defined (PARTITION_SPECIAL)
MODULE MOD_PAR  
   USE MOD_TYPES
   USE MOD_UTILS
   USE LIMS, ONLY : NGL, MGL
   USE MOD_TIME, ONLY : MPI_TIME
   IMPLICIT NONE
   SAVE

!
!--Global Information
!
   INTEGER, POINTER :: EL_PID(:)         !!PROCESSOR OWNER OF GLOBAL ELEMENT

   INTEGER, POINTER :: ELID(:)    !!LOCAL VALUE OF GLOBAL ELEMENT
   INTEGER, POINTER :: NLID(:)    !!LOCAL VALUE OF GLOBAL NODE 
   INTEGER, POINTER :: ELID_X(:)             !!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
   INTEGER, POINTER :: NLID_X(:)             !!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
!
!--Internal Information (Local)
!
   INTEGER, POINTER :: EGID(:)    !!GLOBAL ID OF LOCAL ELEMENT
   INTEGER, POINTER :: NGID(:)    !!GLOBAL ID OF LOCAL NODE 
   INTEGER, POINTER :: EGID_X(:)             !!GLOBAL ID OF LOCAL ELEMENT
   INTEGER, POINTER :: NGID_X(:)             !!GLOBAL ID OF LOCAL NODE 

!
!--Boundary Information: Halo Elements          
!
   INTEGER              :: NHE               !!NUMBER OF HALO ELEMENTS
   INTEGER, POINTER :: HE_LST(:)         !!GLOBAL IDENTITIES OF HALO ELEMENTS
   INTEGER, POINTER :: HE_OWN(:)         !!OWNER OF HALO ELEMENTS

!
!--Boundary Information: Internal Boundary Nodes
!
   INTEGER :: NBN                            !!NUMBER OF BOUNDARY NODES
   INTEGER :: MX_MLT                         !!MAX MULTIPLICITY OF BOUNDARY NODES
   INTEGER, POINTER :: BN_LST(:)         !!GLOBAL IDENTITY OF BOUNDARY NODES
   INTEGER, POINTER :: BN_LOC(:)         !!LOCAL IDENTITY OF BOUNDARY NODES
   INTEGER, POINTER :: BN_MLT(:)         !!MULTIPLICITY OF BOUNDARY NODES
   INTEGER, POINTER :: BN_NEY(:,:)       !!NODE OWNER LIST
   INTEGER, POINTER :: NDE_ID(:)         !! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY

!
!--Boundary Information: Halo Nodes
!
   INTEGER :: NHN  !!NUMBER OF HALO NODES
   INTEGER, POINTER :: HN_LST(:)         !!LIST OF HALO NODES 
   INTEGER, POINTER :: HN_OWN(:)         !!PRIMARY OWNER OF HALO NODES

# if defined (MULTIPROCESSOR)

!
!--Communication Objects    [SIZE: NPROCS] 
!
   TYPE(COMM), POINTER, DIMENSION(:) :: EC,NC,BNC 

!
!--MPI TYPE Objects    [SIZE: NComponents] 
!
   TYPE:: TYPE_DEF
      INTEGER :: NCOMPS
      INTEGER :: OFFSET(100)
      INTEGER :: BLOCKTYPE(100)
      INTEGER :: BLOCKCOUNT(100)
   END TYPE TYPE_DEF

!
!--Maps for Global Array Reconstruction   [SIZE: NPROCS]
!
   TYPE(MAP), POINTER, DIMENSION(:) :: EMAP,NMAP, LSFMAP
   TYPE(MAP), POINTER, DIMENSION(:) :: EXMAP,NXMAP,BCMAP

   TYPE MAPLINK
      TYPE(MAP), POINTER, DIMENSION(:) :: MAP
      TYPE(MAPLINK), POINTER :: NEXT
   END TYPE MAPLINK

   TYPE(MAPLINK) :: HALO_MAPS

   TYPE(MAPLINK) :: INTERNAL_MAPS

!
!--Statistics Calculation   [SIZE: NPROCS]
!

   INTEGER, ALLOCATABLE :: PNE(:)        !!NUMBER OF ELEMENTS IN EACH PROC
   INTEGER, ALLOCATABLE :: PNN(:)        !!NUMBER OF NODES IN EACH PROC
   INTEGER, ALLOCATABLE :: PNHE(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
   INTEGER, ALLOCATABLE :: PNBN(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
   INTEGER, ALLOCATABLE :: PMBM(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
   INTEGER, ALLOCATABLE :: PNHN(:)       !!NUMBER OF HALO NODES IN EACH PROC



   ! DEAL FOR POINTERS
   INTERFACE PDEAL
      MODULE PROCEDURE VEC_INT_PDEAL
      MODULE PROCEDURE ARR_INT_PDEAL
      MODULE PROCEDURE CUB_INT_PDEAL
      MODULE PROCEDURE FDA_INT_PDEAL

      MODULE PROCEDURE VEC_FLT_PDEAL
      MODULE PROCEDURE ARR_FLT_PDEAL
      MODULE PROCEDURE CUB_FLT_PDEAL
      MODULE PROCEDURE FDA_FLT_PDEAL

      MODULE PROCEDURE VEC_DBL_PDEAL
      MODULE PROCEDURE ARR_DBL_PDEAL
      MODULE PROCEDURE CUB_DBL_PDEAL
      MODULE PROCEDURE FDA_DBL_PDEAL
   END INTERFACE
   ! DEAL FOR ALLOCATABLES
   INTERFACE ADEAL
      MODULE PROCEDURE VEC_INT_ADEAL
      MODULE PROCEDURE ARR_INT_ADEAL
      MODULE PROCEDURE CUB_INT_ADEAL
      MODULE PROCEDURE FDA_INT_ADEAL

      MODULE PROCEDURE VEC_FLT_ADEAL
      MODULE PROCEDURE ARR_FLT_ADEAL
      MODULE PROCEDURE CUB_FLT_ADEAL
      MODULE PROCEDURE FDA_FLT_ADEAL

      MODULE PROCEDURE VEC_DBL_ADEAL
      MODULE PROCEDURE ARR_DBL_ADEAL
      MODULE PROCEDURE CUB_DBL_ADEAL
      MODULE PROCEDURE FDA_DBL_ADEAL
   END INTERFACE

   ! COLLECT FOR POINTER
   INTERFACE PCOLLECT
      MODULE PROCEDURE VEC_INT_PCOLLECT
      MODULE PROCEDURE ARR_INT_PCOLLECT
      MODULE PROCEDURE CUB_INT_PCOLLECT
      MODULE PROCEDURE FDA_INT_PCOLLECT

      MODULE PROCEDURE VEC_FLT_PCOLLECT
      MODULE PROCEDURE ARR_FLT_PCOLLECT
      MODULE PROCEDURE CUB_FLT_PCOLLECT
      MODULE PROCEDURE FDA_FLT_PCOLLECT

      MODULE PROCEDURE VEC_DBL_PCOLLECT
      MODULE PROCEDURE ARR_DBL_PCOLLECT
      MODULE PROCEDURE CUB_DBL_PCOLLECT
      MODULE PROCEDURE FDA_DBL_PCOLLECT
   END INTERFACE

   ! COLLECT FOR ALLOCATABLES
   INTERFACE ACOLLECT
      MODULE PROCEDURE VEC_INT_ACOLLECT
      MODULE PROCEDURE ARR_INT_ACOLLECT
      MODULE PROCEDURE CUB_INT_ACOLLECT
      MODULE PROCEDURE FDA_INT_ACOLLECT

      MODULE PROCEDURE VEC_FLT_ACOLLECT
      MODULE PROCEDURE ARR_FLT_ACOLLECT
      MODULE PROCEDURE CUB_FLT_ACOLLECT
      MODULE PROCEDURE FDA_FLT_ACOLLECT

      MODULE PROCEDURE VEC_DBL_ACOLLECT
      MODULE PROCEDURE ARR_DBL_ACOLLECT
      MODULE PROCEDURE CUB_DBL_ACOLLECT
      MODULE PROCEDURE FDA_DBL_ACOLLECT
   END INTERFACE

   INTERFACE PEXCHANGE
      MODULE PROCEDURE VEC_FLT_PEXCHANGE
      MODULE PROCEDURE VEC_INT_PEXCHANGE
      MODULE PROCEDURE VEC_DBL_PEXCHANGE

      MODULE PROCEDURE ARR_FLT_PEXCHANGE
      MODULE PROCEDURE ARR_INT_PEXCHANGE
      MODULE PROCEDURE ARR_DBL_PEXCHANGE
   END INTERFACE

   INTERFACE AEXCHANGE
      MODULE PROCEDURE VEC_FLT_AEXCHANGE
      MODULE PROCEDURE VEC_INT_AEXCHANGE
      MODULE PROCEDURE VEC_DBL_AEXCHANGE

      MODULE PROCEDURE ARR_FLT_AEXCHANGE
      MODULE PROCEDURE ARR_INT_AEXCHANGE
      MODULE PROCEDURE ARR_DBL_AEXCHANGE
   END INTERFACE

   INTERFACE PPRINT
      MODULE PROCEDURE PPRINT_ARR
      MODULE PROCEDURE PPRINT_VEC
   END INTERFACE

   INTERFACE APRINT
      MODULE PROCEDURE APRINT_ARR
      MODULE PROCEDURE APRINT_VEC
   END INTERFACE
!===================================================================================|
   CONTAINS   !!INCLUDED SUBROUTINES FOLLOW
!===================================================================================|

   SUBROUTINE INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID) 
!===================================================================================|
!  INITIALIZE MPI ENVIRONMENT                                                       |
!===================================================================================|
     INTEGER, INTENT(OUT) :: MYID,NPROCS,MSRID
     LOGICAL, INTENT(OUT) :: SERIAL,PAR,MSR
     INTEGER IERR
     
     if(DBG_SET(dbg_sbr)) &
          & write(IPT,*) "STARTING INIT_MPI_ENV"
     
     IERR=0
     
     CALL MPI_INIT(IERR)
     IF(IERR/=0) WRITE(*,*) "BAD MPI_INIT"
     CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
     IF(IERR/=0) WRITE(*,*) "BAD MPI_COMM_RANK"
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
     IF(IERR/=0) WRITE(*,*) "BAD MPI_COMM_SIZE"
     
     MYID = MYID + 1
     MSRID = 1
     IF(NPROCS > 1) SERIAL=.FALSE.
     IF(NPROCS > 1) PAR   =.TRUE.
     IF(MYID /=  1) MSR   =.FALSE.
     
     ! INITIALIZE THE LIST OF MAPS
     nullify(halo_maps%next)
     nullify(halo_maps%map)
     
     nullify(internal_maps%next)
     nullify(internal_maps%map)
     
     ! USE MPI TYPES TO EXCHANGE FVCOM TIME TYPE
     CALL CREATE_MPI_TIME

     if(DBG_SET(dbg_sbr)) &
          & write(IPT,*) "END INIT_MPI_ENV"
     
     RETURN
   END SUBROUTINE INIT_MPI_ENV
!==============================================================================|
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! THE FOLLOWING THREE SUBROUTINES ARE USED TO DEFINE MPI TYPES WHICH
! CAN TRANSMIT FORTRAN TYPE DATA BETWEEN PROCESSORS. USE
! 'ADD_TO_MPI_TYPE' TO ADD COMPONENTS TO AN MPI TYPE. THE ORDER
! MATTERS AND MUST MATCH THE ORDER IN WHICH THE DATA IS DECLARED
! IN THE FORTRAN TYPE. SEE EXAMPLE IN particle.F
!
!
! MPI AND FORTRAN DATA TYPE NAMES
!
!
!$MPI datatype          $FORTRAN datatype
!MPI_INTEGER            INTEGER
!MPI_REAL               REAL
!MPI_REAL8              REAL*8
!MPI_DOUBLE_PRECISION   DOUBLE PRECISION
!MPI_COMPLEX            COMPLEX
!MPI_LOGICAL            LOGICAL
!MPI_CHARACTER          CHARACTER
!MPI_BYTE               -	 
!MPI_PACKED             -
!
! mod_prec.F DECLARES THE FOLLOWING EQUIVELENCE
!
! MPI_F                REAL(SP) (FOR DOUBLE AND SINGLE PRECISION MODELS)
! MPI_DP               REAL(DP) 
!
!==============================================================================|
   FUNCTION INIT_TYPE_DEF() RESULT(DEF)
     IMPLICIT NONE
     TYPE(TYPE_DEF) :: DEF
     DEF%NCOMPS     = 0
     DEF%BLOCKTYPE  = 0
     DEF%BLOCKCOUNT = 0
     DEF%OFFSET     = 0
   END FUNCTION INIT_TYPE_DEF
!==============================================================================|
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
   SUBROUTINE ADD_TO_MPI_TYPE(DEF,MYTYPE,COUNT)
     IMPLICIT NONE

     TYPE(TYPE_DEF) :: DEF
     INTEGER, INTENT(IN) :: MYTYPE
     INTEGER, INTENT(IN) :: COUNT

     INTEGER :: N, IERR
     INTEGER(KIND=MPI_ADDRESS_KIND) :: EXTENT,LBND

     ! INCRIMENT THE NUMBER OF COMPONENTS
     DEF%NCOMPS = DEF%NCOMPS+1

     IF (DEF%NCOMPS .GT. 99) CALL FATAL_ERROR&
          &("ADD_TO_MPI_TYPE: More than 100 components!",&
          & "You must edit the TYPE_DEF and recompile!" )

     N = DEF%NCOMPS

     ! SET THE NEW BLOCKTYPE
     DEF%BLOCKTYPE(N) = MYTYPE 

     ! SET THE NEW BLOCKCOUNT
     DEF%BLOCKCOUNT(N) = COUNT 

     call MPI_TYPE_GET_EXTENT(MYTYPE, LBND, EXTENT, IERR)
     IF(IERR /= 0) CALL FATAL_ERROR&
          & ("ADD_TO_MPI_TYPE: COULD NOT GET EXTENT FOR THE TYPE?")
        
     ! SET THE NEXT OFFSET (LAST ONE IS NOT USED)
     DEF%OFFSET(N+1) = DEF%OFFSET(N)+ COUNT * EXTENT


   END SUBROUTINE ADD_TO_MPI_TYPE
!==============================================================================|
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
 SUBROUTINE CREATE_MPI_TIME
    IMPLICIT NONE
    
    TYPE(TYPE_DEF) :: MPIT

    ! ZERO THE TYPE
    MPIT = INIT_TYPE_DEF()

    !ADD EACH COMPONENT
    CALL ADD_TO_MPI_TYPE(MPIT,MPI_DOUBLE_PRECISION,2) ! 2 long integers

    ! DEFINE AND COMMIT THE TYPE
    MPI_TIME = CREATE_MPI_TYPE(MPIT)

    END SUBROUTINE CREATE_MPI_TIME

!==============================================================================|
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
   FUNCTION CREATE_MPI_TYPE(DEF) RESULT(MYMPITYPE)
     IMPLICIT NONE

     TYPE(TYPE_DEF),INTENT(IN) :: DEF
     INTEGER :: MYMPITYPE
     
     INTEGER :: IERR

     CALL MPI_TYPE_STRUCT(DEF%NCOMPS, DEF%blockcount, DEF%offset, DEF%blocktype, mympitype, ierr) 
     IF(IERR /= 0) CALL FATAL_ERROR&
          & ("CREATE_MPI_TYPE: COULD NOT CREATE MPI_TYPE_STRUCT?")


      CALL MPI_TYPE_COMMIT(MYMPITYPE, ierr)
      IF(IERR /= 0) CALL FATAL_ERROR&
          & ("CREATE_MPI_TYPE: COULD NOT COMMIT MPI TYPE")

   END FUNCTION CREATE_MPI_TYPE
!==============================================================================|
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
  SUBROUTINE DOMDEC(NGL,NVG,NPROCS,EL_PID,MSR)              
!==============================================================================|
# if defined (THIN_DAM)
   USE control, only : casename,input_dir
   USE Mod_Utils, only: pstop
# endif

   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: NGL
   INTEGER, INTENT(IN) :: NVG(0:NGL,4)
   INTEGER, INTENT(IN)  :: NPROCS
   INTEGER, INTENT(OUT) :: EL_PID(NGL)
   LOGICAL, INTENT(IN)  :: MSR
   INTEGER, ALLOCATABLE :: NVT(:)
   INTEGER :: I,J,COUNT,NTEMP,IERR,ii
# if defined (THIN_DAM)
   character(len=120)   :: cellfile,nodefile 
   INTEGER              :: n,corner_proc_id,corner_elem,processor
   logical              :: fexist
   INTEGER ,ALLOCATABLE :: NODE_DAM1(:,:),NODE_DAM2(:,:),NODE_DAM3(:,:)
   INTEGER ,ALLOCATABLE :: CELL_DAM(:,:),CELL_PROC(:),NODE_PROC(:)
   INTEGER ,ALLOCATABLE :: HOST(:)
   REAL(SP),ALLOCATABLE :: D_DAM1(:,:),D_DAM2(:,:),D_DAM3(:,:)
   INTEGER              :: NN_DAM1,NN_DAM2,NN_DAM3,NC_DAM,DAM_NPROCS
   integer(kind=8) ,allocatable :: dxadj(:),dadjncy(:)
!   type(idxtype),allocatable :: dxadj(:),dadjncy(:)
   integer              :: nvertex,nedge,jj,kk,nbnd
   integer              :: edge1(2),edge2(2),edge3(2)
   integer              :: count1,count2,count3
   logical              :: bexist1,bexist2,bexist3

# endif

!==============================================================================|
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING DOMDEC"

!
!----------------CONVERT NVG to 1D-----------------------!
! 

   IF(MSR) THEN
      ALLOCATE(NVT(3*NGL))
      DO I=1,NGL
         DO J = 1,3
            COUNT = (I-1)*3 + J
            ! FVCOM SWITCHES NODE ORDER
            NVT(COUNT) = NVG(I,5-J)
         END DO
      END DO

!# if defined (THIN_DAM)
!     nbnd = 0
!     do i=1,ngl
!        edge1(1)=nvg(i,1);edge1(2)=nvg(i,2)
!        edge2(1)=nvg(i,2);edge2(2)=nvg(i,3)
!        edge3(1)=nvg(i,3);edge3(2)=nvg(i,1)
!        bexist1=.true.; bexist2=.true.; bexist3=.true.
!        do ii=1,ngl          
!          if(i/=ii)then
!             count1=0; count2=0;count3=0
!             do kk=1,3
!                if(nvg(ii,kk)==edge1(1))count1=count1+1
!                if(nvg(ii,kk)==edge1(2))count1=count1+1
!                if(nvg(ii,kk)==edge2(1))count2=count2+1
!                if(nvg(ii,kk)==edge2(2))count2=count2+1
!                if(nvg(ii,kk)==edge3(1))count3=count3+1
!                if(nvg(ii,kk)==edge3(2))count3=count3+1
!             end do
!             if(count1==2)bexist1=.false.
!             if(count2==2)bexist2=.false.
!             if(count3==2)bexist3=.false.
!          end if
!        end do
!        if(bexist1.or.bexist2.or.bexist3)then
!          nbnd = nbnd + 1
!        end if       
!     end do
!     nedge = mgl*10 !(3*ngl+nbnd)/2*2
!     nvertex = mgl*2
!     print*,'nvertex = ',nvertex
!     print*,'nedge   = ',nedge
!     allocate(dxadj(nvertex))
!     allocate(dadjncy(nedge))
!     dxadj = 0
!     dadjncy = 0
!     print*,'before mesh-to-nodal'
!     CALL METIS_MESHTODUAL(NGL,MGL,loc(NVT),1,1,loc(dxadj),loc(dadjncy))
!     call mesh2nodal(NGL,MGL,loc(NVT),1,1,loc(dxadj),loc(dadjncy))
!     print*,'after mesh-to-nodal'
!# endif
      
!
!-------------DECOMPOSE ELEMENTS USING METIS GRAPH PARTITIONING ---------------!
!
      
      
!$      CALL PARTITION(NPROCS,NGL,MAXVAL(NVT),loc(NVT),loc(EL_PID))
!$      EL_PID = EL_PID + 1

# if defined (THIN_DAM)
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
!---read in nodes list. ---------------------
      OPEN(111,FILE=trim(nodefile),status='old')

!---read in type 1 dam. ---------------------
      READ(111,*)
      READ(111,*) NN_DAM1
      print*,'dam1:',nn_dam1 
      ALLOCATE(NODE_DAM1(NN_DAM1,2));    NODE_DAM1 = 0
      ALLOCATE(D_DAM1(NN_DAM1,2));       D_DAM1    = 0
   
      DO I=1,NN_DAM1
        READ(111,*) NODE_DAM1(I,1),NODE_DAM1(I,2),D_DAM1(I,1),D_DAM1(I,2)
      END DO

!---read in type 2 dam. ---------------------   
     READ(111,*)
     READ(111,*) NN_DAM2
     print*,'dam2:',nn_dam2
     ALLOCATE(NODE_DAM2(NN_DAM2,3));    NODE_DAM2 = 0
     ALLOCATE(D_DAM2(NN_DAM2,3));       D_DAM2    = 0
   
     DO I=1,NN_DAM2
       READ(111,*) NODE_DAM2(I,1),NODE_DAM2(I,2),NODE_DAM2(I,3),    &
                 D_DAM2(I,1),D_DAM2(I,2),D_DAM2(I,3)
     END DO
   
!---read in type 3 dam. ---------------------   
     READ(111,*)
     READ(111,*) NN_DAM3
     print*,'dam3:',nn_dam3
     ALLOCATE(NODE_DAM3(NN_DAM2,4));    NODE_DAM3 = 0
     ALLOCATE(D_DAM3(NN_DAM3,4));       D_DAM3    = 0
   
     DO I=1,NN_DAM3
       READ(111,*) NODE_DAM3(I,1),NODE_DAM3(I,2),NODE_DAM3(I,3),NODE_DAM3(I,4), &
                D_DAM3(I,1),D_DAM3(I,2),D_DAM3(I,3),D_DAM3(I,4)
     END DO
   
     CLOSE(111)
!---read in cells list. ---------------------     
      OPEN(111,FILE=trim(cellfile))
      READ(111,*) NC_DAM
   
      ALLOCATE(CELL_DAM(NC_DAM,2));        CELL_DAM = 0
      ALLOCATE(CELL_PROC(NGL));            CELL_PROC = 0
      ALLOCATE(NODE_PROC(MGL));            NODE_PROC = 0

      DO I=1,NC_DAM
        READ(111,*) CELL_DAM(I,1),CELL_DAM(I,2),PROCESSOR
        CELL_PROC(CELL_DAM(I,1))=PROCESSOR
        CELL_PROC(CELL_DAM(I,2))=PROCESSOR
        NODE_PROC(NVG(CELL_DAM(I,1),1))=PROCESSOR
        NODE_PROC(NVG(CELL_DAM(I,1),2))=PROCESSOR
        NODE_PROC(NVG(CELL_DAM(I,1),3))=PROCESSOR
        NODE_PROC(NVG(CELL_DAM(I,2),1))=PROCESSOR
        NODE_PROC(NVG(CELL_DAM(I,2),2))=PROCESSOR
        NODE_PROC(NVG(CELL_DAM(I,2),3))=PROCESSOR
      END DO
      CLOSE(111)
      DAM_NPROCS = MAXVAL(CELL_PROC)

!$      allocate(host(NPROCS));     host = 0
!$!      corner_elem = 1                       !global element # of corner element
!$!      corner_proc_id =  el_pid(corner_elem)  !METIS-assigned owner of corner element

!     partition with decreased CPU number           
      CALL PARTITION(NPROCS-DAM_NPROCS,NGL,MAXVAL(NVT),loc(NVT),loc(EL_PID))
      EL_PID = EL_PID + 1

!     set additional domains for dam cells
      do n=1,ngl
         do i=1,nn_dam1
           if( nvg(n,2)==node_dam1(i,1).or.nvg(n,2)==node_dam1(i,2) &
          &.or.nvg(n,3)==node_dam1(i,1).or.nvg(n,3)==node_dam1(i,2) &
          &.or.nvg(n,4)==node_dam1(i,1).or.nvg(n,4)==node_dam1(i,2) )then
             el_pid(n)=NPROCS-DAM_NPROCS+NODE_PROC(node_dam1(i,1))
!$             host(el_pid(n))=host(el_pid(n))+1
             exit
           end if
         end do
         do i=1,nn_dam2
           if( nvg(n,2)==node_dam2(i,1).or.nvg(n,2)==node_dam2(i,2) &
         & .or.nvg(n,2)==node_dam2(i,3).or.nvg(n,3)==node_dam2(i,1) &
         & .or.nvg(n,3)==node_dam2(i,2).or.nvg(n,3)==node_dam2(i,3) &
         & .or.nvg(n,4)==node_dam2(i,1).or.nvg(n,4)==node_dam2(i,2) &
         & .or.nvg(n,4)==node_dam2(i,3))then
             el_pid(n)=NPROCS-DAM_NPROCS+NODE_PROC(node_dam2(i,1))
!$             host(el_pid(n))=host(el_pid(n))+1
             exit
           end if
         end do
         do i=1,nn_dam3
           if(nvg(n,2)==node_dam3(i,1).or.nvg(n,2)==node_dam3(i,2) &
         &.or.nvg(n,2)==node_dam3(i,3).or.nvg(n,2)==node_dam3(i,4) &
         &.or.nvg(n,3)==node_dam3(i,1).or.nvg(n,3)==node_dam3(i,2) &
         &.or.nvg(n,3)==node_dam3(i,3).or.nvg(n,3)==node_dam3(i,4) &
         &.or.nvg(n,4)==node_dam3(i,1).or.nvg(n,4)==node_dam3(i,2) &
         &.or.nvg(n,4)==node_dam3(i,3).or.nvg(n,4)==node_dam3(i,4) )then
             el_pid(n)=NPROCS-DAM_NPROCS+NODE_PROC(node_dam3(i,1))
!$             host(el_pid(n))=host(el_pid(n))+1
             exit
           end if
         end do
      end do

      deallocate(cell_dam)
      deallocate(NODE_DAM1,NODE_DAM2,NODE_DAM3)
      deallocate(D_DAM1,D_DAM2,D_DAM3)
      deallocate(NODE_PROC,CELL_PROC)

# else
      CALL PARTITION(NPROCS,NGL,MAXVAL(NVT),loc(NVT),loc(EL_PID))
      EL_PID = EL_PID + 1
# endif

      DEALLOCATE(NVT)
   END IF


!---------------------BROADCAST RESULT TO ALL PROCESSORS-----------------------!

   
   CALL MPI_BCAST(EL_PID,NGL,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END DOMDEC"

   END SUBROUTINE DOMDEC
!===================================================================================|

!------------------------------------------------------------------------------
!===================================================================================|
   SUBROUTINE VEC_FLT_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:)

   REAL(SPA), POINTER :: AP(:)
   REAL(SPA), POINTER :: BP(:)
   REAL(SPA), POINTER :: CP(:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE VEC_FLT_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SPA), INTENT(INOUT), POINTER           :: A(:)
   REAL(SPA), INTENT(INOUT), POINTER, OPTIONAL :: B(:)
   REAL(SPA), INTENT(INOUT), POINTER, OPTIONAL :: C(:)
   INTEGER :: NT  
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(SPA), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_FLT_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1)               = ",Size(A,1)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
                & write(IPT,*) "SIZE(B,1)           = ",Size(B,1)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                & write(IPT,*) "SIZE(C,1)           = ",Size(C,1)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("VEC_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS + 1
       LBUF  = NVARS * CM(I)%NRCV
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_REAL,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 

       N1 = N2+1  ; N2 = N1 + LBUF -1
       SBUF(N1:N2) = A(CM(I)%SNDP(:))
       IF(BYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = B(CM(I)%SNDP(:))
       END IF
       IF(CYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = C(CM(I)%SNDP(:))
       END IF

       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_REAL,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP

     N1 = N2+1 ; N2 = N1 + LBUF -1
     A(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     IF(BYES)THEN
        N1 = N2+1; N2 = N1 + LBUF -1
        B(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
     IF(CYES)THEN
        N1 = N2+1 ; N2 = N1 + LBUF -1
        C(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
  END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_FLT_EXCHANGE"

   RETURN
 END SUBROUTINE VEC_FLT_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_INT_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:)

   INTEGER, POINTER :: AP(:)
   INTEGER, POINTER :: BP(:)
   INTEGER, POINTER :: CP(:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL VEC_INT_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL VEC_INT_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL VEC_INT_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL VEC_INT_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE VEC_INT_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_INT_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   INTEGER, INTENT(INOUT), POINTER           :: A(:)
   INTEGER, INTENT(INOUT), POINTER, OPTIONAL :: B(:)
   INTEGER, INTENT(INOUT), POINTER, OPTIONAL :: C(:)
   INTEGER :: NT  
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   INTEGER, ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_INT_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1)               = ",Size(A,1)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
                & write(IPT,*) "SIZE(B,1)           = ",Size(B,1)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                & write(IPT,*) "SIZE(C,1)           = ",Size(C,1)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("VEC_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS + 1
       LBUF  = NVARS * CM(I)%NRCV
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_INTEGER,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 

       N1 = N2+1  ; N2 = N1 + LBUF -1
       SBUF(N1:N2) = A(CM(I)%SNDP(:))
       IF(BYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = B(CM(I)%SNDP(:))
       END IF
       IF(CYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = C(CM(I)%SNDP(:))
       END IF

       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_INTEGER,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP

     N1 = N2+1 ; N2 = N1 + LBUF -1
     A(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     IF(BYES)THEN
        N1 = N2+1; N2 = N1 + LBUF -1
        B(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
     IF(CYES)THEN
        N1 = N2+1 ; N2 = N1 + LBUF -1
        C(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
  END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_INT_EXCHANGE"

   RETURN
 END SUBROUTINE VEC_INT_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_DBL_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:)

   REAL(DP), POINTER :: AP(:)
   REAL(DP), POINTER :: BP(:)
   REAL(DP), POINTER :: CP(:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE VEC_DBL_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(DP), INTENT(INOUT), POINTER           :: A(:)
   REAL(DP), INTENT(INOUT), POINTER, OPTIONAL :: B(:)
   REAL(DP), INTENT(INOUT), POINTER, OPTIONAL :: C(:)
   INTEGER :: NT  
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(DP), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_DBL_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1)               = ",Size(A,1)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
                & write(IPT,*) "SIZE(B,1)           = ",Size(B,1)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                & write(IPT,*) "SIZE(C,1)           = ",Size(C,1)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("VEC_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS + 1
       LBUF  = NVARS * CM(I)%NRCV
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_DP,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 

       N1 = N2+1  ; N2 = N1 + LBUF -1
       SBUF(N1:N2) = A(CM(I)%SNDP(:))
       IF(BYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = B(CM(I)%SNDP(:))
       END IF
       IF(CYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = C(CM(I)%SNDP(:))
       END IF

       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_DP,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP

     N1 = N2+1 ; N2 = N1 + LBUF -1
     A(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     IF(BYES)THEN
        N1 = N2+1; N2 = N1 + LBUF -1
        B(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
     IF(CYES)THEN
        N1 = N2+1 ; N2 = N1 + LBUF -1
        C(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
  END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_DBL_EXCHANGE"

   RETURN
 END SUBROUTINE VEC_DBL_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_FLT_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:,:)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:,:)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:,:)

   REAL(SPA), POINTER :: AP(:,:)
   REAL(SPA), POINTER :: BP(:,:)
   REAL(SPA), POINTER :: CP(:,:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE ARR_FLT_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SPA), INTENT(INOUT), POINTER           :: A(:,:)
   REAL(SPA), INTENT(INOUT), POINTER, OPTIONAL :: B(:,:)
   REAL(SPA), INTENT(INOUT), POINTER, OPTIONAL :: C(:,:)
   INTEGER :: NT,KT
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(SPA), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_FLT_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                  = ",MYID
        write(IPT,*) "NPROCS                = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1),Size(A,2)   = ",Size(A,1),Size(A,2)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
             &write(IPT,*) "SIZE(B,1),Size(B,2)   = ",Size(B,1),Size(B,2)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                &write(IPT,*) "SIZE(C,1),Size(C,2)   = ",Size(C,1),Size(C,2)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("ARR_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     KT = Size(A,2)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1) .OR. KT .NE. Size(B,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1) .OR. KT .NE. Size(C,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS*KT + 1
       LBUF  = NVARS * CM(I)%NRCV *KT
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_REAL,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*KT*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 
       DO J=1,KT
         N1 = N2+1  ; N2 = N1 + LBUF -1
         SBUF(N1:N2) = A(CM(I)%SNDP(:),J)
         IF(BYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = B(CM(I)%SNDP(:),J)
         END IF
         IF(CYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = C(CM(I)%SNDP(:),J)
         END IF
       END DO
       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_REAL,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*KT*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS*KT 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP
     DO J=1,KT
       N1 = N2+1 ; N2 = N1 + LBUF -1
       A(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       IF(BYES)THEN
         N1 = N2+1; N2 = N1 + LBUF -1
         B(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
       IF(CYES)THEN
         N1 = N2+1 ; N2 = N1 + LBUF -1
         C(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
     END DO
   END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_FLT_EXCHANGE"

   RETURN
 END SUBROUTINE ARR_FLT_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_INT_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:,:)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:,:)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:,:)

   INTEGER, POINTER :: AP(:,:)
   INTEGER, POINTER :: BP(:,:)
   INTEGER, POINTER :: CP(:,:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL ARR_INT_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL ARR_INT_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL ARR_INT_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL ARR_INT_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE ARR_INT_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_INT_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   INTEGER, INTENT(INOUT), POINTER           :: A(:,:)
   INTEGER, INTENT(INOUT), POINTER, OPTIONAL :: B(:,:)
   INTEGER, INTENT(INOUT), POINTER, OPTIONAL :: C(:,:)
   INTEGER :: NT,KT
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   INTEGER, ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_INT_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                  = ",MYID
        write(IPT,*) "NPROCS                = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1),Size(A,2)   = ",Size(A,1),Size(A,2)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
             &write(IPT,*) "SIZE(B,1),Size(B,2)   = ",Size(B,1),Size(B,2)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                &write(IPT,*) "SIZE(C,1),Size(C,2)   = ",Size(C,1),Size(C,2)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("ARR_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     KT = Size(A,2)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1) .OR. KT .NE. Size(B,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1) .OR. KT .NE. Size(C,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS*KT + 1
       LBUF  = NVARS * CM(I)%NRCV *KT
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_INTEGER,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*KT*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 
       DO J=1,KT
         N1 = N2+1  ; N2 = N1 + LBUF -1
         SBUF(N1:N2) = A(CM(I)%SNDP(:),J)
         IF(BYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = B(CM(I)%SNDP(:),J)
         END IF
         IF(CYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = C(CM(I)%SNDP(:),J)
         END IF
       END DO
       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_INTEGER,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*KT*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS*KT 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP
     DO J=1,KT
       N1 = N2+1 ; N2 = N1 + LBUF -1
       A(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       IF(BYES)THEN
         N1 = N2+1; N2 = N1 + LBUF -1
         B(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
       IF(CYES)THEN
         N1 = N2+1 ; N2 = N1 + LBUF -1
         C(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
     END DO
   END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_INT_EXCHANGE"

   RETURN
 END SUBROUTINE ARR_INT_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_DBL_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:,:)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:,:)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:,:)

   REAL(DP), POINTER :: AP(:,:)
   REAL(DP), POINTER :: BP(:,:)
   REAL(DP), POINTER :: CP(:,:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE ARR_DBL_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(DP), INTENT(INOUT), POINTER           :: A(:,:)
   REAL(DP), INTENT(INOUT), POINTER, OPTIONAL :: B(:,:)
   REAL(DP), INTENT(INOUT), POINTER, OPTIONAL :: C(:,:)
   INTEGER :: NT,KT
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(DP), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_DBL_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                  = ",MYID
        write(IPT,*) "NPROCS                = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1),Size(A,2)   = ",Size(A,1),Size(A,2)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
             &write(IPT,*) "SIZE(B,1),Size(B,2)   = ",Size(B,1),Size(B,2)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                &write(IPT,*) "SIZE(C,1),Size(C,2)   = ",Size(C,1),Size(C,2)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("ARR_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     KT = Size(A,2)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1) .OR. KT .NE. Size(B,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1) .OR. KT .NE. Size(C,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS*KT + 1
       LBUF  = NVARS * CM(I)%NRCV *KT
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_DP,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*KT*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 
       DO J=1,KT
         N1 = N2+1  ; N2 = N1 + LBUF -1
         SBUF(N1:N2) = A(CM(I)%SNDP(:),J)
         IF(BYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = B(CM(I)%SNDP(:),J)
         END IF
         IF(CYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = C(CM(I)%SNDP(:),J)
         END IF
       END DO
       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_DP,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*KT*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS*KT 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP
     DO J=1,KT
       N1 = N2+1 ; N2 = N1 + LBUF -1
       A(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       IF(BYES)THEN
         N1 = N2+1; N2 = N1 + LBUF -1
         B(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
       IF(CYES)THEN
         N1 = N2+1 ; N2 = N1 + LBUF -1
         C(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
     END DO
   END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_DBL_EXCHANGE"

   RETURN
 END SUBROUTINE ARR_DBL_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE NODE_MATCH(IMATCH,NBN,BN_MLT,BN_LOC,CM,NT,KT,MYID,NPROCS,A,B,C) 
!===================================================================================|
! IMATCH=1:   ENFORCE AGREEMENT OF A,B,C ON BOUNDARY NODES                          |
! IMATCH=0:   ACCUMULATE VALUES OF A,B,C AT BOUNDARY NODES                          |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: IMATCH
   INTEGER, INTENT(IN)             :: NBN
   INTEGER, INTENT(IN)             :: BN_MLT(NBN)
   INTEGER, INTENT(IN)             :: BN_LOC(NBN)
   INTEGER, INTENT(IN)             :: NT,KT,MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SP), INTENT(INOUT)           :: A(0:NT,KT)
   REAL(SP), INTENT(INOUT), OPTIONAL :: B(0:NT,KT)
   REAL(SP), INTENT(INOUT), OPTIONAL :: C(0:NT,KT)
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(SP), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING NODE_MATCH"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "KT                      = ",KT
        write(IPT,*) "SIZE(A,1)            = ",Size(A,1)
        IF(PRESENT(B)) &
             & write(IPT,*) "SIZE(B,1)           = ",Size(B,1)
        IF(PRESENT(C)) &
             & write(IPT,*) "SIZE(C,1)           = ",Size(C,1)

        DO I = 1,NPROCS
        write(IPT,*)"=============MAP INFO=============="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
        write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="

           
     END if



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS*KT + 1
       LBUF  = NVARS * CM(I)%NRCV *KT
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_F,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*KT*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 
       DO J=1,KT
         N1 = N2+1  ; N2 = N1 + LBUF -1
         SBUF(N1:N2) = A(CM(I)%SNDP(:),J)
         IF(BYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = B(CM(I)%SNDP(:),J)
         END IF
         IF(CYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = C(CM(I)%SNDP(:),J)
         END IF
       END DO
       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_F,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*KT*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
!   DO NMSG = 1,TRCV 
   TRCV = 0
   DO LPROC=1,NPROCS
     IF(CM(LPROC)%NRCV > 0)THEN
       TRCV  = TRCV + 1
       CALL MPI_WAIT(IREQR(TRCV),STAT,IERR)
!     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
!     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS*KT 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP
     DO J=1,KT
       N1 = N2+1 ; N2 = N1 + LBUF -1
       A(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) + A(CM(LPROC)%RCVP(:),J) 
       IF(BYES)THEN
         N1 = N2+1; N2 = N1 + LBUF -1
         B(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) + B(CM(LPROC)%RCVP(:),J)
       END IF
       IF(CYES)THEN
         N1 = N2+1 ; N2 = N1 + LBUF -1
         C(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) + C(CM(LPROC)%RCVP(:),J) 
       END IF
     END DO
    END IF
   END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

!===================================================================================|
!  USE MULTIPLICITY OF NODES TO COMPUTE TRUE AVERAGE VALUE                          |
!===================================================================================|
   IF(IMATCH /=1)RETURN

   DO J=1,KT
     DO I=1,NBN
       A( BN_LOC(I),J) = A( BN_LOC(I),J)/FLOAT(BN_MLT(I))
       IF(BYES)B( BN_LOC(I),J) = B( BN_LOC(I),J)/FLOAT(BN_MLT(I))
       IF(CYES)C( BN_LOC(I),J) = C( BN_LOC(I),J)/FLOAT(BN_MLT(I))
     END DO
   END DO

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END NODE_MATCH"

   RETURN
   END SUBROUTINE NODE_MATCH
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_INT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of INTEGERS FROM A GLOBAL VECTOR                                 |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                        |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:) :: AP
   INTEGER,   POINTER,DIMENSION(:) :: AGP
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:) :: A
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN) :: AG
   IF(ALLOCATED(A))  AP => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_INT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE VEC_INT_ADEAL

!===================================================================================|
   SUBROUTINE VEC_INT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of INTEGERS FROM A GLOBAL VECTOR                                 |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                        |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:) :: A
   INTEGER,   POINTER,DIMENSION(:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: SBUF(:),RBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30111 ! UNIQUE TAG FOR VEC_INT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_INT_PDEAL"


   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER")

      
   
   if(DBG_SET(dbg_sbrio)) then
        
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      if(associated(A)) then
         write(IPT,*) "Ubound(A,1)          = ",Ubound(A,1)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "Ubound(AG,1)         = ",Ubound(AG,1)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"=============MAP INFO=============="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)    = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)  = ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if

   if (MYID .EQ. SENDID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR&
           &("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: VEC_INT_PDEAL")
      

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("THE GLOBAL ARRAY UBOUND DOES NOT MATCH THE MAP IN VEC_INT_PDEAL")

         if (NSZE == 0) CYCLE

         
         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: VEC_INT_PDEAL")

            if(Ubound(A,1) < LSZE) & ! ALLOW FOR M index into MT array
                 & CALL FATAL_ERROR("DEALER POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_INT_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I)) = AG(GM(IP)%LOC_2_GL(I))
            END DO

         else

            ALLOCATE(SBUF(NSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("DEALER CAN NOT ALLOCATE MEMORY IN VEC_INT_PDEAL")

            DO I=1,NSZE
               SBUF(I) = AG(GM(IP)%LOC_2_GL(I))
            END DO


            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("Send Error in VEC_INT_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("DEALER CAN NOT DEALLOCATE MEMORY IN VEC_INT_PDEAL")

         end if


      END DO
   else ! IF I AM ONE OF THE RECEIVERS

      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0 ) RETURN

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: VEC_INT_PDEAL")

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_INT_PDEAL")

      ALLOCATE(RBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_INT_PDEAL")
      
      RBUF=0
      CALL MPI_RECV(RBUF,NSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD , STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN VEC_INT_PDEAL")
        
      DO I=1,NSZE
         A(GM(MYID)%LOC_2_GRID(I)) = RBUF(I)
      END DO
      
      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_INT_PDEAL")
      
   end if

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_INT_PDEAL"
 END SUBROUTINE VEC_INT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_INT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF INTEGERS FROM LOCAL VECTORS                                |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN) :: A
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:)   :: AG
   INTEGER,   POINTER,DIMENSION(:) :: AP
   INTEGER,   POINTER,DIMENSION(:) :: AGP
   IF(ALLOCATED(A))   AP  => A
   IF(ALLOCATED(AG))  AGP => AG
   CALL VEC_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE VEC_INT_ACOLLECT
!===================================================================================|
   SUBROUTINE VEC_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF INTEGERS FROM LOCAL VECTORS                                |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:),INTENT(IN) :: A
   INTEGER,   POINTER,DIMENSION(:)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30112 ! UNIQUE TAG FOR VEC_INT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_INT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL COLLECT AS THE RECEIVER: VEC_INT_PCOLLECT")


   if(DBG_SET(dbg_sbrio)) then
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1)            = ",UBOUND(A,1)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1)           = ",UBOUND(AG,1)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
   END if
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: VEC_INT_PCOLLECT")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: VEC_INT_PCOLLECT")

            if(Ubound(A,1)<LSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_INT_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = A(GM(IP)%LOC_2_GRID(I))
            END DO
         else
            ALLOCATE(RBUF(NSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_INT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in VEC_INT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = RBUF(I)
            END DO
            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_INT_PCOLLECT")
         end if
       End DO
   else

      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0)  return
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: VEC_INT_PCOLLECT")

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_INT_PCOLLECT")

      ALLOCATE(SBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_INT_PCOLLECT")

      DO I=1,NSZE
         SBUF(I) = A(GM(MYID)%LOC_2_GRID(I))
      END DO

      DEST = RECVID - 1
      CALL MPI_SEND(SBUF,NSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN VEC_INT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_INT_PCOLLECT")

   end if
   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_INT_PCOLLECT"
 END SUBROUTINE VEC_INT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_INT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:) :: AGP
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:) :: A
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_INT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE ARR_INT_ADEAL
!===================================================================================|
   SUBROUTINE ARR_INT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:) :: A
   INTEGER,   POINTER,DIMENSION(:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: SBUF(:,:),RBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30121 ! UNIQUE TAG FOR ARR_INT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_INT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: ARR_INT_PDEAL")


   if(DBG_SET(dbg_sbrio)) then
         
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2)= ",UBOUND(A,1),UBOUND(A,2)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1), UBOUND(AG,2)= ",UBOUND(AG,1),UBOUND(AG,2)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
   

   if (MYID .EQ. SENDID) then
      if(.not. associated(AG)) CALL FATAL_ERROR&
           &("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: ARR_INT_PDEAL")
      
      PSZE=ubound(AG,2)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("THE GLOBAL ARRAY UBOUND DOES NOT MATCH THE MAP IN ARR_INT_PDEAL")

         if (NSZE == 0) cycle


         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: ARR_INT_PDEAL")
           
            if(Ubound(A,1)<LSZE .or. ubound(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_INT_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I),1:PSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE)
            END DO

         else

            ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_INT_PDEAL")
            
            DO I=1,NSZE
               SBUF(I,1:PSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE)
            END DO

            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in ARR_INT_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_INT_PDEAL")

         end if
      END DO
   else ! I am a receiver
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0 ) RETURN
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: ARR_INT_PDEAL")
      
      PSZE=ubound(A,2)

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_INT_PDEAL")
      
      ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_INT_PDEAL")

      RBUF=0
      CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN ARR_INT_PDEAL")
      
      DO I=1,NSZE
         A(GM(MYID)%LOC_2_Grid(I),1:PSZE) = RBUF(I,1:PSZE)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_INT_PDEAL")
      
     end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_INT_PDEAL"
 END SUBROUTINE ARR_INT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_INT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:) :: AGP
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN) :: A
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE ARR_INT_ACOLLECT
!===================================================================================|
   SUBROUTINE ARR_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:),INTENT(IN) :: A
   INTEGER,   POINTER,DIMENSION(:,:)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: RBUF(:,:),SBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30122 ! UNIQUE TAG FOR ARR_INT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_INT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE RECEIVER: ARR_INT_PCOLLECT")

   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2)= ",UBOUND(A,1),UBOUND(A,2)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1), UBOUND(AG,2)= ",UBOUND(AG,1),UBOUND(AG,2)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   if (MYID .EQ. RECVID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: ARR_INT_PCOLLECT")

     PSZE = UBOUND(AG,2)
      
      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: ARR_INT_PCOLLECT")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_INT_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE) = A(GM(IP)%LOC_2_GRID(I),1:PSZE)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_INT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in ARR_INT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE) = RBUF(I,1:PSZE)
            END DO
           
            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_INT_PCOLLECT")
         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: ARR_INT_PCOLLECT")
      PSZE = UBOUND(A,2)

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_INT_PCOLLECT")

      ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_INT_PCOLLECT")

      DO I=1,NSZE
         SBUF(I,1:PSZE) = A(GM(MYID)%LOC_2_GRID(I),1:PSZE)
      END DO

      CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN ARR_INT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_INT_PCOLLECT")
   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_INT_PCOLLECT"
 END SUBROUTINE ARR_INT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_INT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:,:) :: AGP
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:,:) :: A
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_INT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE CUB_INT_ADEAL
!===================================================================================|
   SUBROUTINE CUB_INT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:) :: A
   INTEGER,   POINTER,DIMENSION(:,:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: SBUF(:,:,:),RBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30131 ! UNIQUE TAG FOR ARR_INT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_INT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: CUB_INT_PDEAL")


   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)   = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
   

   if (MYID .EQ. SENDID) then
      if(.not. associated(AG)) CALL FATAL_ERROR&
           & ("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: CUB_INT_PDEAL")
      
      PSZE=UBOUND(AG,2)
      QSZE=UBOUND(AG,3)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("The global array ubound does not match the map: CUB_INT_PDEAL")
         
         if (NSZE == 0) cycle


         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: CUB_INT_PDEAL")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_INT_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I),1:PSZE,1:QSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE)
            END DO

         else

            ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN CUB_INT_PDEAL")
            
            DO I=1,NSZE
               SBUF(I,1:PSZE,1:QSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE)
            END DO
            
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error: CUB_INT_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_INT_PDEAL")

         end if


      END DO
   else ! I am a receiver
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0 ) RETURN
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: CUB_INT_PDEAL")
      
      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_INT_PDEAL")
      
      PSZE=UBOUND(A,2)
      QSZE=UBOUND(A,3)
      
      ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_INT_PDEAL")
      RBUF=0
      
      CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN CUB_INT_PDEAL")
      
      DO I=1,NSZE
         A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE)=RBUF(I,1:PSZE,1:QSZE)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_INT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_INT_PDEAL"
 END SUBROUTINE CUB_INT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_INT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:,:) :: AGP
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN) :: A
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE CUB_INT_ACOLLECT
!===================================================================================|
   SUBROUTINE CUB_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:),INTENT(IN) :: A
   INTEGER,   POINTER,DIMENSION(:,:,:)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: RBUF(:,:,:),SBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30132 ! UNIQUE TAG FOR CUB_INT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_INT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE RECEIVER: CUB_INT_PCOLLECT")

    if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)   = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if

   if (MYID .EQ. RECVID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: CUB_INT_PCOLLECT")

     PSZE = UBOUND(AG,2)
     QSZE = UBOUND(AG,3)
      
      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: CUB_INT_PCOLLECT")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_INT_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE) = A(GM(IP)%LOC_2_GRID(I),1:PSZE,1:QSZE)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_INT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in CUB_INT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE) = RBUF(I,1:PSZE,1:QSZE)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_INT_PCOLLECT")

         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: CUB_INT_PCOLLECT")
      PSZE = UBOUND(A,2)
      QSZE = UBOUND(A,3)

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_INT_PCOLLECT")
 
      ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN CUB_INT_PCOLLECT")

      DO I=1,NSZE
         SBUF(I,1:PSZE,1:QSZE) = A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE)
      END DO
      
      CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN CUB_INT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_INT_PCOLLECT")

   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_INT_PCOLLECT"
 END SUBROUTINE CUB_INT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE FDA_INT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:,:,:) :: AGP
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:,:,:) :: A
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:,:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL FDA_INT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE FDA_INT_ADEAL
!===================================================================================|
   SUBROUTINE FDA_INT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:,:) :: A
   INTEGER,   POINTER,DIMENSION(:,:,:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: SBUF(:,:,:,:),RBUF(:,:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,RSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30131 ! UNIQUE TAG FOR ARR_INT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING FDA_INT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: FDA_INT_PDEAL")


   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3), UBOUND(A,4)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3),UBOUND(A,4)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)  = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
   

   if (MYID .EQ. SENDID) then
      if(.not. associated(AG)) CALL FATAL_ERROR&
           & ("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: FDA_INT_PDEAL")
      
      PSZE=UBOUND(AG,2)
      QSZE=UBOUND(AG,3)
      RSZE=UBOUND(AG,4)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("The global array ubound does not match the map: FDA_INT_PDEAL")
         
         if (NSZE == 0) cycle


         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: FDA_INT_PDEAL")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE .or. UBOUND(A,4) .NE. RSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_INT_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I),1:PSZE,1:QSZE,1:RSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE)
            END DO

         else

            ALLOCATE(SBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN FDA_INT_PDEAL")
            
            DO I=1,NSZE
               SBUF(I,1:PSZE,1:QSZE,1:RSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE)
            END DO
            
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE*RSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error: FDA_INT_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN FDA_INT_PDEAL")

         end if


      END DO
   else ! I am a receiver
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0 ) RETURN
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: FDA_INT_PDEAL")
      
      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_INT_PDEAL")
      
      PSZE=UBOUND(A,2)
      QSZE=UBOUND(A,3)
      RSZE=UBOUND(A,4)
      
      ALLOCATE(RBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN FDA_INT_PDEAL")
      RBUF=0
      
      CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE*RSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN FDA_INT_PDEAL")
      
      DO I=1,NSZE
         A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE,1:RSZE)=RBUF(I,1:PSZE,1:QSZE,1:RSZE)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN FDA_INT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END FDA_INT_PDEAL"
 END SUBROUTINE FDA_INT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE FDA_INT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:,:,:) :: AGP
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:),INTENT(IN) :: A
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL FDA_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE FDA_INT_ACOLLECT
!===================================================================================|
   SUBROUTINE FDA_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:,:),INTENT(IN) :: A
   INTEGER,   POINTER,DIMENSION(:,:,:,:)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: RBUF(:,:,:,:),SBUF(:,:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,RSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30132 ! UNIQUE TAG FOR FDA_INT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING FDA_INT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE RECEIVER: FDA_INT_PCOLLECT")

    if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3), UBOUND(A,4)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3),UBOUND(A,4)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)  = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if

   if (MYID .EQ. RECVID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: FDA_INT_PCOLLECT")

     PSZE = UBOUND(AG,2)
     QSZE = UBOUND(AG,3)
     RSZE = UBOUND(AG,4)
      
      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: FDA_INT_PCOLLECT")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE .or. UBOUND(A,4) .NE. RSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_INT_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE) = A(GM(IP)%LOC_2_GRID(I),1:PSZE,1:QSZE,1:RSZE)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN FDA_INT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE*RSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in FDA_INT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE) = RBUF(I,1:PSZE,1:QSZE,1:RSZE)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN FDA_INT_PCOLLECT")

         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: FDA_INT_PCOLLECT")
      PSZE = UBOUND(A,2)
      QSZE = UBOUND(A,3)
      RSZE = UBOUND(A,4)

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_INT_PCOLLECT")
 
      ALLOCATE(SBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN FDA_INT_PCOLLECT")

      DO I=1,NSZE
         SBUF(I,1:PSZE,1:QSZE,1:RSZE) = A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE,1:RSZE)
      END DO
      
      CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE*RSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN FDA_INT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN FDA_INT_PCOLLECT")

   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END FDA_INT_PCOLLECT"
 END SUBROUTINE FDA_INT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


!===================================================================================|
   SUBROUTINE VEC_FLT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of REALS FROM A GLOBAL VECTOR                                    |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                         |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE VEC_FLT_ADEAL
!===================================================================================|
   SUBROUTINE VEC_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of REALS FROM A GLOBAL VECTOR                                    |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                         |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:) :: A
   REAL(SPA),   POINTER,DIMENSION(:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: SBUF(:),RBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30211 ! UNIQUE TAG FOR VEC_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_FLT_PDEAL"


   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: VEC_FLT_PDEAL")

      

   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1)          = ",UBOUND(A,1)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1)         = ",UBOUND(AG,1)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"=============MAP INFO=============="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)    = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)  = ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
   END if
   
   if (MYID .EQ. SENDID) then
      if(.not. associated(AG)) CALL FATAL_ERROR&
           &("POINTER (AG) PASSED TO DEAL MUST ASSOCIATED FOR THE DEALER: VEC_FLT_PDEAL")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("THE GLOBAL ARRAY UBOUND DOES NOT MATCH THE MAP IN VEC_FLT_PDEAL")
         
         if (NSZE == 0) cycle


         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: VEC_FLT_PDEAL")

            if(Ubound(A,1)<LSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_FLT_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I)) = AG(GM(IP)%LOC_2_GL(I))
            END DO

         else
            
            ALLOCATE(SBUF(NSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_FLT_PDEAL")
            
            DO I=1,NSZE
               SBUF(I) = AG(GM(IP)%LOC_2_GL(I))
            END DO
            
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("Send Error in VEC_FLT_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_FLT_PDEAL")

         end if
        

      END DO
   else

      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: VEC_FLT_PDEAL")

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_FLT_PDEAL")

      ALLOCATE(RBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_FLT_PDEAL")
      RBUF=0.0_SPA

      CALL MPI_RECV(RBUF,NSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN VEC_FLT_PDEAL")
      
      DO I=1,NSZE
         A(GM(MYID)%LOC_2_GRID(I)) = RBUF(I)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_FLT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_FLT_PDEAL"
 END SUBROUTINE VEC_FLT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_FLT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF REALS FROM LOCAL VECTORS                                   |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:) :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE VEC_FLT_ACOLLECT
!===================================================================================|
   SUBROUTINE VEC_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF REALS FROM LOCAL VECTORS                                   |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:),INTENT(IN) :: A
   REAL(SPA),   POINTER,DIMENSION(:)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30212 ! UNIQUE TAG FOR VEC_FLT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_FLT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL COLLECT AS THE RECEIVER: VEC_FLT_PCOLLECT")

   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1)          = ",UBOUND(A,1)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1)         = ",UBOUND(AG,1)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"=============MAP INFO=============="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)    = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)  = ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
   END if

   if (MYID .EQ. RECVID) then
        if(.not. associated(AG))&
             & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: VEC_FLT_PCOLLECT")
      
      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: VEC_FLT_PCOLLECT")

            if(Ubound(A,1)<LSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_FLT_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = A(GM(IP)%LOC_2_GRID(I))
            END DO
         else
            ALLOCATE(RBUF(NSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_FLT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in VEC_FLT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = RBUF(I)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_FLT_PCOLLECT")
            
         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: VEC_FLT_PCOLLECT")

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_FLT_PCOLLECT")

      ALLOCATE(SBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_FLT_PCOLLECT")

      DO I=1,NSZE
         SBUF(I) = A(GM(MYID)%LOC_2_GRID(I))
      END DO

      CALL MPI_SEND(SBUF,NSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD,IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN VEC_FLT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_FLT_PCOLLECT")
      
   end if
   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_FLT_PCOLLECT"
 END SUBROUTINE VEC_FLT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_FLT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE ARR_FLT_ADEAL
!===================================================================================|
   SUBROUTINE ARR_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: SBUF(:,:),RBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30221 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_FLT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: ARR_FLT_PDEAL")

   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2)= ",UBOUND(A,1),UBOUND(A,2)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1), UBOUND(AG,2)= ",UBOUND(AG,1),UBOUND(AG,2)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if


   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: ARR_FLT_PDEAL")

      PSZE=ubound(AG,2)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("THE GLOBAL ARRAY UBOUND DOES NOT MATCH THE MAP IN ARR_FLT_PDEAL")

         if (NSZE == 0) cycle


         if (IP == MYID) then
            if(.not. associated(A)) CALL FATAL_ERROR&
                 &("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: ARR_FLT_PDEAL")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_FLT_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I),1:PSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE)
            END DO
         else

            ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_FLT_PDEAL")
            
            DO I=1,NSZE
               SBUF(I,1:PSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE)
            END DO
            
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in ARR_FLT_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_FLT_PDEAL")
         end if


      END DO
   else
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: ARR_FLT_PDEAL")

      PSZE=ubound(A,2)

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_FLT_PDEAL")
      
      ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_FLT_PDEAL")
      RBUF=0.0_SPA
     
      CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN ARR_FLT_PDEAL")
      
      DO I=1,NSZE
         A(GM(MYID)%LOC_2_Grid(I),1:PSZE) = RBUF(I,1:PSZE)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_FLT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_FLT_PDEAL"
 END SUBROUTINE ARR_FLT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_FLT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE ARR_FLT_ACOLLECT
!===================================================================================|
   SUBROUTINE ARR_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:),INTENT(IN) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: RBUF(:,:),SBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30222 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_FLT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE RECEIVER: ARR_FLT_PCOLLECT")

   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2)= ",UBOUND(A,1),UBOUND(A,2)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1), UBOUND(AG,2)= ",UBOUND(AG,1),UBOUND(AG,2)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: ARR_FLT_PCOLLECT")
      
      PSZE = UBOUND(AG,2)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: ARR_FLT_PCOLLECT")

            if(Ubound(A,1)<NSZE .or. UBOUND(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_FLT_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE) = A(GM(IP)%LOC_2_GRID(I),1:PSZE)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_FLT_PCOLLECT")
            RBUF = 0.0_SPA

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in ARR_FLT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE) = RBUF(I,1:PSZE)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_FLT_PCOLLECT")

         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

       if (NSZE == 0)  return
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: ARR_FLT_PCOLLECT")
      PSZE = size(A,2)

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_FLT_PCOLLECT")

      ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_FLT_PCOLLECT")

      DO I=1,NSZE
         SBUF(I,1:PSZE) = A(GM(MYID)%LOC_2_GRID(I),1:PSZE)
      END DO

      CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN ARR_FLT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_FLT_PCOLLECT")

   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_FLT_PCOLLECT"
 END SUBROUTINE ARR_FLT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_FLT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE CUB_FLT_ADEAL
!===================================================================================|
   SUBROUTINE CUB_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: SBUF(:,:,:),RBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30231 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_FLT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: CUB_FLT_PDEAL")

   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)   = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if

   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: CUB_FLT_PDEAL")

      PSZE=UBOUND(AG,2)
      QSZE=UBOUND(AG,3)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("The global array ubound does not match the map: CUB_FLT_PDEAL")

         if (NSZE == 0) cycle


         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: CUB_FLT_PDEAL")

            if(Ubound(A,1)<LSZE .or. Size(A,2) .NE. PSZE .or. Size(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_FLT_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I),1:PSZE,1:QSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE)
            END DO

         else

            ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN CUB_FLT_PDEAL")
            
            DO I=1,NSZE
               SBUF(I,1:PSZE,1:QSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE)
            END DO

            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in CUB_FLT_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_FLT_PDEAL")
            
         end if

      END DO
   else
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: CUB_FLT_PDEAL")

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_FLT_PDEAL")

      PSZE=Ubound(A,2)
      QSZE=Ubound(A,3)

      ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_FLT_PDEAL")
      RBUF=0.0_SPA

      CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN CUB_FLT_PDEAL")

      DO I=1,NSZE
         A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE)=RBUF(I,1:PSZE,1:QSZE)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_FLT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_FLT_PDEAL"
 END SUBROUTINE CUB_FLT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_FLT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE CUB_FLT_ACOLLECT
!===================================================================================|
   SUBROUTINE CUB_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:),INTENT(IN) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: RBUF(:,:,:),SBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30232 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_FLT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE RECEIVER: CUB_FLT_PCOLLECT")

      if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)   = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: CUB_FLT_PCOLLECT")
      
      PSZE = UBOUND(AG,2)
      QSZE = UBOUND(AG,3)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: CUB_FLT_PCOLLECT")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_FLT_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE) = A(GM(IP)%LOC_2_GRID(I),1:PSZE,1:QSZE)
            END DO

         else
            ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_FLT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in CUB_FLT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE) = RBUF(I,1:PSZE,1:QSZE)
            END DO
            
            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_FLT_PCOLLECT")
                        
         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0)  return
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: CUB_FLT_PCOLLECT")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_FLT_PCOLLECT")

      PSZE = UBOUND(A,2)
      QSZE = UBOUND(A,3)

      ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_FLT_PCOLLECT")

      DO I=1,NSZE
         SBUF(I,1:PSZE,1:QSZE) = A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE)
      END DO

      CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN CUB_FLT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_FLT_PCOLLECT")
      
   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_FLT_PCOLLECT"
 END SUBROUTINE CUB_FLT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE FDA_FLT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:,:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL FDA_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE FDA_FLT_ADEAL
!===================================================================================|
   SUBROUTINE FDA_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:,:) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:,:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: SBUF(:,:,:,:),RBUF(:,:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,RSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30231 ! UNIQUE TAG FOR FDA_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING FDA_FLT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: FDA_FLT_PDEAL")

   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3), UBOUND(A,4)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3),UBOUND(A,4)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)  = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if

   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: FDA_FLT_PDEAL")

      PSZE=UBOUND(AG,2)
      QSZE=UBOUND(AG,3)
      RSZE=UBOUND(AG,4)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("The global array ubound does not match the map: FDA_FLT_PDEAL")

         if (NSZE == 0) cycle


         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: FDA_FLT_PDEAL")

            if(Ubound(A,1)<LSZE .or. Size(A,2) .NE. PSZE .or. Size(A,3) .NE. QSZE .or. Size(A,4) .NE. RSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_FLT_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I),1:PSZE,1:QSZE,1:RSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE)
            END DO

         else

            ALLOCATE(SBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN FDA_FLT_PDEAL")
            
            DO I=1,NSZE
               SBUF(I,1:PSZE,1:QSZE,1:RSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE)
            END DO

            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE*RSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in FDA_FLT_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN FDA_FLT_PDEAL")
            
         end if

      END DO
   else
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: FDA_FLT_PDEAL")

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_FLT_PDEAL")

      PSZE=Ubound(A,2)
      QSZE=Ubound(A,3)
      RSZE=Ubound(A,4)

      ALLOCATE(RBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN FDA_FLT_PDEAL")
      RBUF=0.0_SPA

      CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE*RSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN FDA_FLT_PDEAL")

      DO I=1,NSZE
         A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE,1:RSZE)=RBUF(I,1:PSZE,1:QSZE,1:RSZE)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN FDA_FLT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END FDA_FLT_PDEAL"
 END SUBROUTINE FDA_FLT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE FDA_FLT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:,:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:),INTENT(IN) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL FDA_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE FDA_FLT_ACOLLECT
!===================================================================================|
   SUBROUTINE FDA_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:,:),INTENT(IN) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:,:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: RBUF(:,:,:,:),SBUF(:,:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,RSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30232 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING FDA_FLT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE RECEIVER: FDA_FLT_PCOLLECT")

      if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3), UBOUND(A,4)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3),UBOUND(A,4)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)  = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: FDA_FLT_PCOLLECT")
      
      PSZE = UBOUND(AG,2)
      QSZE = UBOUND(AG,3)
      RSZE = UBOUND(AG,4)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: FDA_FLT_PCOLLECT")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE .or. UBOUND(A,4) .NE. RSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_FLT_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE) = A(GM(IP)%LOC_2_GRID(I),1:PSZE,1:QSZE,1:RSZE)
            END DO

         else
            ALLOCATE(RBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN FDA_FLT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE*RSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in FDA_FLT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE) = RBUF(I,1:PSZE,1:QSZE,1:RSZE)
            END DO
            
            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN FDA_FLT_PCOLLECT")
                        
         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0)  return
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: FDA_FLT_PCOLLECT")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_FLT_PCOLLECT")

      PSZE = UBOUND(A,2)
      QSZE = UBOUND(A,3)
      RSZE = UBOUND(A,4)

      ALLOCATE(SBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN FDA_FLT_PCOLLECT")

      DO I=1,NSZE
         SBUF(I,1:PSZE,1:QSZE,1:RSZE) = A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE,1:RSZE)
      END DO

      CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE*RSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN FDA_FLT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN FDA_FLT_PCOLLECT")
      
   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END FDA_FLT_PCOLLECT"
 END SUBROUTINE FDA_FLT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_DBL_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of DOUBLES FROM A GLOBAL VECTOR                                  |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                         |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:) :: AP
   REAL(DP),   POINTER,DIMENSION(:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE VEC_DBL_ADEAL
!===================================================================================|
   SUBROUTINE VEC_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of DOUBLES FROM A GLOBAL VECTOR                                  |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                         |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:) :: A
   REAL(DP),   POINTER,DIMENSION(:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: SBUF(:),RBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30311 ! UNIQUE TAG FOR VEC_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_DBL_PDEAL"


   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: VEC_DBL_PDEAL")

   if(DBG_SET(dbg_sbrio)) then
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1)            = ",UBOUND(A,1)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1)           = ",UBOUND(AG,1)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
   END if
   

   if (MYID .EQ. SENDID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREADY BE ASSOCIATED: VEC_DBL_PDEAL")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("THE GLOBAL ARRAY UBOUND DOES NOT MATCH THE MAP IN VEC_DBL_PDEAL")

         if (NSZE == 0) cycle
         

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: VEC_DBL_PDEAL")

            if(Ubound(A,1)<LSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_DBL_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I)) = AG(GM(IP)%LOC_2_GL(I))
            END DO

         else

            ALLOCATE(SBUF(NSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_DBL_PDEAL")
            
            DO I=1,NSZE
               SBUF(I) = AG(GM(IP)%LOC_2_GL(I))
            END DO
            
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("Send Error in VEC_DBL_PDEAL")
            
            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_DBL_PDEAL")

         end if

      END DO
   else
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: VEC_DBL_PDEAL")

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_DBL_PDEAL")
      
      ALLOCATE(RBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_DBL_PDEAL")
      RBUF=0.0_DP

      CALL MPI_RECV(RBUF,NSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN VEC_DBL_PDEAL")
      
      DO I=1,NSZE
         A(GM(MYID)%LOC_2_GRID(I)) = RBUF(I)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_DBL_PDEAL")
      
   end if
   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_DBL_PDEAL"
 END SUBROUTINE VEC_DBL_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_DBL_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF DOUBLES FROM LOCAL VECTORS                                 |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:) :: AP
   REAL(DP),   POINTER,DIMENSION(:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE VEC_DBL_ACOLLECT
!===================================================================================|
   SUBROUTINE VEC_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF DOUBLES FROM LOCAL VECTORS                                 |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:),INTENT(IN) :: A
   REAL(DP),   POINTER,DIMENSION(:)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30312 ! UNIQUE TAG FOR VEC_FLT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_DBL_PCOLLECT"
   
  if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL COLLECT AS THE RECEIVER: VEC_DBL_PCOLLECT")

   if(DBG_SET(dbg_sbrio)) then
        
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      if(associated(A)) then
         write(IPT,*) "Ubound(A,1)          = ",Ubound(A,1)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "Ubound(AG,1)         = ",Ubound(AG,1)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"=============MAP INFO=============="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)    = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)  = ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if

   if (MYID .EQ. RECVID) then
        if(.not. associated(AG))&
             & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: VEC_DBL_PCOLLECT")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: VEC_DBL_PCOLLECT")

            if(Ubound(A,1)<LSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_DBL_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = A(GM(IP)%LOC_2_GRID(I))
            END DO
         else
            ALLOCATE(RBUF(NSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_DBL_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in VEC_DBL_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = RBUF(I)
            END DO
            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_DBL_PCOLLECT")
         end if
       End DO
   else
      
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: VEC_DBL_PCOLLECT")

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: VEC_DBL_PCOLLECT")

      ALLOCATE(SBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_DBL_PCOLLECT")

      DO I=1,NSZE
         SBUF(I) = A(GM(MYID)%LOC_2_GRID(I))
      END DO

      CALL MPI_SEND(SBUF,NSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD,IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN VEC_DBL_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_DBL_PCOLLECT")
   end if
   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_DBL_PCOLLECT"
 END SUBROUTINE VEC_DBL_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
 SUBROUTINE ARR_DBL_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE ARR_DBL_ADEAL
!===================================================================================|
 SUBROUTINE ARR_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:) :: A
   REAL(DP),   POINTER,DIMENSION(:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: SBUF(:,:),RBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30321 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_DBL_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: ARR_DBL_PDEAL")

   if(DBG_SET(dbg_sbrio)) then
         
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2)= ",UBOUND(A,1),UBOUND(A,2)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1), UBOUND(AG,2)= ",UBOUND(AG,1),UBOUND(AG,2)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if

   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: ARR_DBL_PDEAL")

      PSZE=size(AG,2)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE
         
         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("THE GLOBAL ARRAY UBOUND DOES NOT MATCH THE MAP IN ARR_FLT_PDEAL")
         
         if (NSZE == 0)  cycle
         
         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: ARR_DBL_PDEAL")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_DBL_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I),1:PSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE)
            END DO

         else
         
            ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PDEAL")
            
            DO I=1,NSZE
               SBUF(I,1:PSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE)
            END DO
            
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in ARR_DBL_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PDEAL")

         end if

      END DO
   else

      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: ARR_DBL_PDEAL")
      PSZE=UBOUND(A,2)

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_DBL_PDEAL")

      ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PDEAL")
      RBUF=0.0_DP

      CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD,STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN ARR_DBL_PDEAL")

      DO I=1,NSZE
         A(GM(MYID)%LOC_2_Grid(I),1:PSZE) = RBUF(I,1:PSZE)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PDEAL")
   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_DBL_PDEAL"
 END SUBROUTINE ARR_DBL_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_DBL_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE ARR_DBL_ACOLLECT
!===================================================================================|
   SUBROUTINE ARR_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:),INTENT(IN) :: A
   REAL(DP),   POINTER,DIMENSION(:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: RBUF(:,:), SBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30322 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_DBL_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE RECEIVER: ARR_DBL_PCOLLECT")
   
   if(DBG_SET(dbg_sbrio)) then
         
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2)= ",UBOUND(A,1),UBOUND(A,2)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1), UBOUND(AG,2)= ",UBOUND(AG,1),UBOUND(AG,2)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
  
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: ARR_DBL_PCOLLECT")
      
      PSZE = UBOUND(AG,2)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: ARR_DBL_PCOLLECT")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_DBL_PCOLLECT")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE) = A(GM(IP)%LOC_2_GRID(I),1:PSZE)
            END DO

         else
            ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in ARR_DBL_PCOLLECT: ARR_DBL_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE) = RBUF(I,1:PSZE)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PCOLLECT")
         end if
       End DO
   else

      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")
      PSZE = size(A,2)

      if(Ubound(A,1)<LSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: ARR_DBL_PCOLLECT")
      
      ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PCOLLECT")

      DO I=1,NSZE
         SBUF(I,1:PSZE) = A(GM(MYID)%LOC_2_GRID(I),1:PSZE)
      END DO

      CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN ARR_DBL_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PCOLLECT")

   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_DBL_PCOLLECT"
 END SUBROUTINE ARR_DBL_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
 SUBROUTINE CUB_DBL_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE CUB_DBL_ADEAL
!===================================================================================|
 SUBROUTINE CUB_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: A
   REAL(DP),   POINTER,DIMENSION(:,:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: SBUF(:,:,:),RBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30331 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_DBL_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: CUB_DBL_PDEAL")

   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)   = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if

   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: CUB_DBL_PDEAL")

      PSZE=UBOUND(AG,2)
      QSZE=UBOUND(AG,3)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("The global array ubound does not match the map: CUB_DBL_PDEAL")

         if (NSZE == 0)  cycle
         
         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: CUB_DBL_PDEAL")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_DBL_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I),1:PSZE,1:QSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE)
            END DO

         else
         
            ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN CUB_DBL_PDEAL")
            
            DO I=1,NSZE
               SBUF(I,1:PSZE,1:QSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE)
            END DO
            
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in CUB_DBL_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_DBL_PDEAL")

         end if
      END DO
   else

      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      
      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: CUB_DBL_PDEAL")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_DBL_PDEAL")
      
      PSZE=Ubound(A,2)
      QSZE=Ubound(A,3)

      ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_DBL_PDEAL")
      RBUF=0.0_DP

      CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN CUB_DBL_PDEAL")

      DO I=1,NSZE
         A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE)=RBUF(I,1:PSZE,1:QSZE)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_DBL_PDEAL")
   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_DBL_PDEAL"
 END SUBROUTINE CUB_DBL_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_DBL_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE CUB_DBL_ACOLLECT
!===================================================================================|
   SUBROUTINE CUB_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:),INTENT(IN) :: A
   REAL(DP),   POINTER,DIMENSION(:,:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: RBUF(:,:,:),SBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30332 ! UNIQUE TAG FOR CUB_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_DBL_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE RECEIVER: CUB_DBL_PCOLLECT")

   
   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)   = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: CUB_DBL_PCOLLECT")
      
      PSZE = UBOUND(AG,2)
      QSZE = UBOUND(AG,3)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: CUB_DBL_PCOLLECT")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_DBL_PCOLLECT")
                        
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE) = A(GM(IP)%LOC_2_GRID(I),1:PSZE,1:QSZE)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_DBL_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in CUB_DBL_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE) = RBUF(I,1:PSZE,1:QSZE)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_DBL_PCOLLECT")
         end if
       End DO
   else

      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: CUB_DBL_PCOLLECT")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: CUB_DBL_PCOLLECT")

      PSZE = UBOUND(A,2)
      QSZE = UBOUND(A,3)
      ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN CUB_DBL_PCOLLECT")

      DO I=1,NSZE
         SBUF(I,1:PSZE,1:QSZE) = A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE)
      END DO

      CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN CUB_DBL_PCOLLECT")
      
      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_DBL_PCOLLECT")
   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_DBL_PCOLLECT"
 END SUBROUTINE CUB_DBL_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
 SUBROUTINE FDA_DBL_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:,:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL FDA_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE FDA_DBL_ADEAL
!===================================================================================|
 SUBROUTINE FDA_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:,:) :: A
   REAL(DP),   POINTER,DIMENSION(:,:,:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: SBUF(:,:,:,:),RBUF(:,:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,RSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG = 30331 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING FDA_DBL_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE SENDER: FDA_DBL_PDEAL")

   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "SENDID                  = ",SENDID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3), UBOUND(A,4)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3),UBOUND(A,4)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)  = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if

   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST BE ASSOCIATED FOR THE DEALER: FDA_DBL_PDEAL")

      PSZE=UBOUND(AG,2)
      QSZE=UBOUND(AG,3)
      RSZE=UBOUND(AG,4)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         IF(UBOUND(AG,1) .NE. GSZE) CALL FATAL_ERROR&
              &("The global array ubound does not match the map: FDA_DBL_PDEAL")

         if (NSZE == 0)  cycle
         
         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: FDA_DBL_PDEAL")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE .or. UBOUND(A,4) .NE. RSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_DBL_PDEAL")

            DO I=1,NSZE
               A(GM(IP)%LOC_2_Grid(I),1:PSZE,1:QSZE,1:RSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE)
            END DO

         else
         
            ALLOCATE(SBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN FDA_DBL_PDEAL")
            
            DO I=1,NSZE
               SBUF(I,1:PSZE,1:QSZE,1:RSZE) = AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE)
            END DO
            
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE*RSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in FDA_DBL_PDEAL")

            DEALLOCATE(SBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN FDA_DBL_PDEAL")

         end if
      END DO
   else

      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE
      
      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED: FDA_DBL_PDEAL")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_DBL_PDEAL")
      
      PSZE=Ubound(A,2)
      QSZE=Ubound(A,3)
      RSZE=Ubound(A,4)

      ALLOCATE(RBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN FDA_DBL_PDEAL")
      RBUF=0.0_DP

      CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE*RSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN FDA_DBL_PDEAL")

      DO I=1,NSZE
         A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE,1:RSZE)=RBUF(I,1:PSZE,1:QSZE,1:RSZE)
      END DO

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN FDA_DBL_PDEAL")
   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END FDA_DBL_PDEAL"
 END SUBROUTINE FDA_DBL_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE FDA_DBL_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:,:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:),INTENT(IN) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL FDA_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE FDA_DBL_ACOLLECT
!===================================================================================|
   SUBROUTINE FDA_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:,:),INTENT(IN) :: A
   REAL(DP),   POINTER,DIMENSION(:,:,:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: RBUF(:,:,:,:),SBUF(:,:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE,RSZE,LSZE,GSZE
   INTEGER, PARAMETER :: TAG =  30332 ! UNIQUE TAG FOR CUB_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING FDA_DBL_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE RECEIVER: FDA_DBL_PCOLLECT")

   
   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "UBOUND(A,1), UBOUND(A,2), UBOUND(A,3), UBOUND(A,4)   = ",UBOUND(A,1),UBOUND(A,2),UBOUND(A,3),UBOUND(A,4)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)  = ",UBOUND(AG,1),UBOUND(AG,2),UBOUND(AG,3),UBOUND(AG,4)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "GM(I)%LSIZE             = ",GM(I)%LSIZE
         write(IPT,*) "GM(I)%GSIZE             = ",GM(I)%GSIZE
         write(IPT,*) "UBOUND(GM(I)%LOC_2_GL)  = ",UBOUND(GM(I)%LOC_2_GL)
         if(associated(GM(I)%LOC_2_GRID)) write(IPT,*) "UBOUND(GM(I)%LOC_2_GRID)= ",UBOUND(GM(I)%LOC_2_GRID)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED: FDA_DBL_PCOLLECT")
      
      PSZE = UBOUND(AG,2)
      QSZE = UBOUND(AG,3)
      RSZE = UBOUND(AG,4)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         LSZE = GM(IP)%LSIZE
         GSZE = GM(IP)%GSIZE

         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: FDA_DBL_PCOLLECT")

            if(Ubound(A,1)<LSZE .or. UBOUND(A,2) .NE. PSZE .or. UBOUND(A,3) .NE. QSZE .or. UBOUND(A,4) .NE. RSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_DBL_PCOLLECT")
                        
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE) = A(GM(IP)%LOC_2_GRID(I),1:PSZE,1:QSZE,1:RSZE)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN FDA_DBL_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE*RSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in FDA_DBL_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),1:PSZE,1:QSZE,1:RSZE) = RBUF(I,1:PSZE,1:QSZE,1:RSZE)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN FDA_DBL_PCOLLECT")
         end if
       End DO
   else

      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      LSZE = GM(MYID)%LSIZE
      GSZE = GM(MYID)%GSIZE

      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED: FDA_DBL_PCOLLECT")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP: FDA_DBL_PCOLLECT")

      PSZE = UBOUND(A,2)
      QSZE = UBOUND(A,3)
      RSZE = UBOUND(A,4)

      ALLOCATE(SBUF(NSZE,PSZE,QSZE,RSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN FDA_DBL_PCOLLECT")

      DO I=1,NSZE
         SBUF(I,1:PSZE,1:QSZE,1:RSZE) = A(GM(MYID)%LOC_2_GRID(I),1:PSZE,1:QSZE,1:RSZE)
      END DO

      CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE*RSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN FDA_DBL_PCOLLECT")
      
      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN FDA_DBL_PCOLLECT")
   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END FDA_DBL_PCOLLECT"
 END SUBROUTINE FDA_DBL_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE SORT(N,M,SZE)
!===================================================================================|
!    SORT ELEMENTS IN N AND RETURN ORDER IN M                                       |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN)             :: SZE
   INTEGER, INTENT(INOUT)          :: N(SZE)
   INTEGER, INTENT(OUT)            :: M(SZE)
!------------------------------------------------------------------------------
   INTEGER I,LAST,J(1)
   REAL(SP)  VALM
   INTRINSIC MINLOC
!------------------------------------------------------------------------------

   M = 0
   LAST = 0
   DO I=1,SZE
     J = MINLOC(N,MASK=N > LAST)
     M(I) = J(1)
     LAST = N(J(1))
   END DO

   RETURN
   END SUBROUTINE SORT
!===================================================================================|


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   INTEGER FUNCTION GETLOC(GLOC,MAP,N)
!===================================================================================|
!    DETERMINE LOCAL IDENTITY OF ELEMENT/NODE I USING MAP                           |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: GLOC,N
   INTEGER, INTENT(IN) :: MAP(N)
   INTEGER  I,FOUND
!------------------------------------------------------------------------------
   FOUND = 0
   DO I=1,N
     IF(MAP(I)==GLOC) FOUND = I
   END DO

   GETLOC = FOUND
   RETURN 
   END FUNCTION GETLOC
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
   SUBROUTINE ADD_MAP2LIST(HEAD,MP)
     IMPLICIT NONE

     TYPE(MAPLINK), target :: HEAD
     TYPE(MAP),POINTER ::MP(:)
     TYPE(MAPLINK), POINTER :: current, previous
     integer CNT
     
     CNT = 1
     previous => HEAD
     current => previous%next
     
     DO
        IF(.NOT. ASSOCIATED(CURRENT)) EXIT

        CURRENT  => CURRENT%NEXT
        PREVIOUS => PREVIOUS%NEXT
        CNT = CNT +1
        IF(CNT > 100) CALL FATAL_ERROR&
             &("ADD_MAP_TO_LIST: LOOP COUNT EXCEEDED 100; STOP!")

     END DO
     
     allocate(previous%next)
     previous%next%next =>current
     previous%next%MAP => MP

   END SUBROUTINE ADD_MAP2LIST
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
! MAKE A GLOBAL COLLECT/DEAL MAP AND EXCHANGE TO ALL PROCS
!===================================================================================|
    FUNCTION MAKE_MAP(MYID,NPROCS,GSIZE,LSIZE,N2G,N2L) RESULT(MYMAP)
   ! BECAUSE THE HALO IS AT THE END OF THE ARRAY, THERE IS NO NEED
   ! FOR AN EXPLICIT LOCAL MAP. THE DATA ARRAY SENT IN THE MPI CALL
   ! IS ONE TO ONE WITH DATA IN THE LOCAL ARRAY. MAKE_MAP ACCEPTS A
   ! LOCAL INDEX ARRAY THAT MAPS THE DATA SENT IN THE MPI CALL TO THE
   ! LOCAL ARRAY FOR MORE COMPLEX CASES, I.E. SUBDOMAINS.
      IMPLICIT NONE
      
      TYPE(MAP), POINTER, DIMENSION(:) :: MYMAP

      INTEGER,INTENT(IN) :: MYID,NPROCS,GSIZE,LSIZE

      INTEGER, POINTER :: N2G(:)

      INTEGER, OPTIONAL, POINTER :: N2L(:)
      
      INTEGER, POINTER :: TEMP(:)

      INTEGER :: I, J, SENDER,IERR,NSIZE
      
      if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "START MAKE_MAP"

      if(DBG_SET(dbg_sbrio)) then
         write(IPT,*) "myid",myid
         write(IPT,*) "nprocs",nprocs
         write(IPT,*) "gsize",gsize
         write(IPT,*) "lsize",lsize
         write(IPT,*) "N2G-UBOUND",UBOUND(N2G,1)
         IF(PRESENT(N2L)) THEN
            write(IPT,*) "N2L-UBOUND",UBOUND(N2L,1)
         END if
      end if


      ALLOCATE(MYMAP(NPROCS))
      
      ! MAKE SAFE FOR ALLGATHER USING MPI_IO_MODE
      ALLOCATE(TEMP(NPROCS+1))

      MYMAP(:)%NSIZE = 0
      MYMAP(:)%GSIZE = GSIZE
      MYMAP(:)%LSIZE = 0
      
      NSIZE = 0
      IF(associated(N2G)) NSIZE = UBOUND(N2G,1)


      !--Determine Number of Elements for Each Processor
      CALL MPI_ALLGATHER(NSIZE,1,MPI_INTEGER,TEMP,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      MYMAP(:)%NSIZE = TEMP(1:NPROCS)

!!$      DO I=1,NPROCS    
!!$         IF(MYID == I) MYMAP(I)%NSIZE = NSIZE
!!$         SENDER = I - 1
!!$         CALL MPI_BCAST(MYMAP(I)%NSIZE,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
!!$      END DO
      

      IF(NSIZE > LSIZE) CALL FATAL_ERROR &
           &("MAKE_MAP: WHEN CREATING AN MPI MAP TYPE IF THERE IS NO LOCAL_2_GRID INDEX",&
           & "THEN THE Local Size MUST BE GREATER THAN NSIZE, THE UBOUND OF N2G INDEX ARRAY")
      
      ! MAP IS ONE TO ONE SO THE LOCAL SIZE IS THE SAME AS THE LOCAL DATA SIZE
      CALL MPI_ALLGATHER(LSIZE,1,MPI_INTEGER,TEMP,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      MYMAP(:)%LSIZE = TEMP(1:NPROCS)

!!$      DO I=1,NPROCS    
!!$         IF(MYID == I) MYMAP(I)%LSIZE = LSIZE
!!$         SENDER = I - 1
!!$         CALL MPI_BCAST(MYMAP(I)%LSIZE,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
!!$      END DO
      
      !--Construct Mapping Array for Each Processor 
      DO I=1,NPROCS
                  
         NSIZE = MYMAP(I)%NSIZE
         IF(NSIZE > 0) THEN
            ALLOCATE(MYMAP(I)%LOC_2_GL(0:NSIZE))
            MYMAP(I)%LOC_2_GL=0
            
            
            if(myid == I) MYMAP(I)%LOC_2_GL(1:NSIZE) =  N2G(1:NSIZE)
            SENDER = I - 1
            CALL MPI_BCAST(MYMAP(I)%LOC_2_GL(1:NSIZE),NSIZE,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
         ELSE
            NULLIFY(MYMAP(I)%LOC_2_GL)
         END IF

         NULLIFY(MYMAP(I)%LOC_2_GRID)
      END DO


      ! THIS ARRAY IS ONLY NEEDED ON THE LOCAL PROCESSOR - DO NOT USE
      ! MPI TO BROADCAST IT!

      NSIZE = 0
      IF(MYID <= NPROCS) NSIZE = MYMAP(MYID)%NSIZE

      IF(NSIZE > 0) THEN
         IF(PRESENT(N2L)) THEN
            
            ALLOCATE(MYMAP(MYID)%LOC_2_GRID(0:NSIZE))
            MYMAP(MYID)%LOC_2_GRID = 0
            
            MYMAP(MYID)%LOC_2_GRID(1:NSIZE) = N2L(1:NSIZE)
            
         ELSE
            ALLOCATE(MYMAP(MYID)%LOC_2_GRID(0:NSIZE))
            MYMAP(MYID)%LOC_2_GRID = 0
            
            DO J= 1,NSIZE
               MYMAP(MYID)%LOC_2_GRID(J)=J
            END DO
            
         END IF

      END IF

      if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END MAKE_MAP"


    END FUNCTION MAKE_MAP
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   FUNCTION FIND_MAP(HEAD,Gsize,Lsize,FOUND) RESULT(MP)
     IMPLICIT NONE
     TYPE(MAPLINK),target :: HEAD
     TYPE(MAP), POINTER :: MP(:)
     LOGICAL :: FOUND 
     integer, intent(IN) :: Gsize
     integer,POINTER, intent(IN) :: Lsize(:)
     TYPE(MAPLINK), POINTER :: current, previous
     INTEGER CNT

     previous => HEAD
     current => previous%next
     CNT = 0
     FOUND = .FALSE.
     nullify(mp)

     DO
        IF(.NOT. ASSOCIATED(CURRENT)) RETURN

        IF(.NOT. ASSOCIATED(CURRENT%MAP))&
             CALL FATAL_ERROR("FIND_MAP: MAP LIST LINK HAS UNASSOCIATED MAP!")

        IF(DBG_SET(DBG_MPI))THEN
           WRITE(IPT,*) "========================================================="
           WRITE(IPT,*) "SEARCHING MAPS: MAP GLOBAL_SIZEs= ",CURRENT%MAP%GSIZE
           WRITE(IPT,*) "SEARCHING MAPS for Gsize= ",Gsize
           WRITE(IPT,*) "SEARCHING MAPS: MAP LOCAL_SIZEs= ",CURRENT%MAP%LSIZE
           WRITE(IPT,*) "SEARCHING MAPS for Lsize= ",Lsize
        END IF

        if(ALL(CURRENT%MAP%GSIZE == GSIZE) .and. SUM(abs(Lsize-CURRENT%MAP%LSIZE)) == 0 ) then
           MP => CURRENT%MAP
           FOUND = .TRUE.
           RETURN
        end if

        CURRENT  => CURRENT%NEXT
        
        IF (CNT > 100)CALL FATAL_ERROR&
             &("FIND_MAP: LOOP COUNT EXCEEDED 100; STOP!")

     END DO

   END FUNCTION FIND_MAP

   SUBROUTINE PRINT_MAP(MP,MSG)
     IMPLICIT NONE
     TYPE(MAP),POINTER :: MP(:)
     CHARACTER(LEN=*),OPTIONAL :: MSG
     INTEGER :: I, NN,LL

     IF(DBG_SET(DBG_LOG)) THEN
        
        IF(PRESENT(MSG)) THEN
           WRITE(IPT,*) "! ==== PRINT_MAP: "//TRIM(MSG)//" ===="
        ELSE
           WRITE(IPT,*) "! ==== PRINT_MAP ===="
        END IF

        IF(.not.associated(MP)) THEN
           WRITE(IPT,*) "! ==== MAP NOT ASSOCIATED! ===="
           RETURN
        END IF


        DO I =1,size(MP)
           WRITE(IPT,*) "! *** proc map#",I
           WRITE(IPT,*) "! SIZES: G,N,L",MP(I)%GSIZE,MP(I)%NSIZE,MP(I)%LSIZE
           NN = -1 
           LL = -1
           IF(ASSOCIATED(MP(I)%LOC_2_GL)) NN = ubound(MP(I)%LOC_2_GL,1)
           IF(ASSOCIATED(MP(I)%LOC_2_GRID)) LL = ubound(MP(I)%LOC_2_GRID,1)

           WRITE(IPT,*) "! LOC_2_GL Ubound=",NN
           WRITE(IPT,*) "! LOC_2_GRID Ubound=",LL

        END DO

           WRITE(IPT,*) "! ==== END PRINT_MAP ===="

     END IF

   END SUBROUTINE PRINT_MAP

# else

   INTERFACE PPRINT
      MODULE PROCEDURE PPRINT_ARR
      MODULE PROCEDURE PPRINT_VEC
   END INTERFACE

   INTERFACE APRINT
      MODULE PROCEDURE APRINT_ARR
      MODULE PROCEDURE APRINT_VEC
   END INTERFACE
!===================================================================================|
   CONTAINS   !!INCLUDED SUBROUTINES FOLLOW
!===================================================================================|
# endif

!==============================================================================|
!  WRITE OUT VARIABLE INFORMATION TO LOCAL FILES                               |
!									       |
!  USAGE EXAMPLES                                                              |
!									       |
!  write u velocity at surface in triangle 256 to file fort.306 with iteration |
!  I1 = LBOUND(U,1) ; I2 = UBOUND(U,1)                                         |
!  CALL PPRINT(306,I1,I2,KB,U,"element",256,1,1,FLOAT(IINT))                   |
!									       |
!  I1 = LBOUND(EL,1) ; I2 = UBOUND(EL,1)                                       |
!  write surface elevation at node 233 to file fort.409 with time in hours     |
!  CALL PPRINT(406,I1,I2,1,EL,"node",233,1,1,THOUR)                            |
!									       |
!  I1 = LBOUND(T1,1) ; I2 = UBOUND(T1,1)                                       |
!  write vertical distribution of salinity at node 422 to file fort.433        |
!  CALL PPRINT(433,I1,I2,KB,T1,"node",422,1,KBM1,THOUR)                        |
!									       |
!  ARGUMENT LIST                                                               |
!       PPRINT(IUNIT,LB1,UB1,UB2,VARP,VART,ILOC,K1,K2,REF)                     |
!    1.) IUNIT - UNIT NUMBER FOR OUTPUT FILE (MUST BE >= 300 .and. <7000)      |
!    2.) LB1   - LOWER BOUND OF 1ST ARGUMENT OF ARRAY TO PRINT (USUAlLY 0)     |
!    3.) LB2   - UPPER BOUND OF 1ST ARRAY DIMENSION (USUALLY NT OR MT)         |
!        NOTE: LB1/LB2 CAN BE DETERMINE AUTOMATICALLY WITH LBOUND/UBOUND       |
!    4.) UB2   - UPPER BOUND OF SECOND ARRAY DIMENSION                         | 
!        UB2   = 1 FOR SURFACE ARRAYS LIKE EL,UA                               |
!        UB2   = KB FOR 3D ARRAYS LIKE U/V                                     |
!    5.) VARP  = VARIABLE TO PRINT (ARRAY NAME = U,V,WW,EL,T1,RHO1, etc)       |
!    6.) VART  = VARIABLE LOCATION ("element" or "node")                       |
!    7.) ILOC  = INDEX OF ELEMENT/NODE TO PRINT                                | 
!    8.) K1    = LOWER RANGE OF SIGMA LEVEL TO PRINT                           |
!    9.) K2    = UPPER RANGE OF SIGMA LEVEL TO PRINT                           |
!        K1 = 1,K2 = 1 FOR SURFACE VALUES ONLY                                 |
!        K1 = 1,K2 = KBM1 FOR ALL LEVELS                                       |
!   10.) REF   = REFERENCE VALUE FOR DATA (MUST BE FLOAT)                      |
!        REF = THOURS FOR CALCULATION TIME IN HOURS                            |
!        REF = FLOAT(IINT) FOR ITERATION NUMBER                                |
!   11.) IPT = UNIT TO WRITE ERRORS TO (USE IPT)                               |
!==============================================================================|

!==============================================================================|
   SUBROUTINE APRINT_VEC(IUNIT,VARP,VART,NOW,ILOC,MSG)
     USE MOD_TIME
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) :: NOW
     INTEGER,  INTENT(IN) :: IUNIT,ILOC
     REAL(SP), ALLOCATABLE, INTENT(IN),TARGET :: VARP(:)
     CHARACTER(LEN=*), INTENT(IN) :: VART
     CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: MSG
     REAL(SP), POINTER:: VARP_O(:)
     
     VARP_O =>  VARP
     
     IF(PRESENT(MSG))THEN
        CALL PPRINT_VEC(IUNIT,VARP_O,VART,NOW,ILOC,MSG)
     ELSE
        CALL PPRINT_VEC(IUNIT,VARP_O,VART,NOW,ILOC)
     END IF

   END SUBROUTINE APRINT_VEC
   
!==============================================================================|
   SUBROUTINE PPRINT_VEC(IUNIT,VARP,VART,NOW,ILOC,MSG)
     USE CONTROL
     USE LIMS

     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) :: NOW
     INTEGER,  INTENT(IN) :: IUNIT,ILOC
     REAL(SP), POINTER, INTENT(IN) :: VARP(:)
     CHARACTER(LEN=*), INTENT(IN) :: VART
     CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: MSG

     CHARACTER(LEN=80),parameter :: VAR_E = "element"
     CHARACTER(LEN=80),parameter :: VAR_N = "node"

     CHARACTER(LEN=100) :: STRNG
     CHARACTER(LEN=20) :: short

     INTEGER :: I,J,K,PROCMAX,II,IBND,Kopt,IERR
     LOGICAL :: PRINT_PROC

     !==============================================================================|

     !------------------------------------------------------------------------------|
     !  Process Iunit for Errors                                                    |
     !------------------------------------------------------------------------------|
     IF(IUNIT /= IPT .and. (IUNIT < 300 .or. IUNIT > 7000) )THEN
        CALL FATAL_ERROR('ERROR IN PPRINT',&
             & 'FILE UNIT < 300 AND UNIT > 7000 ARE RESERVED FOR FVCOM I/O',&
             & 'PLEASE INCREASE IUNIT TO 300+')
     END IF

     !------------------------------------------------------------------------------|
     !  Process Vartype for Errors                                                  |
     !------------------------------------------------------------------------------|
     IF(VART /= VAR_E .AND. VART /= VAR_N)THEN
        CALL FATAL_ERROR('VART IN PPRINT NOT CORRECT :'//TRIM(VART),&
             & 'SHOULD BE "'//trim(var_e)//'" or "'//trim(var_n)//'"')
     END IF


     !------------------------------------------------------------------------------|
     !  Process string output                                                       |
     !------------------------------------------------------------------------------|
     IF(PRESENT(MSG)) STRNG=TRIM(MSG)//"; IINT"


     IF(abs(IINT) .lt. 1000) THEN
        WRITE(SHORT,'(I5)') IINT
     ELSE IF(abs(IINT) .lt. 1000000) THEN
        WRITE(SHORT,'(I8)') IINT
     ELSE
        WRITE(SHORT,*) IINT
     END IF

     IF(USE_REAL_WORLD_TIME) THEN
        STRNG = TRIM(STRNG)//TRIM(SHORT)//", Date/Time:"&
             &//TRIM(WRITE_DATETIME(NOW,3,TIMEZONE))//"; ILOC= "
     ELSE
        STRNG = TRIM(STRNG)//TRIM(SHORT)//", Time(s):"
        WRITE(SHORT,'(f16.8)') SECONDS(NOW)
        STRNG = TRIM(STRNG)//TRIM(SHORT)//"; ILOC="
     END IF

     WRITE(SHORT,'(I8)') ILOC
     STRNG = TRIM(STRNG)//TRIM(SHORT)//"; VALUE="

     !------------------------------------------------------------------------------|
     !  Single Processor Case                                                       |
     !------------------------------------------------------------------------------|
     IF(NPROCS == 1)THEN
        WRITE(IUNIT,*) TRIM(STRNG),VARP(ILOC)
     END IF

     !------------------------------------------------------------------------------|
     !  Multi Processor Case with Element Based Variable (u,v,ww, etc)              |
     !      Transform to Local Element ID with "ELID"                               |
     !------------------------------------------------------------------------------|

     IF(NPROCS /= 1 .AND. VART == var_e .AND. ELID(ILOC) /= 0)THEN

        WRITE(IUNIT,*) TRIM(STRNG),VARP(ELID(ILOC))

     END IF

     !------------------------------------------------------------------------------|
     !  Multi Processor Case with Node Based Variable (s1,t1,rho1,e1, etc)          |
     !      Transform to Local Node ID with "NLID"                                  |
     !      If Node is Interprocessor Boundary Node, Choose Processor with Highest  |
     !      ID Number to Write Values to File                                       |
     !------------------------------------------------------------------------------|

     IF(NPROCS /= 1 .AND. VART == var_n .AND. NLID(ILOC) > 0)THEN

        PRINT_PROC = .TRUE.
        IF(NDE_ID(NLID(ILOC)) == 1)THEN   !!BOUNDARY NODE

           DO II=1,NBN
              IF(BN_LST(II) == ILOC) IBND = II
           END DO

           PROCMAX = 10000
           DO J=1,NPROCS
              IF(BN_NEY(IBND,J)==1) THEN
                 IF(J < PROCMAX) PROCMAX = J
              END IF
           END DO

           IF(PROCMAX /=  MYID) PRINT_PROC = .FALSE.  !!NOT RESPONSIBLE FOR OUTPUT 
        END IF

        IF(PRINT_PROC)THEN
           WRITE(IUNIT,*) TRIM(STRNG),VARP(NLID(ILOC))       
        END IF

     END IF


# if defined (MULTIPROCESSOR)
     IF(.NOT. IOPROC) CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
# endif

     RETURN
   END SUBROUTINE PPRINT_VEC
!==============================================================================|
!==============================================================================|
   SUBROUTINE APRINT_ARR(IUNIT,VARP,VART,NOW,ILOC,K1,K2,MSG)
     USE MOD_TIME
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) :: NOW
     INTEGER,  INTENT(IN) :: IUNIT,ILOC,K1
     INTEGER,  INTENT(IN),OPTIONAL :: K2
     REAL(SP), ALLOCATABLE, INTENT(IN),TARGET :: VARP(:,:)
     CHARACTER(LEN=*), INTENT(IN) :: VART
     CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: MSG
     REAL(SP), POINTER:: VARP_O(:,:)
     
     VARP_O =>  VARP
     
     IF (PRESENT(K2)) THEN
        IF(PRESENT(MSG))THEN
           CALL PPRINT_ARR(IUNIT,VARP_O,VART,NOW,ILOC,K1,K2,MSG)
        ELSE
           CALL PPRINT_ARR(IUNIT,VARP_O,VART,NOW,ILOC,K1,K2)
        END IF
     ELSE
        IF(PRESENT(MSG))THEN
           CALL PPRINT_ARR(IUNIT,VARP_O,VART,NOW,ILOC,K1,K1,MSG)
        ELSE
           CALL PPRINT_ARR(IUNIT,VARP_O,VART,NOW,ILOC,K1)
        END IF
     END IF
   END SUBROUTINE APRINT_ARR
   
!==============================================================================|
   SUBROUTINE PPRINT_ARR(IUNIT,VARP,VART,NOW,ILOC,K1,K2,MSG)
     USE CONTROL
     USE LIMS

     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) :: NOW
     INTEGER,  INTENT(IN) :: IUNIT,ILOC,K1
     INTEGER,  INTENT(IN),OPTIONAL :: K2
     REAL(SP), POINTER, INTENT(IN) :: VARP(:,:)
     CHARACTER(LEN=*), INTENT(IN) :: VART
     CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: MSG

     CHARACTER(LEN=80),parameter :: VAR_E = "element"
     CHARACTER(LEN=80),parameter :: VAR_N = "node"

     CHARACTER(LEN=100) :: STRNG
     CHARACTER(LEN=20) :: short

     INTEGER :: I,J,K,PROCMAX,II,IBND,Kopt,IERR
     LOGICAL :: PRINT_PROC

     !==============================================================================|

     !------------------------------------------------------------------------------|
     !  Process Iunit for Errors                                                    |
     !------------------------------------------------------------------------------|
     IF(IUNIT /= IPT .and. (IUNIT < 300 .or. IUNIT > 7000) )THEN
        CALL FATAL_ERROR('ERROR IN PPRINT',&
             & 'FILE UNIT < 300 AND UNIT > 7000 ARE RESERVED FOR FVCOM I/O',&
             & 'PLEASE INCREASE IUNIT TO 300+')
     END IF

     !------------------------------------------------------------------------------|
     !  Process Vartype for Errors                                                  |
     !------------------------------------------------------------------------------|
     IF(VART /= VAR_E .AND. VART /= VAR_N)THEN
        CALL FATAL_ERROR('VART IN PPRINT NOT CORRECT :'//TRIM(VART),&
             & 'SHOULD BE "'//trim(var_e)//'" or "'//trim(var_n)//'"')
     END IF

     !------------------------------------------------------------------------------|
     !  Process optional sigma level range                                          |
     !------------------------------------------------------------------------------|

     IF(PRESENT(K2)) THEN
        KOPT=K2
     ELSE
        KOPT=K1
     END IF

     !------------------------------------------------------------------------------|
     !  Process string output                                                       |
     !------------------------------------------------------------------------------|
     IF(PRESENT(MSG)) STRNG=TRIM(MSG)//"; IINT"


     IF(abs(IINT) .lt. 1000) THEN
        WRITE(SHORT,'(I5)') IINT
     ELSE IF(abs(IINT) .lt. 1000000) THEN
        WRITE(SHORT,'(I8)') IINT
     ELSE
        WRITE(SHORT,*) IINT
     END IF

     IF(USE_REAL_WORLD_TIME) THEN
        STRNG = TRIM(STRNG)//TRIM(SHORT)//", Date/Time:"&
             &//TRIM(WRITE_DATETIME(NOW,3,TIMEZONE))//"; ILOC "
     ELSE
        STRNG = TRIM(STRNG)//TRIM(SHORT)//", Time(s):"
        WRITE(SHORT,'(f16.8)') SECONDS(NOW)
        STRNG = TRIM(STRNG)//TRIM(SHORT)//"; ILOC="
     END IF

     WRITE(SHORT,'(I8)') ILOC
     STRNG = TRIM(STRNG)//TRIM(SHORT)//"; VALUES="
     !------------------------------------------------------------------------------|
     !  Single Processor Case                                                       |
     !------------------------------------------------------------------------------|
     IF(NPROCS == 1)THEN
        WRITE(IUNIT,*) TRIM(STRNG)
        WRITE(IUNIT,*) (VARP(ILOC,K),K=K1,KOPT)
     END IF

     !------------------------------------------------------------------------------|
     !  Multi Processor Case with Element Based Variable (u,v,ww, etc)              |
     !      Transform to Local Element ID with "ELID"                               |
     !------------------------------------------------------------------------------|

     IF(NPROCS /= 1 .AND. VART == var_e .AND. ELID(ILOC) /= 0)THEN

        WRITE(IUNIT,*) TRIM(STRNG)
        WRITE(IUNIT,*) (VARP(ELID(ILOC),K),K=K1,KOPT)

     END IF

     !------------------------------------------------------------------------------|
     !  Multi Processor Case with Node Based Variable (s1,t1,rho1,e1, etc)          |
     !      Transform to Local Node ID with "NLID"                                  |
     !      If Node is Interprocessor Boundary Node, Choose Processor with Highest  |
     !      ID Number to Write Values to File                                       |
     !------------------------------------------------------------------------------|

     IF(NPROCS /= 1 .AND. VART == var_n .AND. NLID(ILOC) > 0)THEN

        PRINT_PROC = .TRUE.
        IF(NDE_ID(NLID(ILOC)) == 1)THEN   !!BOUNDARY NODE

           DO II=1,NBN
              IF(BN_LST(II) == ILOC) IBND = II
           END DO

           PROCMAX = 10000
           DO J=1,NPROCS
              IF(BN_NEY(IBND,J)==1) THEN
                 IF(J < PROCMAX) PROCMAX = J
              END IF
           END DO

           IF(PROCMAX /=  MYID) PRINT_PROC = .FALSE.  !!NOT RESPONSIBLE FOR OUTPUT 
        END IF

        IF(PRINT_PROC)THEN
           WRITE(IUNIT,*) TRIM(STRNG)
           WRITE(IUNIT,*) (VARP(NLID(ILOC),K),K=K1,KOPT)        
        END IF

     END IF

# if defined (MULTIPROCESSOR)
     IF(.NOT. IOPROC) CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
# endif

     RETURN
   END SUBROUTINE PPRINT_ARR
!==============================================================================|

END MODULE MOD_PAR  
# endif
