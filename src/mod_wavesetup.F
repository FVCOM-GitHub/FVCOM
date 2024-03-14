!--------------------------------------------------------------------|
!--------------------------------------------------------------------|
  MODULE MOD_WAVESETUP
# if defined (WAVE_SETUP)

  IMPLICIT NONE

! must include all necessary headers in petsc fortran interface
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  PRIVATE
  PUBLIC :: WAVE_INDUCED_FORCE
  PUBLIC :: WAVE_INDUCED_SETUP
  PUBLIC :: ALLOC_VARS_WSU

  REAL, PUBLIC, ALLOCATABLE  ::  RHS(:)
  INTEGER :: IDXMAX
  REAL, ALLOCATABLE,PUBLIC :: SETUP2(:),SETUP2_TMP(:)

  CONTAINS
!==========================================================================|
  SUBROUTINE WAVE_INDUCED_FORCE(DEP2,DEPSAV)
!  USE ALL_VARS
# if defined (MULTIPROCESSOR)
  USE MOD_PAR
# endif
  USE MOD_USGRID    !, ONLY : MDC,MSC
  
  IMPLICIT NONE
  
  REAL :: DEP2(MT),DEPSAV(MT)
  REAL :: SIGLOC,KLOC,CGLOC
  
  REAL :: DEPMAX,RSXX,RSXY,RSYY,SXXIJ,SXYIJ,SYYIJ,FXIJ,FYIJ
  REAL :: DEPLOC,CK,ELOC,RHS_TMP
  REAL, ALLOCATABLE :: PSXXPX(:),PSXYPX(:),PSXYPY(:),PSYYPY(:)
  REAL, ALLOCATABLE :: FX(:),FY(:)
  REAL, ALLOCATABLE :: SETPDA(:,:)
  INTEGER :: I,ISS,ID,IA,IB,J1,J2,I1,J
  REAL    :: NR,ND
  REAL    :: DLTXE_TMP,DLTYE_TMP
  REAL    :: DEPSAV_TMP(MGL)
  INTEGER :: IERR
  
  DEPMAX = 0.
  IDXMAX = 0
!
!---------initializing SETPDA array
  ALLOCATE(SETPDA(MT,5));    SETPDA = 0.

# if defined(MULTIPROCESSOR)
  IF(PAR) CALL EXCHANGE(NC,MT,1,MYID,NPROCS,DEP2)
# endif  

  IF(SERIAL)THEN
    DO I=1,M
!
!---------seek deepest point
!
      IF(DEPSAV(I) > DEPMAX)THEN
        DEPMAX = DEPSAV(I)
        IDXMAX = I
      END IF
    END DO  
  END IF  

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    CALL GATHER(LBOUND(DEPSAV,1),UBOUND(DEPSAV,1),M,MGL,1,MYID,NPROCS,   &
                NMAP,DEPSAV,DEPSAV_TMP)    
    IF(MSR)THEN
      DO I=1,MGL
!
!---------seek deepest point
!
        IF(DEPSAV_TMP(I) > DEPMAX)THEN
          DEPMAX = DEPSAV_TMP(I)
          IDXMAX = I
        END IF
      END DO  
    END IF
    CALL MPI_BCAST(DEPMAX,1,MPI_F,0,MPI_FVCOM_GROUP,IERR)      
    CALL MPI_BCAST(IDXMAX,1,MPI_F,0,MPI_FVCOM_GROUP,IERR)      
  END IF  
# endif
  
  DO I = 1,MT
    IF(DEP2(I) > DEPMIN)THEN
!JQI!
!JQI!---------seek deepest point
!JQI!
!JQI      IF(DEPSAV(I) > DEPMAX)THEN
!JQI        DEPMAX = DEPSAV(I)
!JQI	IDXMAX = I
!JQI      END IF
!
!---------compute radiation stress components RSXX,RSXY and RSYY 
!
      RSXX = 0.
      RSXY = 0.
      RSYY = 0.
      DEPLOC = DEP2(I)
      DO ISS = 1,MSC
        SIGLOC = SPCSIG(ISS)
	CALL KSCIP1(1,SIGLOC,DEPLOC,KLOC,CGLOC,NR,ND)
	CK = CGLOC*KLOC
	DO ID = 1, MDC
	  ELOC = SIGLOC*AC2(ID,ISS,I)
!                          -
!                          / (cos(Theta))^2         for i = 4
!         SPCDIR(ID,i) is </  sin(Theta)cos(Theta)  for i = 5
!                          / (sin(Theta))^2         for i = 6
!                          -
          RSXX = RSXX+(CK*SPCDIR(ID,4)+CK-SIGLOC/2.)*ELOC
	  RSXY = RSXY+CK*SPCDIR(ID,5)*ELOC
	  RSYY = RSYY+(CK*SPCDIR(ID,6)+CK-SIGLOC/2.)*ELOC
	END DO
      END DO
!
!---------store radiation stress components in array SETPDA
!
!     DDIR   is width of directional band
!     FRINTF is frequency integration factor df/f
!
      SETPDA(I,3) = RSXX*DDIR*FRINTF
      SETPDA(I,4) = RSXY*DDIR*FRINTF
      SETPDA(I,5) = RSYY*DDIR*FRINTF
    END IF
  END DO
      print*, depmax,idxmax,'CHECKCHECKCHECK'
  
  ALLOCATE(PSXXPX(0:NT));     PSXXPX = 0.0
  ALLOCATE(PSXYPX(0:NT));     PSXYPX = 0.0
  ALLOCATE(PSXYPY(0:NT));     PSXYPY = 0.0
  ALLOCATE(PSYYPY(0:NT));     PSYYPY = 0.0
  ALLOCATE(FX(0:NT));           FX     = 0.0
  ALLOCATE(FY(0:NT));           FY     = 0.0
  
  DO I = 1, NE
    IA = IEC(I,1)
    IB = IEC(I,2)
    J1 = IENODE(I,1)
    J2 = IENODE(I,2)
    SXXIJ = 0.5*(SETPDA(J1,3)+SETPDA(J2,3))
    SXYIJ = 0.5*(SETPDA(J1,4)+SETPDA(J2,4))
    SYYIJ = 0.5*(SETPDA(J1,5)+SETPDA(J2,5))
    
    PSXXPX(IA) = PSXXPX(IA)+SXXIJ*DLTYC(I)
    PSXYPX(IA) = PSXYPX(IA)+SXYIJ*DLTYC(I)
    PSXYPY(IA) = PSXYPY(IA)-SXYIJ*DLTXC(I)
    PSYYPY(IA) = PSYYPY(IA)-SYYIJ*DLTXC(I)
    
    PSXXPX(IB) = PSXXPX(IB)-SXXIJ*DLTYC(I)
    PSXYPX(IB) = PSXYPX(IB)-SXYIJ*DLTYC(I)
    PSXYPY(IB) = PSXYPY(IB)+SXYIJ*DLTXC(I)
    PSYYPY(IB) = PSYYPY(IB)+SYYIJ*DLTXC(I)
  END DO
 
  DO I = 1, NT
    PSXXPX(I) = PSXXPX(I)/ART(I)
    PSXYPX(I) = PSXYPX(I)/ART(I)
    PSXYPY(I) = PSXYPY(I)/ART(I)
    PSYYPY(I) = PSYYPY(I)/ART(I)
    
    FX(I) = -PSXXPX(I)-PSXYPY(I)
    FY(I) = -PSXYPX(I)-PSYYPY(I)
  END DO    
    
  IF(.NOT. ALLOCATED(RHS)) ALLOCATE(RHS(0:MT))
  RHS  = 0.0
  
  DO I = 1, NCV
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)
    
    FXIJ = FX(I1)
    FYIJ = FY(I1)
    RHS_TMP = -FXIJ*DLTYE(I)+FYIJ*DLTXE(I)
    RHS(IA) = RHS(IA)-RHS_TMP
    RHS(IB) = RHS(IB)+RHS_TMP
  END DO
  
  DO I = 1, M
    IF(ISONB(I) /= 0)THEN
      FXIJ = 0.0
      FYIJ = 0.0
      DO J=1,NTVE(I)
        FXIJ = FXIJ + FX(NBVE(I,J))
	FYIJ = FYIJ + FY(NBVE(I,J))
      END DO
      FXIJ = FXIJ/NTVE(I)
      FYIJ = FYIJ/NTVE(I)

      IF(NBSN(I,NTSN(I)-1) > M)THEN
        DLTXE_TMP = VX(I)-VX(NBSN(I,2))
        DLTYE_TMP = VY(I)-VY(NBSN(I,2))
      ELSE IF(NBSN(I,2) > M)THEN
        DLTXE_TMP = VX(NBSN(I,NTSN(I)-1))-VX(I)
        DLTYE_TMP = VY(NBSN(I,NTSN(I)-1))-VY(I)
      ELSE 
        DLTXE_TMP = VX(NBSN(I,NTSN(I)-1))-VX(NBSN(I,2))
        DLTYE_TMP = VY(NBSN(I,NTSN(I)-1))-VY(NBSN(I,2))
      END IF

      RHS_TMP = -FXIJ*DLTYE_TMP+FYIJ*DLTXE_TMP
      RHS(I) = RHS(I)-RHS_TMP*0.5
    END IF
  END DO
  
  DO I = 1, M    !T
    RHS(I) = RHS(I)/ART1(I)
  END DO
      
  DEALLOCATE(PSXXPX,PSXYPX,PSXYPY,PSYYPY)
  DEALLOCATE(FX,FY)
  DEALLOCATE(SETPDA)

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
!JQI    CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,RHS)
    CALL EXCHANGE(NC,MT,1,MYID,NPROCS,RHS)
  END IF    
# endif

  RETURN
  END SUBROUTINE WAVE_INDUCED_FORCE
!=========================================================================|
!
!=========================================================================|  
  SUBROUTINE WAVE_INDUCED_SETUP(DEP2,DEPSAV)
  USE MOD_USGRID
!  USE MOD_OBCS
  USE ALL_VARS
# if defined (MULTIPROCESSOR)
  USE MOD_PAR
# endif
  USE MOD_PETSC, ONLY : ALO_2_PLO_NODE,BL_EL,B_EL,N_VERTS,L2G_EL,A_EL,  &
                        PLO_2_PGO_NODE,X_EL,XL_EL,G2L_EL,XVALS_EL,      &
			PLO_2_ALO_NODE,PETSc_SETICS_EL,PETSc_SOLVER_EL
  
  IMPLICIT NONE
  
  REAL :: DEP2(MT),DEPSAV(MT)
  
  INTEGER :: I,IA,IB,J,I1,LN,NNZ2,LOC,JJ,IK,NODE
  INTEGER :: PROW1,PROW2,PCOL1,PCOL2,PCOL3,TMP2,TMP3,PETSc_POS
  REAL    :: AA1,AB1,AA2,AB2,AA3,AB3,DIJ
  REAL    :: S_UPDP,S_UPCOR
  REAL    :: DLTXE_TMP,DLTYE_TMP
  REAL    :: SETUP2_TEMP(MGL)
!  REAL, ALLOCATABLE :: PCOEF(:,:),PCOEF1(:,:)
  
  PetscReal :: STERM
  PetscReal :: VCOL1,VCOL2
  PetscInt:: IERR

  CALL VecSet(BL_EL,ZERO,IERR);CHKERRQ(IERR)
  CALL VecSet(B_EL,ZERO,IERR);CHKERRQ(IERR)

  DO I=1,M
    PETSc_POS = ALO_2_PLO_NODE(I)
    IF(PETSc_POS > N_VERTS) CYCLE
    
    IF(ISONB(I) /= 2)THEN
      STERM = -RHS(I)
    ELSE
      STERM = 0.0
    END IF    
    CALL VecSetValues(BL_EL,1,PETSc_POS-1,STERM,INSERT_VALUES,IERR);CHKERRQ(IERR)
  END DO
  
# if defined (OLD_PETSC)
  CALL VecScatterBegin(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(L2G_EL,BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(L2G_EL,BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif


  CALL VecAssemblyBegin(B_EL,IERR);CHKERRQ(IERR)
  CALL VecAssemblyEnd(B_EL,IERR);CHKERRQ(IERR)

  CALL MatZeroEntries(A_EL,IERR);CHKERRQ(IERR)   
    
# if defined(MULTIPROCESSOR)
!  IF(PAR) CALL EXCHANGE(NC,MT,1,MYID,NPROCS,DEP2)
!  IF(PAR) CALL EXCHANGE(NC,MT,1,MYID,NPROCS,DEPSAV)
# endif  

  DO I = 1, NCV_I
    I1 = NTRG(I)
    IA = NIEC(I,1)
    IB = NIEC(I,2)
    DIJ = DEP2(NV(I1,1))+DEP2(NV(I1,2))+DEP2(NV(I1,3))
    DIJ = DIJ/3.0
    
    IF(ISONB(IA)+ISONB(IB) < 4)THEN
      DO J = 1, 3
        IF(NV(I1,J) == IA)THEN
	  TMP2 = NV(I1,J+1-INT((J+1)/4)*3)
	  TMP3 = NV(I1,J+2-INT((J+2)/4)*3)
	  AA1  = -0.5*(VX(TMP3)-VX(TMP2))/ART(I1)
          AB1  =  0.5*(VY(TMP3)-VY(TMP2))/ART(I1)
          AA2  = -0.5*(VX(IA)  -VX(TMP3))/ART(I1)
          AB2  =  0.5*(VY(IA)  -VY(TMP3))/ART(I1)
          AA3  = -0.5*(VX(TMP2)-VX(IA)  )/ART(I1)
          AB3  =  0.5*(VY(TMP2)-VY(IA)  )/ART(I1)
        END IF
      END DO
     
      PROW1 = ALO_2_PLO_NODE(IA)
      PROW2 = ALO_2_PLO_NODE(IB)
      
      PCOL1 = ALO_2_PLO_NODE(IA)
      PCOL2 = ALO_2_PLO_NODE(IB)
      PCOL3 = ALO_2_PLO_NODE(TMP3)
      
      VCOL1 = (-AA1*DLTXE(I)+AB1*DLTYE(I))*DIJ/ART1(IA)
      VCOL2 = -(-AA1*DLTXE(I)+AB1*DLTYE(I))*DIJ/ART1(IB)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)	  
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL1-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
      
      VCOL1 = (-AA2*DLTXE(I)+AB2*DLTYE(I))*DIJ/ART1(IA)
      VCOL2 = -(-AA2*DLTXE(I)+AB2*DLTYE(I))*DIJ/ART1(IB)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL2-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)	  
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL2-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
      
      VCOL1 = (-AA3*DLTXE(I)+AB3*DLTYE(I))*DIJ/ART1(IA)
      VCOL2 = -(-AA3*DLTXE(I)+AB3*DLTYE(I))*DIJ/ART1(IB)
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL3-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)	  
      CALL MatSetValuesLocal(A_EL,1,PROW2-1,1,PCOL3-1,VCOL2,ADD_VALUES,IERR);CHKERRQ(IERR)
    END IF
  END DO

  DO I=1,M
    IF(ISONB(I) == 1)THEN
      DIJ = DEP2(I)
      DO J=1,NTVE(I)
        DO JJ = 1, 3
          IF(NV(NBVE(I,J),JJ) == I)THEN
	    TMP2 = NV(NBVE(I,J),JJ+1-INT((JJ+1)/4)*3)
	    TMP3 = NV(NBVE(I,J),JJ+2-INT((JJ+2)/4)*3)
	    AA1  = -0.5*(VX(TMP3)-VX(TMP2))/ART(NBVE(I,J))/NTVE(I)
            AB1  =  0.5*(VY(TMP3)-VY(TMP2))/ART(NBVE(I,J))/NTVE(I)
            AA2  = -0.5*(VX(I)  -VX(TMP3))/ART(NBVE(I,J))/NTVE(I)
            AB2  =  0.5*(VY(I)  -VY(TMP3))/ART(NBVE(I,J))/NTVE(I)
            AA3  = -0.5*(VX(TMP2)-VX(I)  )/ART(NBVE(I,J))/NTVE(I)
            AB3  =  0.5*(VY(TMP2)-VY(I)  )/ART(NBVE(I,J))/NTVE(I)
          END IF
        END DO
     
        PROW1 = ALO_2_PLO_NODE(I)
      
        PCOL1 = ALO_2_PLO_NODE(I)
        PCOL2 = ALO_2_PLO_NODE(TMP2)
        PCOL3 = ALO_2_PLO_NODE(TMP3)
      
        IF(NBSN(I,NTSN(I)-1) > M)THEN
          DLTXE_TMP = VX(I)-VX(NBSN(I,2))
          DLTYE_TMP = VY(I)-VY(NBSN(I,2))
        ELSE IF(NBSN(I,2) > M)THEN
          DLTXE_TMP = VX(NBSN(I,NTSN(I)-1))-VX(I)
          DLTYE_TMP = VY(NBSN(I,NTSN(I)-1))-VY(I)
        ELSE 
          DLTXE_TMP = VX(NBSN(I,NTSN(I)-1))-VX(NBSN(I,2))
          DLTYE_TMP = VY(NBSN(I,NTSN(I)-1))-VY(NBSN(I,2))
        END IF
      
        VCOL1 = (-AA1*DLTXE_TMP+AB1*DLTYE_TMP)*DIJ*0.5/ART1(I)
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)	  
      
        VCOL1 = (-AA2*DLTXE_TMP+AB2*DLTYE_TMP)*DIJ*0.5/ART1(I)
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL2-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)	  
      
        VCOL1 = (-AA3*DLTXE_TMP+AB3*DLTYE_TMP)*DIJ*0.5/ART1(I)
        CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL3-1,VCOL1,ADD_VALUES,IERR);CHKERRQ(IERR)	  
      END DO
      !!!      
    END IF
  END DO

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)

  DO I=1, IOBCN

    PROW1 = ALO_2_PLO_NODE(I_OBC_N(I))
    IF(PROW1>N_VERTS) CYCLE

    VCOL1 = 0.0D0
    NODE  = I_OBC_N(I)
    DO J=1, NTSN(NODE)-1

!      PCOL1 = ALO_2_PLO3(NBSN(NODE,J))
      PCOL1 = ALO_2_PLO_NODE(NBSN(NODE,J))
      CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PCOL1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)

    ENDDO

    VCOL1 = 1.0D0
    PROW1 = ALO_2_PLO_NODE(NODE)
    CALL MatSetValuesLocal(A_EL,1,PROW1-1,1,PROW1-1,VCOL1,INSERT_VALUES,IERR);CHKERRQ(IERR)

  ENDDO

  CALL MatAssemblyBegin(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
  CALL MatAssemblyEnd(A_EL,MAT_FINAL_ASSEMBLY,IERR);CHKERRQ(IERR)
!  CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
     print*,'STOP HEREEEEEEEEEEEEEEEEEEEEEE111111111',MYID  

  CALL PETSc_SETICS_EL(SETUP2)
     print*,'STOP HEREEEEEEEEEEEEEEEEEEEEEE222222222',MYID  
  CALL PETSc_SOLVER_EL

# if defined (OLD_PETSC)
  CALL VecScatterBegin(X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,G2L_EL,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,G2L_EL,IERR);CHKERRQ(IERR)
# else
  CALL VecScatterBegin(G2L_EL,X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
  CALL VecScatterEnd(G2L_EL,X_EL,XL_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
# endif


  CALL VecGetArrayF90(XL_EL,XVALS_EL,IERR);CHKERRQ(IERR)


   IF(LSETUP == 1)THEN                                               
!    set set-up to 0 for deepest point (This is allowed because the    
!    solution of a Poisson equation + constant is again a solution of  
!    the same Poisson equation)                                        
     SETUP2 = 0.0
     DO I=1,N_VERTS
       IK = PLO_2_ALO_NODE(I)
       SETUP2(IK) = XVALS_EL(I)
     END DO
       
     print*,'STOP HEREEEEEEEEEEEEEEEEEEEEEE',MYID  
       
     IF(SERIAL)THEN
       S_UPDP = SETUP2(IDXMAX)                                           
       S_UPCOR = S_UPDP - PSETUP(2) 
     END IF   
         
# if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,SETUP2)
       CALL EXCHANGE(NC,MT,1,MYID,NPROCS,SETUP2)
       CALL GATHER(LBOUND(SETUP2,1),UBOUND(SETUP2,1),M,MGL,1,MYID,NPROCS, &
                   NMAP,SETUP2,SETUP2_TEMP)    
       IF(MSR)THEN
         S_UPDP = SETUP2_TEMP(IDXMAX)                                           
         S_UPCOR = S_UPDP - PSETUP(2)  
       END IF	 
       CALL MPI_BCAST(S_UPDP,1,MPI_F,0,MPI_FVCOM_GROUP,IERR)      
       CALL MPI_BCAST(S_UPCOR,1,MPI_F,0,MPI_FVCOM_GROUP,IERR)      
     END IF  
# endif

!JQI     S_UPDP = SETUP2(IDXMAX)                                           
!JQI     S_UPCOR = S_UPDP - PSETUP(2)  
     
     print*,"S_UPDP=",S_UPDP,IDXMAX,S_UPCOR,PSETUP(2)
                                         
     DO I = 1, MT                                                    
       IF(DEP2(I) > DEPMIN)THEN                           
         SETUP2(I) = SETUP2(I) - S_UPCOR  
       ELSE
!           IF(ABS(SETUP2(I)) > 1.E-7)THEN                  
!             CHARS(1) = INTSTR(I)                            
!             CALL TXPBLA(CHARS(1),IF1,IL1)                      
!             MSGSTR = 'Set-up in dry point with index '//            &
!                        CHARS(1)(IF1:IL1)                         
!             CALL MSGERR ( 2, MSGSTR )                          
!           END IF                                                
       END IF                                                   
     END DO                                                         
   END IF                                                              

   SETUP2_TMP = 0.0
!   DO I=1, MCGRD
!     write(401,*) i,ndd(i),SETUP2_TMP(NDD(I)),SETUP2(I)
!     IF(NDD(I) /= 0)       &
!     SETUP2_TMP(NDD(I)) = SETUP2(I)
!   ENDDO
!
!  --- include computed set-up to depth
!
   DO I = 1, MT
     DEP2(I) = DEPSAV(I) + SETUP2(I)
!     write(100,*)i,setup2(i),depsav(i),dep2(i),s_upcor
   END DO

# if defined (MULTIPROCESSOR)
    IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,DEP2)
    IF(PAR) CALL EXCHANGE(NC,MT,1,MYID,NPROCS,DEP2)
!    IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,SETUP2)
!    IF(PAR) CALL EXCHANGE(NC,MT,1,MYID,NPROCS,SETUP2)
# endif

  RETURN
  END SUBROUTINE WAVE_INDUCED_SETUP
!==========================================================================
!
!========================================================================|
!
  SUBROUTINE ALLOC_VARS_WSU
  USE ALL_VARS
  
  IMPLICIT NONE
  
  ALLOCATE(SETUP2_TMP(0:MT));     SETUP2_TMP = 0.0
  ALLOCATE(SETUP2(0:MT));         SETUP2     = 0.0
  ALLOCATE(RHS(0:MT));            RHS        = 0.0

  RETURN    	  	    			
  END SUBROUTINE ALLOC_VARS_WSU
!==============================================================
# endif
END MODULE MOD_WAVESETUP

