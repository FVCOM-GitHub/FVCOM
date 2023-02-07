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

!==============================================================
SUBROUTINE SETUP_DOMAIN
!==============================================================
! IF(PAR): 
!          Uses METIS to determine local domains
!          Create global to local and local to global maps
!==============================================================
  USE MOD_PAR
  USE ALL_VARS
  USE MOD_NCDIO
  USE MOD_SETUP
  
  !==============================================================================|
  !  GENERATE LOCAL NODE CONNECTIVITY (NV) FROM GLOBAL NODE CONNECTIVITY (NVG)   |
  !  USING LOCAL TO GLOBAL MAPPING FOR INTERIOR ELEMENTS (EGID)                  |
  !  AND LOCAL TO GLOBAL MAPPING FOR HALO ELEMENTS (HE_LST)                      |
  !  TRANSFER GLOBAL COORDINATES TO LOCAL DOMAINS
  !==============================================================================|

  IMPLICIT NONE
  INTEGER :: SENDER,STATUS
  INTEGER :: nvals, IERR, I, EGL
  

     if (dbg_set(dbg_log)) &
          & write(IPT,*) "! SETTING UP FVCOM PROCESSOR ENVIRONMENT"
  
  ! SETUP OUTPUT MODE OR TURN IT OFF
  CALL SETUP_MPI_IO_MODE(USE_MPI_IO_MODE,MPI_FVCOM_GROUP)
  
  IF(SERIAL) THEN
     
     ! DO NOT ALLOCATE EL_PID FOR NON-MULTIPROCESSOR CASES
          
     ! SHIFT GRID DIMENSIONS TO LOCAL FOR ANY SERIAL CASE
     MT=MGL
     NT=NGL
     
     M = MGL
     N = NGL
     
     ! TRANSFER CONNECTIVTY TO LOCAL DOMAIN NAME VARIABLES
     ALLOCATE(NV(0:NT,4))
     NV   = NVG
     NV(0,:) = 0
     !DEALLOCATE(NVG)
     

     ! MUST ALLOCATE EL_PID FOR OUTPUT OF DATA
     ALLOCATE(EL_PID(NGL))  ;  EL_PID = 1


     ! SETUP NGID
     ALLOCATE(NGID(0:MGL))  ; 
     DO i=0,MGL
        NGID(i)=i
     END DO

     NGID_X => NGID

     ! SETUP NLID
     ALLOCATE(NLID(0:MGL))
     DO i=0,MGL
        NLID(i) = i
     END DO

     NLID_X => NLID

     ! SETUP EGID
     ALLOCATE(EGID(0:NGL))  ; 
     DO i=0,NGL
        EGID(i)=i
     END DO

     EGID_X => EGID

     ! SETUP ELID
     ALLOCATE(ELID(0:NGL))
     DO i=0,NGL
        ELID(i) = i
     END DO

     ELID_X => ELID

     ALLOCATE(NDE_ID(MGL)); NDE_ID=0

  ELSE ! IF PAR

#  if defined (MULTIPROCESSOR) 
     
     !
     !  DECOMPOSE DOMAIN BY ELEMENTS USING METIS
     !
     ALLOCATE(EL_PID(NGL))  ;  EL_PID = 1
     CALL DOMDEC(NGL,NVG,NPROCS,EL_PID,MSR)
     
     ! GENERATE MAP TO EXCHANGE DATA BETWEEN GLOBAL AND LOCAL DOMAINS
     CALL GENMAP

     IF(.NOT. IOPROC) THEN
        ! TRANSFER CONNECTIVTY TO LOCAL DOMAIN NAME VARIABLES
        ALLOCATE(NV(0:NT,4))
        ! DUPLICATE FOR NC OUTPUT - HAS GLOBAL NODE NUMBERING
!        ALLOCATE(NVGL(0:NT,3))
        
        DO I=1,NT
           !(MUST HAVE GLOBAL NODE NUMBERING)
           NV(I,1:4) = NLID_X(NVG(EGID_X(I),1:4))! NV is the local connectivity
        END DO
        NV(0,:) = 0

     END IF

# endif
  END IF
  

  ! PASS OBC TO LOCAL VARIABLES AND GENERATE MAP FROM GLOBAL

  IF(OBC_LONGSHORE_FLOW_ON) CALL GENMAP_LSF

  IF(OBC_ON) CALL GENMAP_OBC


  IF (USE_MPI_IO_MODE) then
# if defined (MULTIPROCESSOR)
     CALL MPI_IO_LOOP
# else
     CALL FATAL_ERROR("'USE_MPI_IO_MODE' SHOULD NEVER BE TRUE WHEN C&
          &OMPILED WITHOUT MULTIPROCESSOR")
# endif
  END IF
  
  if(DBG_SET(dbg_log))  write(IPT,*) "! Finished SETUP_LOCAL_DOMAINS"
  

END SUBROUTINE SETUP_DOMAIN

