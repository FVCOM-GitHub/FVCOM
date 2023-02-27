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

MODULE MOD_WD

   USE MOD_PREC
   IMPLICIT NONE
   SAVE
!
!--Parameters for Wet/Dry Treatment                 
!

!-----variables controlling porosities through wet/dry determination----------------!
   
   INTEGER , ALLOCATABLE :: ISWETN(:)           !!NODE POROSITY AT NODES FOR TIME N
   INTEGER , ALLOCATABLE :: ISWETC(:)           !!CELL POROSITY AT CELLS FOR TIME N
   INTEGER , ALLOCATABLE :: ISWETNT(:)          !!NODE POROSITY AT NODES FOR TIME N-1 INTERNAL
   INTEGER , ALLOCATABLE :: ISWETCT(:)          !!CELL POROSITY AT CELLS FOR TIME N-1 INTERNAL
   INTEGER , ALLOCATABLE :: ISWETCE(:)          !!CELL POROSITY AT CELLS FOR TIME N-1 EXTERNAL

!qxu{for inundation maps
   INTEGER , ALLOCATABLE :: INUNDATION_MAPS(:)         !!is flooding cell for Time N
!qxu}
!!   REAL(SP), ALLOCATABLE :: UAS(:)          !!VERT-AVGD X-VELOC USED FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: VAS(:)          !!VERT-AVGD Y-VELOC USED FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: UARDS(:)        !!UA AVGD OVER EXT INT FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: VARDS(:)        !!VA AVGD OVER EXT INT FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: US(:,:)         !!X-VELOCITY FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: VS(:,:)         !!Y-VELOCITY FOR MASS CONS IN FLOODING/DRYING PROCESS

   CONTAINS !------------------------------------------------------------------!
            ! ALLOC_WD_DATA       :   ALLOCATE AND INITIALIZE WET/DRY ARRAYS   !
            ! SET_WD_DATA         :   SET VALUES IN WET/DRY ARRAYS             !
            ! WET_JUDGE           :   DETERMINE IF NODES/ELEMENTS ARE WET/DRY  ! 
            ! WD_UPDATE           :   SWAP WET/DRY VARIABLES BETWEEN TIME LEVS ! 
            ! WD_DUMP             :   DUMP WET/DRY FLAGS FOR RESTART           !
            ! WD_READ             :   READ WET/DRY FLAGS FOR RESTART           !
            ! -----------------------------------------------------------------!

!==============================================================================|
!==============================================================================|

!!$   SUBROUTINE SETUP_WETDRY 
!!$!------------------------------------------------------------------------------|
!!$!  READ IN PARAMETERS CONTROLLING WET/DRY TREATMENT                            |
!!$!------------------------------------------------------------------------------|
!!$   IMPLICIT NONE
!!$   INTEGER  :: ISCAN
!!$   CHARACTER(LEN=120) :: FNAME
!!$
!!$   ! NO PARAMETERS NECISARRY FOR WET DRY?
!!$        ONLY MAKE FILE OPTIONS
!!$!------------------------------------------------------------------------------|
!!$!   READ IN VARIABLES AND SET VALUES                                           |
!!$!------------------------------------------------------------------------------|
!!$ END SUBROUTINE SETUP_WETDRY


!==============================================================================|
!==============================================================================|
   SUBROUTINE SET_WD_DATA 
!------------------------------------------------------------------------------|
!  INITIALIZE ARRAYS USED FOR WET/DRY TREATMENT                                |
!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE
   INTEGER :: I

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: SET_WD_DATA"

   IF(STARTUP_TYPE == STARTUP_TYPE_COLDSTART) THEN

!-------- SET WET/DRY FLAGS AND MODIFY WATER SURFACE ELEVATION-----------------!

     CALL WET_JUDGE

!-------- EXCHANGE MODIFIED FREE SURFACE ELEVATION ACROSS PROCESSOR BOUNDS-----!

#  if defined (MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,ELF)
#  endif

!-------- TRANSFER ELEVATION FIELD TO DEPTH AND OLD TIME LEVELS----------------!
#   if defined (ICE_EMBEDDING)
  DO I = 1, NT
    QTHICE1(I) = (QTHICE(NV(I,1))+QTHICE(NV(I,2))+QTHICE(NV(I,3)))/3.0_SP
  END DO

     EL1 = ELF1
     D1  = H1 + EL1 - QTHICE1
     EL = ELF
     ET = EL
     D  = EL + H - QTHICE
     DT = D
     DTFA = D
     ET1 = EL1
     DT1 = D1


#   else
     EL1 = ELF1
     D1  = H1 + EL1
     EL = ELF
     ET = EL
     D  = EL + H
     DT = D
     DTFA = D
     ET1 = EL1
     DT1 = D1

#   endif

   END IF 

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SET_WD_DATA"
   RETURN
   END SUBROUTINE SET_WD_DATA

!==============================================================================|
!==============================================================================|

   SUBROUTINE ALLOC_WD_DATA  

!------------------------------------------------------------------------------|
!  ALLOCATE AND INITIALIZE WET/DRY TREATMENT ARRAYS                            |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE
   INTEGER NCT,NDB

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: ALLOC_WD_DATA"

# if !defined (DOUBLE_PRECISION)
   NDB = 1       !!GWC BASE THIS ON KIND
# else
   NDB = 2
# endif

!-----variables controlling porosities through wet/dry determination----------------!
                                                                                                                          
   ALLOCATE(ISWETN(0:MT))        ; ISWETN     = 1
   ALLOCATE(ISWETC(0:NT))        ; ISWETC     = 1
   ALLOCATE(ISWETNT(0:MT))       ; ISWETNT    = 1
   ALLOCATE(ISWETCT(0:NT))       ; ISWETCT    = 1
   ALLOCATE(ISWETCE(0:NT))       ; ISWETCE    = 1
!qxu{for inundation maps
   ALLOCATE(INUNDATION_MAPS(0:NT))      ; INUNDATION_MAPS   = 1
!qxu}

   memcnt = memcnt + 3*NT + 2*MT + NT

!!   ALLOCATE(US(0:NT,KB))         ;US    = ZERO   !!X-VELOCITY FOR MASS CONSERVATION
!!   ALLOCATE(VS(0:NT,KB))         ;VS    = ZERO   !!Y-VELOCITY FOR MASS CONSERVATION

!!   ALLOCATE(UAS(0:NT))           ;UAS       = ZERO  !!VERT AVGD X-VELOC FOR MASS CONSERVATION
!!   ALLOCATE(VAS(0:NT))           ;VAS       = ZERO  !!VERT AVGD Y-VELOC FOR MASS CONSERVATION
!!   ALLOCATE(UARDS(0:NT))         ;UARDS     = ZERO  !!UA AVGD OVER EXTERNAL INT FOR MASS CONSERVATION
!!   ALLOCATE(VARDS(0:NT))         ;VARDS     = ZERO  !!VA AVGD OVER EXTERNAL INT FOR MASS CONSERVATION

!!   memcnt = memcnt + KB*NT*NDB + KB*MT*NDB +NT*4*NDB

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: ALLOC_WD_DATA"
   RETURN
   END SUBROUTINE ALLOC_WD_DATA

!==============================================================================|
!==============================================================================|

   SUBROUTINE WET_JUDGE

!------------------------------------------------------------------------------|
!  DETERMINE IF NODES/ELEMENTS ARE WET OR DRY                                  |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE
   REAL(SP) :: DTMP
   INTEGER  :: ITA_TEMP
   INTEGER  :: I,IL,IA,IB,K1,K2,K3,K4,K5,K6

   integer :: KT

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: WET_JUDGE"
!
!--Determine If Node Points Are Wet/Dry Based on Depth Threshold---------------!
!
   ISWETN = 1
   DO I = 1, M
#   if defined (ICE_EMBEDDING)
     DTMP = H(I) + ELF(I)  - QTHICE(I)
#   else
     DTMP = H(I) + ELF(I)
#   endif

     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETN(I) = 0
   END DO

!
!--Determine if Cells are Wet/Dry Based on Depth Threshold---------------------!
!
   ISWETC = 1
   DO I = 1, N
#   if defined (ICE_EMBEDDING)
     DTMP =  MAX(ELF(NV(I,1)),ELF(NV(I,2)),ELF(NV(I,3)))  + &
             MIN(  H(NV(I,1)),  H(NV(I,2)),  H(NV(I,3)))  - &
             MAX(QTHICE(NV(I,1)),QTHICE(NV(I,2)),QTHICE(NV(I,3)))      !!yding
#   else
     DTMP =  MAX(ELF(NV(I,1)),ELF(NV(I,2)),ELF(NV(I,3)))  + &
             MIN(  H(NV(I,1)),  H(NV(I,2)),  H(NV(I,3)))
#   endif

     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETC(I) = 0
   END DO


   ! ------ Karsten Lettmann, 2016, May ----------------
   ! exchange the wet/dry information of elements between CPU
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     CALL AEXCHANGE(EC,MYID,NPROCS,ISWETC)
   END IF
#  endif
   ! --------------------------------------------------

!
!--A Secondary Condition for Nodal Dryness-(All Elements Around Node Are Dry)--!
!
   DO I = 1, M
     IF(SUM(ISWETC(NBVE(I,1:NTVE(I)))) == 0)  ISWETN(I) = 0
   END DO

!
!--Adjust Water Surface So It Does Not Go Below Minimum Depth------------------!
!
   ELF = MAX(ELF,-H + MIN_DEPTH)

!
!--Recompute Element Based Depths----------------------------------------------!
!
   DO I = 1, N
     ELF1(I) = ONE_THIRD*(ELF(NV(I,1))+ELF(NV(I,2))+ELF(NV(I,3)))
   END DO

!qxu{for inundation maps
!--Flood cells-----------------------------------------------------------------!
!
   INUNDATION_MAPS = 0
   WHERE (H1 < 0.0 )
     INUNDATION_MAPS = ISWETC
   END WHERE
!qxu}       
!
!--Extend Element/Node Based Wet/Dry Flags to Domain Halo----------------------!
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN

!JQI     CALL AEXCHANGE(EC,MYID,NPROCS,ISWETC)
     CALL AEXCHANGE(NC,MYID,NPROCS,ISWETN)
!qxu{for inundation maps
     CALL AEXCHANGE(EC,MYID,NPROCS,INUNDATION_MAPS)
!qxu}

   END IF
#  endif

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: WET_JUDGE"
   RETURN
   END SUBROUTINE WET_JUDGE

!==============================================================================|
!==============================================================================|

   SUBROUTINE WD_UPDATE(INCASE)

!------------------------------------------------------------------------------|
!  SHIFT WET/DRY VARIABLES TO NEW TIME LEVELS                                  |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: INCASE
   INTEGER :: I

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: WD_UPDATE"

   SELECT CASE(INCASE)

!------------------------------------------------------------------------------!
   CASE(1)    !! SHIFT AT END OF EXTERNAL MODE
!------------------------------------------------------------------------------!
   ISWETCE=ISWETC
!------------------------------------------------------------------------------!
   CASE(2)    !! UPDATE NODE WET/DRY AFTER DEPTH ADJUSTMENT
!------------------------------------------------------------------------------!
   DO I = 1,M
     IF(DTFA(I)-MIN_DEPTH <= 1.0E-5_SP) THEN
       ISWETN(I) = 0
     END IF
   END DO
   
! ------- New: Karsten Lettmann, 2017 Jan --------------------------   
#   if defined(MULTIPROCESSOR)
     ! ------- New Karsten Lettmann, Jan 2017 ------------
     IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ISWETN) !Interprocessor Exchange   !
#    endif
! ----------------- end new ----------------------------------------

!------------------------------------------------------------------------------!
   CASE(3)    !! SHIFT VARIABLES AT END OF INTERNAL MODE
!------------------------------------------------------------------------------!

   ISWETCT=ISWETC
   ISWETNT=ISWETN

   END SELECT

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: WD_UPDATE"
   RETURN
   END SUBROUTINE WD_UPDATE

END MODULE MOD_WD
