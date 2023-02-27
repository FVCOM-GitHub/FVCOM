!==============================================================================|

MODULE VARS_WAVE

  USE ALL_VARS
  USE MOD_OBCS

  IMPLICIT NONE
  SAVE

  INTEGER :: INP_CUR_NTIME,INP_WI_NTIME,INP_FR_NTIME,INP_WLEV_NTIME

  REAL, ALLOCATABLE :: AC2(:,:,:),COMPDA(:,:)
!  REAL, ALLOCATABLE :: ALPHA(:)

  LOGICAL :: NESTING

  CHARACTER(LEN=80) UGSWAN_VERSION !!STRING DESCRIBING VERSION
  
  INTEGER, ALLOCATABLE :: CROSS(:,:)                                  
  INTEGER, ALLOCATABLE :: BGRIDP(:)                                   
  REAL   , ALLOCATABLE :: BSPECS(:,:,:,:)                             
  REAL   , ALLOCATABLE :: AC1(:,:,:)                     
  REAL, ALLOCATABLE :: Sice(:,:,:)
  
  REAL   , ALLOCATABLE :: BLKND(:), BLKNDC(:), OURQT(:)   
  INTEGER :: ITW,IT0            

  REAL(SP)   , ALLOCATABLE :: HSC1(:), DIRDEG1(:), TPEAK(:), WLEN(:),QB1(:)
  REAl(SP)   , ALLOCATABLE :: CP(:)      ! Wave phase speed (m/s) ! Siqi Li, 2021-01-27
  REAL, ALLOCATABLE ::  OBC_HS(:,:),OBC_DIR(:,:),OBC_TPEAK(:,:)

  REAL   , ALLOCATABLE :: Pwave_bot(:),Ub_swan(:)
  REAL   , ALLOCATABLE :: Dwave(:) !Surface wind induced wave direction (radians)
  REAL(SP), ALLOCATABLE :: DIRBOT(:)
  REAL(SP), ALLOCATABLE :: SPEC_DENSITY(:,:)

  REAL(SP), ALLOCATABLE :: UWWIND(:)     !SURFACE X-WIND FOR WAVE MODEL
  REAL(SP), ALLOCATABLE :: VWWIND(:)     !SURFACE Y-WIND FOR WAVE MODEL

  TYPE(TIME) :: WaveTime,IMDTW
  
  REAL :: SOURCE_DTMAX,SOURCE_DTMIN
  CHARACTER(LEN=4) COMPUT

  REAL(SP), ALLOCATABLE :: HS_WIND(:),DIRDEG_WIND(:),TPEAK_WIND(:)
  INTEGER,  ALLOCATABLE :: TPEAK_WIND_POS(:)
  REAL(SP), ALLOCATABLE :: HS_SWELL_ALL(:,:),DIRDEG_SWELL_ALL(:,:),TPEAK_SWELL_ALL(:,:)
  INTEGER,  ALLOCATABLE :: TPEAK_SWELL_POS_ALL(:,:)

END MODULE VARS_WAVE
