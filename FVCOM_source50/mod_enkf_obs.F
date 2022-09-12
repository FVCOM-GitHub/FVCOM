!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

MODULE MOD_ENKF_OBS
# if defined (ENKF)
  USE ENKFVAL
  USE CONTROL
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_PREC
  USE MOD_WD
  USE ALL_VARS, only: vxmin,vymin
  USE MOD_STARTUP
  USE MOD_INPUT

  IMPLICIT NONE
  SAVE

  REAL(SP),ALLOCATABLE :: AW0G(:,:),AWXG(:,:),AWYG(:,:)

!--Data Assimilation Parameters for Current Nudging  Assimilation
!

  INTEGER  :: CUR_MAX_LAYER          !!MAXIMUM NUMBER OF VERTICAL DATA FROM ANY OBS POINT


!
  !
!
!
!--Current Assimilation Object Type
!
  TYPE ASSIM_OBJ_CUR
   INTEGER  :: N_TIMES                     !!NUMBER OF DATA TIMES
   INTEGER  :: N_INTPTS                    !!NUMBER OF INTERPOLATION POINTS 
   INTEGER  :: N_T_WEIGHT                  !!DATA TIME FOR CURRENT OBSERVATION WEIGHTING
   INTEGER  :: N_LAYERS                    !!NUMBER OF OBSERVATIONS IN THE VERTICAL
   INTEGER  :: N_CELL                      !!CELL NUMBER OF OBSERVATION
   REAL(SP) :: X,Y                         !!X AND Y COORDINATES OF OBSERVATION
   REAL(SP) :: T_WEIGHT                    !!TIME WEIGHT
   REAL(SP) :: DEPTH                       !!DEPTH AT OBSERVATION STATION (MOORING)
   REAL(SP) :: SITA                        !!LOCAL ISOBATH ANGLE AT OBSERVATION
   REAL(SP),ALLOCATABLE,DIMENSION(:)  :: ODEPTH      !!OBSERVATION DEPTHS
   INTEGER                        :: chose_flag      !!0 IS NO, 1 IS YES
   INTEGER                        :: INNODE     !! OBS NEAREST NODE
   REAL(DP),ALLOCATABLE,DIMENSION(:)  :: obs_days       !!DATA TIMES
   TYPE(TIME),ALLOCATABLE,DIMENSION(:)  :: TIMES       !!DATA TIMES
   REAL(SP),ALLOCATABLE,DIMENSION(:,:):: UO,VO       !!OBSERVATION DATA FOR X,Y VELOCITY COMPONENTS
   INTEGER,ALLOCATABLE,DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC 
   INTEGER,ALLOCATABLE,DIMENSION(:,:) :: S_INT       !!SIGMA INTERVALS SURROUNDING CURRENT MEASUREMENT
   REAL(SP),ALLOCATABLE,DIMENSION(:,:):: S_WEIGHT    !!SIGMA WEIGHTING                               
   REAL(SP),ALLOCATABLE,DIMENSION(:)  :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
  END TYPE ASSIM_OBJ_CUR

!
!--Temperature/Salinity Assimilation Object Type
!
   TYPE ASSIM_OBJ_TS  
     INTEGER  :: N_TIMES                                  !!NUMBER OF DATA TIMES
     INTEGER  :: N_INTPTS                                 !!NUMBER OF INTERPOLATION POINTS
     INTEGER  :: N_T_WEIGHT                               !!DATA TIME FOR CURRENT OBSERVATION WEIGHTING
     INTEGER  :: N_LAYERS                                 !!NUMBER OF OBSERVATIONS IN THE VERTICAL
     INTEGER  :: N_CELL                                   !!CELL NUMBER OF OBSERVATION
     REAL(SP) :: X,Y                                      !!X AND Y COORDINATES OF OBSERVATION
     REAL(SP) :: T_WEIGHT                                 !!TIME WEIGHT
     REAL(SP) :: DEPTH                                    !!DEPTH AT OBSERVATION STATION (MOORING)
     REAL(SP) :: SITA                                     !!LOCAL ISOBATH ANGLE AT OBSERVATION
   REAL(SP),ALLOCATABLE,DIMENSION(:)  :: ODEPTH      !!OBSERVATION DEPTHS
   INTEGER                        :: chose_flag      !!0 IS NO, 1 IS YES
   INTEGER                        :: INNODE      !! OBS NEAREST NODE
   REAL(DP),ALLOCATABLE,DIMENSION(:)  :: obs_days       !!DATA TIMES
     TYPE(TIME),ALLOCATABLE,DIMENSION(:)  :: TIMES        !!DATA TIMES
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP        !!OBSERVATION DATA FOR TEMPERATURE 
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: SAL         !!OBSERVATION DATA FOR SALINITY 
     INTEGER,  ALLOCATABLE, DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC
     INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: S_INT       !!SIGMA INTERVALS SURROUNDING CURRENT MEASUREMENT
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: S_WEIGHT    !!SIGMA WEIGHTING
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
   END TYPE ASSIM_OBJ_TS    
!--EL ASSIMILATION Assimilation Object Type
!
   TYPE ASSIM_OBJ_EL
     INTEGER  :: N_TIMES                                  !!NUMBER OF DATA TIMES
     INTEGER  :: N_INTPTS                                 !!NUMBER OF INTERPOLATION POINTS
     INTEGER  :: N_T_WEIGHT                               !!DATA TIME FOR CURRENT OBSERVATION WEIGHTING
     INTEGER  :: N_LAYERS                                 !!NUMBER OF OBSERVATIONS IN THE VERTICAL
     INTEGER  :: N_CELL                                   !!CELL NUMBER OF OBSERVATION
     REAL(SP) :: X,Y                                      !!X AND Y COORDINATES OF OBSERVATION
     REAL(SP) :: T_WEIGHT                                 !!TIME WEIGHT
     REAL(SP) :: DEPTH                                    !!DEPTH AT OBSERVATION STATION (MOORING)
     REAL(SP) :: SITA                                     !!LOCAL ISOBATH ANGLE AT OBSERVATION
   REAL(SP),ALLOCATABLE,DIMENSION(:)  :: ODEPTH      !!OBSERVATION DEPTHS
   INTEGER                        :: chose_flag      !!0 IS NO, 1 IS YES
   INTEGER                        :: INNODE      !! OBS NEAREST NODE
   REAL(DP),ALLOCATABLE,DIMENSION(:)  :: obs_days       !!DATA TIMES
     TYPE(TIME),ALLOCATABLE,DIMENSION(:)  :: TIMES        !!DATA TIMES
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: EL        !!OBSERVATION DATA FOR TEMPERATURE
     INTEGER,  ALLOCATABLE, DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC
     INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: S_INT       !!SIGMA INTERVALS SURROUNDING CURRENT MEASUREMENT
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: S_WEIGHT    !!SIGMA WEIGHTING
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
   END TYPE ASSIM_OBJ_EL

!  el assimilation variables
   INTEGER                          :: N_ASSIM_EL !!NUMBER OF CURRENT OBSERVATIONS
   TYPE(ASSIM_OBJ_EL), ALLOCATABLE :: EL0_OBS(:)  !!CU ASSIMILATION DATA OBJECTS

                                                                                                                       
!
!--Current Data Assimilation Variables
!
   INTEGER                          :: N_ASSIM_CUR !!NUMBER OF CURRENT OBSERVATIONS 
   TYPE(ASSIM_OBJ_CUR), ALLOCATABLE :: CUR_OBS(:)  !!CURRENT ASSIMILATION DATA OBJECTS

!
!--Salinity/Temperature Data Assimilation Variables
!
   INTEGER                                N_ASSIM_T   !!NUMBER OF TEMPERATURE OBSERVATIONS
   TYPE(ASSIM_OBJ_TS), ALLOCATABLE     :: T0_OBS(:)    !!TEMP ASSIMILATION DATA OBJECTS


   INTEGER                                N_ASSIM_S   !!NUMBER OF TEMPERATURE OBSERVATIONS
   TYPE(ASSIM_OBJ_TS), ALLOCATABLE     :: S0_OBS(:)    !!TEMP ASSIMILATION DATA OBJECTS
CONTAINS                                                                                                                       
   subroutine deallocate_obs_data
   implicit none
   IF (EL_ASSIM) THEN
      DEALLOCATE(EL0_OBS)
   END IF
   IF (UV_ASSIM) THEN
      DEALLOCATE(CUR_OBS)
   END IF
   IF (T_ASSIM) THEN
       DEALLOCATE(T0_OBS)
   END IF
   IF  (S_ASSIM)  THEN
      DEALLOCATE(S0_OBS)
   END IF

   end subroutine deallocate_obs_data
   subroutine filter_obs_data
   IMPLICIT NONE
   INTEGER :: I,J,K
   IF (EL_ASSIM) THEN
      CALL FILTER_EL_OBS_DATA(INTTIME)
   END IF
   IF (UV_ASSIM) THEN
      CALL FILTER_UV_OBS_DATA(INTTIME)
   END IF
   IF (T_ASSIM) THEN
      CALL FILTER_T_OBS_DATA(INTTIME)
   END IF
   IF  (S_ASSIM)  THEN
      CALL FILTER_S_OBS_DATA(INTTIME)
   END IF 
   end subroutine filter_obs_data
   subroutine filter_el_obs_data(current_time)
   IMPLICIT NONE
   real(dp) :: time1
   TYPE(TIME) :: current_Time
!!-------------------------------------------------
   INTEGER I,J,K,NN,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   LOGICAL :: FEXIST

   INTEGER ::IERR,Nsite_tmp
   REAL(DP) :: ODAYS
!   current_time=days2time(time1)
      FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_el.xy"
!--Make Sure Current Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'EL OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
!
!!
!!--Read Number of Current Measurement Stations---------------------------------!
!!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_EL

   ALLOCATE(EL0_OBS(N_ASSIM_EL))

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!

   DO I=1,N_ASSIM_EL
     READ(1,*)ITMP,EL0_OBS(I)%X,EL0_OBS(I)%Y !,EL0_OBS(I)%DEPTH,NLAY,CUR_OBS(I)%SITA
     NLAY=1
     EL0_OBS(I)%N_LAYERS = NLAY

     ALLOCATE(EL0_OBS(I)%ODEPTH(NLAY))
     EL0_OBS(I)%chose_flag=0
     DO J=1,NLAY
       READ(1,*)EL0_OBS(I)%ODEPTH(J)
     END DO
   END DO
   CLOSE(1)
!------------------------------------------------------------------------------!
!  Open Current Observation Files for Each Observation Point and Read Data     !
!------------------------------------------------------------------------------!
   ONAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_el.dat"
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'EL OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_EL
!----Count Number of Data Entries for Observation I----------------------------!

     READ(1,*,IOSTAT=IOS) Nsite_tmp,NCNT     
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_el.dat at site number:',I
       CALL PSTOP
     END IF
     EL0_OBS(I)%N_TIMES = NCNT
     
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR EL OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Current (UA,VA) and Time (TIME)-------------------!
     ALLOCATE(EL0_OBS(I)%TIMES(EL0_OBS(I)%N_TIMES))
     ALLOCATE(EL0_OBS(I)%EL(EL0_OBS(I)%N_TIMES , EL0_OBS(I)%N_LAYERS ))
     ALLOCATE(EL0_OBS(I)%obs_days(EL0_OBS(I)%N_TIMES))
!----Read in Current Data for Observation I------------------------------------!
     NLAY = EL0_OBS(I)%N_LAYERS
     DO J=1,EL0_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS) ODAYS,(EL0_OBS(I)%EL(J,K),K=1,NLAY)
       IF(IOS < 0)THEN         
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_el.dat at site:',I,'Time No:',J
         CALL PSTOP
       END IF  
       EL0_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
       EL0_OBS(I)%obs_days(J) = ODAYS
       if (current_Time == EL0_OBS(I)%TIMES(J)) then !  .and. ABS(CUR_OBS(I)%UO(J,K)) + ABS(CUR_OBS(I)%VO(J,K)) > .0001) then
           EL0_OBS(I)%chose_flag=1
       end if
     END DO

  END DO
  CLOSE(1)
!---------------------------------------write out obs data for current time
  itmp=0
  IF (MSR) THEN
  ONAME = "enkf_el.xy"
   OPEN(1,FILE=ONAME,STATUS='replace')  ; REWIND(1)
    write(1,*) sum(EL0_OBS(:)%chose_flag)
   DO I=1,N_ASSIM_EL
     if (EL0_OBS(I)%chose_flag==1) then
     itmp=itmp+1
     write(1,'(i8,3f18.6,i8,f12.6)')ITMP,EL0_OBS(I)%X,EL0_OBS(I)%Y,0.0,EL0_OBS(I)%N_LAYERS
     DO J=1,EL0_OBS(I)%N_LAYERS
       write(1,*)EL0_OBS(I)%ODEPTH(J)
     END DO
     end if
   END DO
   CLOSE(1)
     itmp=0
   ONAME = "enkf_el.dat"
   OPEN(1,FILE=ONAME,STATUS='replace')  ; REWIND(1)
   DO I=1,N_ASSIM_EL 
     if (EL0_OBS(I)%chose_flag==1) then
         itmp=itmp+1
         write(1,*) itmp,1
         DO J=1,EL0_OBS(I)%N_TIMES
          if (EL0_OBS(I)%TIMES(J)==current_time) then
              WRITE(1,'(f18.8,200f12.6)') EL0_OBS(I)%obs_days(j),(EL0_OBS(I)%EL(J,K),K=1,EL0_OBS(I)%N_LAYERS)
          end if
         END DO
     end if
   END DO
  CLOSE(1)
  END IF
   DEALLOCATE(EL0_OBS)
   end subroutine filter_el_obs_data
!-------------------------------------------------------------------------------



   subroutine filter_uv_obs_data(current_time)
   IMPLICIT NONE
   real(dp) :: time1
   TYPE(TIME) :: current_Time
!!-------------------------------------------------
   INTEGER I,J,K,NN,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   LOGICAL :: FEXIST

   INTEGER ::IERR,Nsite_tmp
   REAL(DP) :: ODAYS
!   current_time=days2time(time1)
      FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_cur.xy"
!--Make Sure Current Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
!
!!
!!--Read Number of Current Measurement Stations---------------------------------!
!!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_CUR

   ALLOCATE(CUR_OBS(N_ASSIM_CUR))

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!

   DO I=1,N_ASSIM_CUR
     READ(1,*)ITMP,CUR_OBS(I)%X,CUR_OBS(I)%Y,CUR_OBS(I)%DEPTH,NLAY,CUR_OBS(I)%SITA
     CUR_OBS(I)%N_LAYERS = NLAY

     ALLOCATE(CUR_OBS(I)%ODEPTH(NLAY))
     CUR_OBS(I)%chose_flag=0
     DO J=1,NLAY
       READ(1,*)CUR_OBS(I)%ODEPTH(J)
       IF(CUR_OBS(I)%ODEPTH(J) > CUR_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'CURRENT OBSERVATION LAYER',J,'OF CURRENT MOORING',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH=',CUR_OBS(I)%ODEPTH(J),CUR_OBS(I)%DEPTH
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO
   CLOSE(1)
!------------------------------------------------------------------------------!
!  Open Current Observation Files for Each Observation Point and Read Data     !
!------------------------------------------------------------------------------!
   ONAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_cur.dat"
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_CUR
!----Count Number of Data Entries for Observation I----------------------------!

     READ(1,*,IOSTAT=IOS) Nsite_tmp,NCNT     
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_cur.dat at site number:',I
       CALL PSTOP
     END IF
     CUR_OBS(I)%N_TIMES = NCNT
     
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR CURRENT OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Current (UA,VA) and Time (TIME)-------------------!
     ALLOCATE(CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES))
     ALLOCATE(CUR_OBS(I)%UO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))
     ALLOCATE(CUR_OBS(I)%VO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))
     ALLOCATE(CUR_OBS(I)%obs_days(CUR_OBS(I)%N_TIMES))
!----Read in Current Data for Observation I------------------------------------!
     NLAY = CUR_OBS(I)%N_LAYERS
     DO J=1,CUR_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS) ODAYS,(CUR_OBS(I)%UO(J,K),CUR_OBS(I)%VO(J,K),K=1,NLAY)
       DO K=1,NLAY
          if ( ABS(CUR_OBS(I)%UO(J,K)) + ABS(CUR_OBS(I)%VO(J,K)) < .0001) then
              PRINT *, 'CUR_OBS' ,i, 'VELOCITY', J,K ,'IS ZERO WRONG'
              CALL PSTOP
          end if
       END DO
       IF(IOS < 0)THEN         
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_cur.dat at site:',I,'Time No:',J
         CALL PSTOP
       END IF  
       CUR_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
       CUR_OBS(I)%obs_days(J) = ODAYS
       if (current_Time == CUR_OBS(I)%TIMES(J)) then !  .and. ABS(CUR_OBS(I)%UO(J,K)) + ABS(CUR_OBS(I)%VO(J,K)) > .0001) then
           CUR_OBS(I)%chose_flag=1
       end if
     END DO

  END DO
  CLOSE(1)
!---------------------------------------write out obs data for current time
  itmp=0
  ONAME = "enkf_cur.xy"
  IF (MSR) THEN
   OPEN(1,FILE=ONAME,STATUS='replace')  ; REWIND(1)
    write(1,*) sum(CUR_OBS(:)%chose_flag)
   DO I=1,N_ASSIM_CUR
     if (CUR_OBS(I)%chose_flag==1) then
     itmp=itmp+1
     write(1,'(i8,3f18.6,i8,f12.6)')ITMP,CUR_OBS(I)%X,CUR_OBS(I)%Y,CUR_OBS(I)%DEPTH,CUR_OBS(I)%N_LAYERS,CUR_OBS(I)%SITA
     DO J=1,CUR_OBS(I)%N_LAYERS
       write(1,*)CUR_OBS(I)%ODEPTH(J)
     END DO
     end if
   END DO
   CLOSE(1)
     itmp=0
   ONAME = "enkf_cur.dat"
   OPEN(1,FILE=ONAME,STATUS='replace')  ; REWIND(1)
   DO I=1,N_ASSIM_CUR 
     if (CUR_OBS(I)%chose_flag==1) then
         itmp=itmp+1
         write(1,*) itmp,1
         DO J=1,CUR_OBS(I)%N_TIMES
          if (CUR_OBS(I)%TIMES(J)==current_time) then
              WRITE(1,'(f18.8,200f12.6)') CUR_OBS(I)%obs_days(j),(CUR_OBS(I)%UO(J,K),CUR_OBS(I)%VO(J,K),K=1,CUR_OBS(I)%N_LAYERS)
          end if
         END DO
     end if
   END DO
   CLOSE(1)
   END IF
   DEALLOCATE(CUR_OBS)
   end subroutine filter_uv_obs_data
!-------------------------------------------------------------------------------

   subroutine filter_t_obs_data(current_time)
   IMPLICIT NONE
   real(dp) :: time1
   TYPE(TIME) :: current_Time
!!-------------------------------------------------
   INTEGER I,J,K,NN,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   LOGICAL :: FEXIST

   INTEGER ::IERR,Nsite_tmp
   REAL(DP) :: ODAYS
!   current_time=days2time(time1)
      FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_t.xy"
!--Make Sure Current Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'TEMPERATURE OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
!
!!
!!--Read Number of Current Measurement Stations---------------------------------!
!!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_T

   ALLOCATE(T0_OBS(N_ASSIM_T))

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!

   DO I=1,N_ASSIM_T
     READ(1,*)ITMP,T0_OBS(I)%X,T0_OBS(I)%Y,T0_OBS(I)%DEPTH,NLAY,T0_OBS(I)%SITA
     T0_OBS(I)%N_LAYERS = NLAY

     ALLOCATE(T0_OBS(I)%ODEPTH(NLAY))
     T0_OBS(I)%chose_flag=0
     DO J=1,NLAY
       READ(1,*)T0_OBS(I)%ODEPTH(J)
       IF(T0_OBS(I)%ODEPTH(J) > T0_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'TEMPERATURE OBSERVATION LAYER',J,'OF TEMPERATURE MOORING',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH=',T0_OBS(I)%ODEPTH(J),T0_OBS(I)%DEPTH
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO
   CLOSE(1)
!------------------------------------------------------------------------------!
!  Open Current Observation Files for Each Observation Point and Read Data     !
!------------------------------------------------------------------------------!
   ONAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_t.dat"
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'TEMPERATURE OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_T
!----Count Number of Data Entries for Observation I----------------------------!

     READ(1,*,IOSTAT=IOS) Nsite_tmp,NCNT     
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_t.dat at site number:',I
       CALL PSTOP
     END IF
     T0_OBS(I)%N_TIMES = NCNT
     
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR CURRENT OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Current (UA,VA) and Time (TIME)-------------------!
     ALLOCATE(T0_OBS(I)%TIMES(T0_OBS(I)%N_TIMES))
     ALLOCATE(T0_OBS(I)%TEMP( T0_OBS(I)%N_TIMES , T0_OBS(I)%N_LAYERS ))
     ALLOCATE(T0_OBS(I)%obs_days(T0_OBS(I)%N_TIMES))
!----Read in Current Data for Observation I------------------------------------!
     NLAY = T0_OBS(I)%N_LAYERS
     DO J=1,T0_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS) ODAYS,(T0_OBS(I)%TEMP(J,K),K=1,NLAY)
       DO K=1,NLAY
          if ( T0_OBS(I)%TEMP(J,K)< -90.) then
              PRINT *, 'T0_OBS' ,i, 'VELOCITY', J,K ,'IS ZERO WRONG'
              CALL PSTOP
          end if
       END DO
       IF(IOS < 0)THEN         
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_t.dat at site:',I,'Time No:',J
         CALL PSTOP
       END IF  
       T0_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
       T0_OBS(I)%obs_days(J) = ODAYS
       if (current_Time == T0_OBS(I)%TIMES(J)) then !  .and. ABS(CUR_OBS(I)%UO(J,K)) + ABS(CUR_OBS(I)%VO(J,K)) > .0001) then
           T0_OBS(I)%chose_flag=1
       end if
     END DO

  END DO
  CLOSE(1)
!---------------------------------------write out obs data for current time
  itmp=0
  ONAME = "enkf_t.xy"
  IF (MSR) THEN
   OPEN(1,FILE=ONAME,STATUS='replace')  ; REWIND(1)
    write(1,*) sum(T0_OBS(:)%chose_flag)
   DO I=1,N_ASSIM_T
     if (T0_OBS(I)%chose_flag==1) then
     itmp=itmp+1
     write(1,'(i8,3f18.6,i8,f12.6)')ITMP,T0_OBS(I)%X,T0_OBS(I)%Y,T0_OBS(I)%DEPTH,T0_OBS(I)%N_LAYERS,T0_OBS(I)%SITA
     DO J=1,T0_OBS(I)%N_LAYERS
       write(1,*)T0_OBS(I)%ODEPTH(J)
     END DO
     end if
   END DO
   CLOSE(1)
     itmp=0
   ONAME = "enkf_t.dat"
   OPEN(1,FILE=ONAME,STATUS='replace')  ; REWIND(1)
   DO I=1,N_ASSIM_T 
     if (T0_OBS(I)%chose_flag==1) then
         itmp=itmp+1
         write(1,*) itmp,1
         DO J=1,T0_OBS(I)%N_TIMES
          if (T0_OBS(I)%TIMES(J)==current_time) then
              WRITE(1,'(f18.8,200f12.6)') T0_OBS(I)%obs_days(j),(T0_OBS(I)%TEMP(J,K),K=1,T0_OBS(I)%N_LAYERS)
          end if
         END DO
     end if
   END DO
   CLOSE(1)
   END IF
   DEALLOCATE(T0_OBS)
   end subroutine filter_t_obs_data
!-------------------------------------------------------------------------------

   subroutine filter_s_obs_data(current_time)
   IMPLICIT NONE
   real(dp) :: time1
   TYPE(TIME) :: current_Time
!!-------------------------------------------------
   INTEGER I,J,K,NN,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   LOGICAL :: FEXIST

   INTEGER ::IERR,Nsite_tmp
   REAL(DP) :: ODAYS
!   current_time=days2time(time1)
      FNAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_s.xy"
!--Make Sure Current Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'SALINITY OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
!
!!
!!--Read Number of Current Measurement Stations---------------------------------!
!!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_S

   ALLOCATE(S0_OBS(N_ASSIM_S))

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!

   DO I=1,N_ASSIM_S
     READ(1,*)ITMP,S0_OBS(I)%X,S0_OBS(I)%Y,S0_OBS(I)%DEPTH,NLAY,S0_OBS(I)%SITA
     S0_OBS(I)%N_LAYERS = NLAY

     ALLOCATE(S0_OBS(I)%ODEPTH(NLAY))
     S0_OBS(I)%chose_flag=0
     DO J=1,NLAY
       READ(1,*)S0_OBS(I)%ODEPTH(J)
       IF(S0_OBS(I)%ODEPTH(J) > S0_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'SALINITY OBSERVATION LAYER',J,'OF SALINITY MOORING',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH=',S0_OBS(I)%ODEPTH(J),S0_OBS(I)%DEPTH
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO
   CLOSE(1)
!------------------------------------------------------------------------------!
!  Open Current Observation Files for Each Observation Point and Read Data     !
!------------------------------------------------------------------------------!
   ONAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_s.dat"
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'SALINITY OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_S
!----Count Number of Data Entries for Observation I----------------------------!

     READ(1,*,IOSTAT=IOS) Nsite_tmp,NCNT     
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_s.dat at site number:',I
       CALL PSTOP
     END IF
     S0_OBS(I)%N_TIMES = NCNT
     
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR SALINITY OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Current (UA,VA) and Time (TIME)-------------------!
     ALLOCATE(S0_OBS(I)%TIMES(S0_OBS(I)%N_TIMES))
     ALLOCATE(S0_OBS(I)%SAL( S0_OBS(I)%N_TIMES , S0_OBS(I)%N_LAYERS ))
     ALLOCATE(S0_OBS(I)%obs_days(S0_OBS(I)%N_TIMES))
!----Read in Current Data for Observation I------------------------------------!
     NLAY = S0_OBS(I)%N_LAYERS
     DO J=1,S0_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS) ODAYS,(S0_OBS(I)%SAL(J,K),K=1,NLAY)
       DO K=1,NLAY
          if ( S0_OBS(I)%SAL(J,K) <-90) then
              PRINT *, 'SAL_OBS' ,i, 'VELOCITY', J,K ,'IS ZERO WRONG'
              CALL PSTOP
          end if
       END DO
       IF(IOS < 0)THEN         
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_s.dat at site:',I,'Time No:',J
         CALL PSTOP
       END IF  
       S0_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
       S0_OBS(I)%obs_days(J) = ODAYS
       if (current_Time == S0_OBS(I)%TIMES(J)) then !  .and. ABS(CUR_OBS(I)%UO(J,K)) + ABS(CUR_OBS(I)%VO(J,K)) > .0001) then
           S0_OBS(I)%chose_flag=1
       end if
     END DO

  END DO
  CLOSE(1)
!---------------------------------------write out obs data for current time
  itmp=0
  ONAME = "enkf_s.xy"
  IF (MSR) THEN
   OPEN(1,FILE=ONAME,STATUS='replace')  ; REWIND(1)
    write(1,*) sum(S0_OBS(:)%chose_flag)
   DO I=1,N_ASSIM_S
     if (S0_OBS(I)%chose_flag==1) then
     itmp=itmp+1
     write(1,'(i8,3f18.6,i8,f12.6)')ITMP,S0_OBS(I)%X,S0_OBS(I)%Y,S0_OBS(I)%DEPTH,S0_OBS(I)%N_LAYERS,S0_OBS(I)%SITA
     DO J=1,S0_OBS(I)%N_LAYERS
       write(1,*)S0_OBS(I)%ODEPTH(J)
     END DO
     end if
   END DO
   CLOSE(1)
     itmp=0
   ONAME = "enkf_s.dat"
   OPEN(1,FILE=ONAME,STATUS='replace')  ; REWIND(1)
   DO I=1,N_ASSIM_S 
     if (S0_OBS(I)%chose_flag==1) then
         itmp=itmp+1
         write(1,*) itmp,1
         DO J=1,S0_OBS(I)%N_TIMES
          if (S0_OBS(I)%TIMES(J)==current_time) then
              WRITE(1,'(f18.8,200f12.6)') S0_OBS(I)%obs_days(j),(S0_OBS(I)%SAL(J,K),K=1,S0_OBS(I)%N_LAYERS)
          end if
         END DO
     end if
   END DO
   CLOSE(1)
   END IF
   DEALLOCATE(S0_OBS)
   end subroutine filter_s_obs_data
!-------------------------------------------------------------------------------
   subroutine set_enkf_assim_data
   IMPLICIT NONE
   INTEGER :: I,J,K
   NLOC=0
   IF (EL_ASSIM) THEN
       CALL SET_ENKF_EL_ASSIM_DATA
       NLOC=NLOC+N_ASSIM_EL
   END IF
   IF (UV_ASSIM) THEN
       CALL SET_ENKF_CUR_ASSIM_DATA
      DO I=1,N_ASSIM_CUR
          NLOC=NLOC+2*CUR_OBS(I)%N_LAYERS
      END DO
   END IF
   IF (T_ASSIM) THEN
       CALL SET_ENKF_T_ASSIM_DATA
      DO I=1,N_ASSIM_T
          NLOC=NLOC+T0_OBS(I)%N_LAYERS
      END DO
   END IF
   IF (S_ASSIM) THEN
       print *, 'start to call set_enkf_s_assim_data'
       CALL SET_ENKF_S_ASSIM_DATA
      DO I=1,N_ASSIM_S
          NLOC=NLOC+S0_OBS(I)%N_LAYERS
      END DO
   END IF
   end subroutine set_enkf_assim_data
!======================================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
   SUBROUTINE SET_ENKF_EL_ASSIM_DATA 
!
!!------------------------------------------------------------------------------!
!!  SET UP ASSIMILATION DATA FOR TEMP/SAL OBSERVATIONS                          |
!!------------------------------------------------------------------------------!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)   
   USE MOD_PAR  
# endif
   IMPLICIT NONE
   INTEGER I,J,K,NN,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   CHARACTER(LEN= 2 ) :: NAC   
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
   REAL(SP), PARAMETER :: ALF = 0.05_SP
   LOGICAL :: FEXIST
   INTEGER :: MAXEL,NBD_CNT
   REAL(SP) :: TS_RADIUS 
   INTEGER :: JMIN,INTMP(1) 
   REAL(SP):: LMIN

   REAL(SP), DIMENSION(1:NGL,1) :: RDLIST
   REAL(SP), DIMENSION(3) :: XTRI,YTRI
   REAL(SP) :: RDLAST
   INTEGER :: LOCIJ(2),MIN_LOC,JJ,IERR,Nsite_tmp
   INTEGER :: ND1,ND2,ND3
   REAL(SP) :: DELTA,COFA,COFB,COFC
   REAL(SP) ::S11,S22,S33,RTMP,RRTMP
   REAL(SP), DIMENSION(KB) :: ZZ_OB
!#  if defined (MULTIPROCESSOR)   
   REAL(SP), ALLOCATABLE :: ZZ_G(:,:)  
!# endif
   REAL(SP) :: X0C,Y0C,HX,HY,H0,TEMP_H
   REAL(DP) :: ODAYS
      
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_TS_ASSIM_DATA" 

!------------------------------------------------------------------------------!
!  Read Number of Scalar Observations and Coordinates of Each                  !
!------------------------------------------------------------------------------!
       
   FNAME = "enkf_el.xy"
!
!--Make Sure elevation Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'EL OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
     
!
!--Read Number of T/S Measurement Stations---------------------------------!
!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD') ;REWIND(1)
   READ(1,*) N_ASSIM_EL                  !nomber of TS Obs station
   ALLOCATE(EL0_OBS(N_ASSIM_EL))          !Type for TS_OBS

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!
   DO I=1,N_ASSIM_EL  
     READ(1,*)ITMP,EL0_OBS(I)%X,EL0_OBS(I)%Y,EL0_OBS(I)%DEPTH,NLAY,EL0_OBS(I)%SITA
     EL0_OBS(I)%N_LAYERS = NLAY
     ALLOCATE(EL0_OBS(I)%ODEPTH(NLAY))
     DO J=1,NLAY
       READ(1,*)EL0_OBS(I)%ODEPTH(J)
       IF(EL0_OBS(I)%ODEPTH(J) > EL0_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'OBSERVATION DEPTH',J,'OF EL OBS',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH=',EL0_OBS(I)%ODEPTH(J),EL0_OBS(I)%DEPTH
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO

!
!--Shift Coordinates-----------------------------------------------------------!
!
# if !defined (SPHERICAL) 
   EL0_OBS(:)%X = EL0_OBS(:)%X - VXMIN 
   EL0_OBS(:)%Y = EL0_OBS(:)%Y - VYMIN 
# endif
   
!# if defined (SPHERICAL) 
   IF(SERIAL)THEN
     XCG = XC
     YCG = YC
     XG = VX
     YG = VY
   END IF
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,XC,XCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,YC,YCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VX,XG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VY,YG)
   IF(PAR)CALL MPI_BCAST(XCG,NGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(YCG,NGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(XG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(YG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#  endif
!# endif
!  
!
!--find the cell number (TS_OBS(:)%N_CELL) of Obs station---------------------
!
   RRTMP = 100000.0     !100km
   DO J= 1,N_ASSIM_EL 
     X0 = EL0_OBS(J)%X
     Y0 = EL0_OBS(J)%Y
     DO I=1,NGL
#    if !defined (SPHERICAL)
       Rtmp = SQRT((XCG(I)-X0)*(XCG(I)-X0)+(YCG(I)-Y0)*(YCG(I)-Y0))
#    else
       CALL ARC(X0,Y0,XCG(I),YCG(I),Rtmp)
#    endif
       IF(Rtmp < RRTMP) THEN
	 S11 = (XG(NVG(I,2))-X0)*(YG(NVG(I,3))-Y0)-&
	       (XG(NVG(I,3))-X0)*(YG(NVG(I,2))-Y0)
         S22 = (XG(NVG(I,3))-X0)*(YG(NVG(I,1))-Y0)-&
	       (XG(NVG(I,1))-X0)*(YG(NVG(I,3))-Y0)
         S33 = (XG(NVG(I,1))-X0)*(YG(NVG(I,2))-Y0)-&
	       (XG(NVG(I,2))-X0)*(YG(NVG(I,1))-Y0)
	 IF(S11 <= 0._SP .AND. S22 <= 0._SP .AND. S33 <= 0._SP) THEN
	   EL0_OBS(J)%N_CELL = I
           S11=SQRT((XG(NVG(I,1))-X0)**2+(YG(NVG(I,1))-Y0)**2) ! CACLULATE DISTANCE
           S22=SQRT((XG(NVG(I,2))-X0)**2+(YG(NVG(I,2))-Y0)**2) 
           S33=SQRT((XG(NVG(I,3))-X0)**2+(YG(NVG(I,3))-Y0)**2)
           intmp=minloc([s11,s22,s33])
           EL0_OBS(J)%innode=NVG(I,intmp(1))
	   GOTO 300
	 ELSE
	   EL0_OBS(J)%N_CELL = 0
	 ENDIF
       ELSE
	 EL0_OBS(J)%N_CELL = -1
       ENDIF
     ENDDO	
300  CONTINUE
     IF(EL0_OBS(J)%N_CELL <= 0) THEN
       IF(MSR) WRITE(IPT,*)'ERROR--EL OBS SITE:',J,' OUT OF DOMAN',&
	                    EL0_OBS(J)%N_CELL 
       CALL PSTOP
     ENDIF 	 

   ENDDO             
!--Gather AW0G,AWXG & AWYG use for interp grid to Obs station
   IF(.NOT. ALLOCATED(AW0G)) ALLOCATE(AW0G(0:NGL,3))
   IF(.NOT. ALLOCATED(AWXG)) ALLOCATE(AWXG(0:NGL,3))
   IF(.NOT. ALLOCATED(AWYG)) ALLOCATE(AWYG(0:NGL,3))
   
   IF(SERIAL)THEN
     AW0G = AW0
     AWXG = AWX
     AWYG = AWY
   END IF
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN     
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AW0,AW0G)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWX,AWXG)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWY,AWYG)
     CALL MPI_BCAST(AW0G,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
     CALL MPI_BCAST(AWXG,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
     CALL MPI_BCAST(AWYG,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   END IF
#  endif     
!
!--Close Current Observation Global File---------------------------------------!
!
   CLOSE(1)
       
!------------------------------------------------------------------------------!
!  Open Temp/Sal Observation Files for Each Observation Point and Read Data    !
!------------------------------------------------------------------------------!

!----Make Sure Temperature/Salinity Observation File Exists--------------------!
   ONAME = 'enkf_el.dat'
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'ELEVATION OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
!----Open Temp/Salinity Observation File for Read------------------------------!
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_EL
     READ(1,*,IOSTAT=IOS) nsite_tmp,NCNT
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_el.dat at site number:',I
       CALL PSTOP
     ENDIF	
     EL0_OBS(I)%N_TIMES = NCNT
      
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR T OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold EL and Time (TIME)----------!
     ALLOCATE(EL0_OBS(I)%TIMES(EL0_OBS(I)%N_TIMES))
     ALLOCATE(EL0_OBS(I)%EL( EL0_OBS(I)%N_TIMES , EL0_OBS(I)%N_LAYERS ))

!----Read in Current Data for Observation I------------------------------------!
     NLAY = EL0_OBS(I)%N_LAYERS
     DO J=1,EL0_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS)ODAYS,(EL0_OBS(I)%EL(J,K),K=1,NLAY)
       IF(IOS < 0) THEN         ! ios=0 if all goes ok.
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_el.dat at site:',I,  &
                     'Time No:',J
         CALL PSTOP
       ENDIF  
       EL0_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
     END DO
  END DO
  CLOSE(1)

!------------------------------------------------------------------------------!
!  Count Number of Points with Bad Data (TEMP = 0. .OR. SAL = 0.)         
!------------------------------------------------------------------------------!
  NBD_CNT = 0
  DO I=1,N_ASSIM_EL
    DO J=1,EL0_OBS(I)%N_TIMES
      DO K=1,EL0_OBS(I)%N_LAYERS
        IF(EL0_OBS(I)%EL(J,K) < -9990.) THEN
          NBD_CNT = NBD_CNT + 1
        END IF
      END DO
    END DO
  END DO
 
!------------------------------------------------------------------------------!
!  Compute Spatial Interpolation Weights for each Mooring Location 
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  Compute Sigma Layer Weights for Vertical Interpolation            
!------------------------------------------------------------------------------!
    ALLOCATE(ZZ_G(0:MGL,KB))  
   IF(.NOT. ALLOCATED(HG)) ALLOCATE(HG(0:MGL))
   IF(SERIAL)THEN
     ZZ_G = ZZ
     HG = H
   END IF      
#  if defined (MULTIPROCESSOR)   
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,ZZ,ZZ_G)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,H,HG) 
   IF(PAR)CALL MPI_BCAST(ZZ_G,(MGL+1)*KB,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(HG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#  endif

   DO I=1,N_ASSIM_EL
     NLAY = EL0_OBS(I)%N_LAYERS
     ALLOCATE(EL0_OBS(I)%S_INT(NLAY,2))
     ALLOCATE(EL0_OBS(I)%S_WEIGHT(NLAY,2))
     
     X0 = EL0_OBS(I)%X
     Y0 = EL0_OBS(I)%Y
# if defined (SPHERICAL)
     DO NN=1,NGL
       CALL ARC(X0,Y0,XCG(NN),YCG(NN),RDLIST(NN,1))
     END DO
# else
     RDLIST(1:NGL,1) = SQRT((XCG(1:NGL)-X0)**2 + (YCG(1:NGL)-Y0)**2)
# endif
     RDLAST = -1.0_SP
in:  DO WHILE(.TRUE.)
       LOCIJ = MINLOC(RDLIST,RDLIST>RDLAST)
       MIN_LOC = LOCIJ(1)
       IF(MIN_LOC == 0)THEN
         EXIT in
       END IF
       XTRI = XG(NVG(MIN_LOC,1:3))
       YTRI = YG(NVG(MIN_LOC,1:3))
       RDLAST = RDLIST(MIN_LOC,1)
       IF(ISINTRIANGLE1(XTRI,YTRI,X0,Y0))THEN
         JJ = MIN_LOC
	 EXIT IN
       END IF
       RDLAST = RDLIST(MIN_LOC,1)
     END DO IN  	 	      

     ND1=NVG(JJ,1)
     ND2=NVG(JJ,2)
     ND3=NVG(JJ,3)

# if defined (SPHERICAL)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),X0C)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),Y0C)
# else
     X0C = X0 - XCG(JJ)
     Y0C = Y0 - YCG(JJ)
# endif
     !----Linear Interpolation of Bathymetry------------------------------------------!
     H0 = AW0G(JJ,1)*HG(ND1)+AW0G(JJ,2)*HG(ND2)+AW0G(JJ,3)*HG(ND3)
     HX = AWXG(JJ,1)*HG(ND1)+AWXG(JJ,2)*HG(ND2)+AWXG(JJ,3)*HG(ND3)
     HY = AWYG(JJ,1)*HG(ND1)+AWYG(JJ,2)*HG(ND2)+AWYG(JJ,3)*HG(ND3)
     TEMP_H = H0 + HX*X0C + HY*Y0C !here is the interpolated observation depth

     DO K=1,KBM1
      !----Linear Interpolation of Bathymetry------------------------------------------!
      H0 = AW0G(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AW0G(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AW0G(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HX = AWXG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWXG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWXG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HY = AWYG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWYG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWYG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      ZZ_OB(K) = (H0 + HX*X0C + HY*Y0C)/TEMP_H
     END DO  

     DO J=1,NLAY
       SIGMA_C = -EL0_OBS(I)%ODEPTH(J)/TEMP_H         !TS_OBS(I)%DEPTH
       DO K=2,KBM1
         IF(ZZ_OB(K) <= SIGMA_C .AND. ZZ_OB(K-1) > SIGMA_C)THEN 
           EL0_OBS(I)%S_INT(J,1) = K-1
           EL0_OBS(I)%S_INT(J,2) = K
           EL0_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ_OB(K))/(ZZ_OB(K-1)-ZZ_OB(K))
           EL0_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - EL0_OBS(I)%S_WEIGHT(J,1) 
         END IF  
       END DO
       IF(ZZ_OB(1) <= SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
         EL0_OBS(I)%S_INT(J,1) = 1
         EL0_OBS(I)%S_INT(J,2) = 1
         EL0_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         EL0_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       IF(ZZ_OB(KBM1) > SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
         EL0_OBS(I)%S_INT(J,1) = KBM1
         EL0_OBS(I)%S_INT(J,2) = KBM1
         EL0_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         EL0_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF

     END DO
   END DO

   DEALLOCATE(ZZ_G)  

!------------------------------------------------------------------------------!
!  Report Number of Interpolation Points, Location and Number of Data 
!------------------------------------------------------------------------------!
   IF(.NOT. MSR)RETURN

   WRITE(IPT,*)
   WRITE(IPT,*)'!            TEMP/SALINITY OBSERVATION DATA           '
   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_NODE  SITA"
   DO I=1,N_ASSIM_EL
!     MAXEL = MAXLOC(EL0_OBS(I)%X_WEIGHT,DIM=1)
!     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
!     I,EL0_OBS(I)%X/1000.,EL0_OBS(I)%Y/1000., &
!       EL0_OBS(I)%N_INTPTS,EL0_OBS(I)%N_TIMES,EL0_OBS(I)%INTPTS(MAXEL),&
!       EL0_OBS(I)%SITA
   END DO
   WRITE(IPT,*)
   WRITE(IPT,*)'NUMBER OF BAD TS DATA POINTS: ',NBD_CNT
   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
   DO I=1,N_ASSIM_EL
   WRITE(IPT,*)I,EL0_OBS(I)%TIMES(1)/(24.*3600.),&
       EL0_OBS(I)%TIMES(EL0_OBS(I)%N_TIMES)/(24.*3600.)
   END DO

   RETURN
   END SUBROUTINE SET_ENKF_EL_ASSIM_DATA
!==============================================================================|




!==============================================================================!
   SUBROUTINE SET_ENKF_CUR_ASSIM_DATA 
!==============================================================================!

!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR CURRENT OBSERVATIONS                           |
!------------------------------------------------------------------------------!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   INTEGER I,J,K,NN,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   CHARACTER(LEN= 2 ) :: NAC   
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
   REAL(SP), PARAMETER :: ALF = 0.05_SP
   LOGICAL :: FEXIST
   INTEGER :: MAXEL,NBD_CNT
   REAL(SP) :: CUR_RADIUS 
   INTEGER :: JJ

   REAL(SP), DIMENSION(3) :: XTRI,YTRI 
   REAL(SP) :: RDLAST 
   INTEGER :: LOCIJ(2),MIN_LOC,IERR,Nsite_tmp 
   INTEGER :: ND1,ND2,ND3 
   REAL(SP) :: DELTA,COFA,COFB,COFC 
   REAL(SP) ::S11,S22,S33,RTMP,RRTMP 
   REAL(SP), DIMENSION(1:NGL,1) :: RDLIST 
   REAL(SP), DIMENSION(KB) :: ZZ_OB 
   REAL(SP), ALLOCATABLE :: ZZ_G(:,:)   
   REAL(SP) :: X0C,Y0C,HX,HY,H0,TEMP_H
   REAL(DP) :: ODAYS

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: SET_CUR_ASSIM_DATA"

!------------------------------------------------------------------------------!
!  Read Number of Current Observations and Coordinates of Each                 !
!------------------------------------------------------------------------------!
    
   FNAME = "enkf_cur.xy"
!
!--Make Sure Current Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
     
!
!--Read Number of Current Measurement Stations---------------------------------!
!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD') ;rewind(1)
   READ(1,*) N_ASSIM_CUR

   ALLOCATE(CUR_OBS(N_ASSIM_CUR))

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!

   DO I=1,N_ASSIM_CUR
     READ(1,*)ITMP,CUR_OBS(I)%X,CUR_OBS(I)%Y,CUR_OBS(I)%DEPTH,NLAY,CUR_OBS(I)%SITA
     print *, N_ASSIM_CUR,CUR_OBS(I)%X,CUR_OBS(I)%Y,CUR_OBS(I)%DEPTH,NLAY,CUR_OBS(I)%SITA
     CUR_OBS(I)%N_LAYERS = NLAY

     ALLOCATE(CUR_OBS(I)%ODEPTH(NLAY))
     DO J=1,NLAY
       READ(1,*)CUR_OBS(I)%ODEPTH(J)
       IF(CUR_OBS(I)%ODEPTH(J) > CUR_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'CURRENT OBSERVATION LAYER',J,'OF CURRENT MOORING',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH=',CUR_OBS(I)%ODEPTH(J),CUR_OBS(I)%DEPTH
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
     DO J=2,NLAY
       IF(CUR_OBS(I)%ODEPTH(J-1) > CUR_OBS(I)%ODEPTH(J))THEN
         IF(MSR)WRITE(IPT,*)'CURRENT OBSERVATION LAYER',J-1,'OF CURRENT MOORING',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS OBSERVATION LAYER',J,'OF CURRENT MOORING',I,CUR_OBS(I)%ODEPTH(J-1),CUR_OBS(I)%ODEPTH(J)
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO
   CUR_MAX_LAYER = MAXVAL(CUR_OBS(1:N_ASSIM_CUR)%N_LAYERS)

    PRINT *, 'FINISH READING ENKF_CUR.XY'
IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
!
!--Shift Coordinates-----------------------------------------------------------!
!
# if !defined (SPHERICAL) 
   CUR_OBS(:)%X = CUR_OBS(:)%X - VXMIN 
   CUR_OBS(:)%Y = CUR_OBS(:)%Y - VYMIN 
# endif

!# if defined (SPHERICAL) 
   IF(SERIAL)THEN
     XCG = XC
     YCG = YC
     XG = VX
     YG = VY
   END IF
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,XC,XCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,YC,YCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VX,XG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VY,YG)
   IF(PAR)CALL MPI_BCAST(XCG,NGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(YCG,NGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(XG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(YG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#  endif
!# endif
!  


!
!--find cell number (CUR_OBS(J)%N_CELL) of Obs.
!
   RRTMP = 100000.0     !100km
   DO J= 1,N_ASSIM_CUR 
     X0 = CUR_OBS(J)%X
     Y0 = CUR_OBS(J)%Y
     DO I=1,NGL
#    if !defined (SPHERICAL)
       Rtmp = SQRT((XCG(I)-X0)*(XCG(I)-X0)+(YCG(I)-Y0)*(YCG(I)-Y0))
#    else
       CALL ARC(X0,Y0,XCG(I),YCG(I),Rtmp)
#    endif
       IF(Rtmp < RRTMP) THEN
	 S11 = (XG(NVG(I,2))-X0)*(YG(NVG(I,3))-Y0)-&
	       (XG(NVG(I,3))-X0)*(YG(NVG(I,2))-Y0)
         S22 = (XG(NVG(I,3))-X0)*(YG(NVG(I,1))-Y0)-&
	       (XG(NVG(I,1))-X0)*(YG(NVG(I,3))-Y0)
         S33 = (XG(NVG(I,1))-X0)*(YG(NVG(I,2))-Y0)-&
	       (XG(NVG(I,2))-X0)*(YG(NVG(I,1))-Y0)
	 IF(S11 <= 0._SP .AND. S22 <= 0._SP .AND. S33 <= 0._SP) THEN
	   CUR_OBS(J)%N_CELL = I
	   GOTO 300
	 ELSE
	   CUR_OBS(J)%N_CELL = 0
	 ENDIF
       ELSE
	 CUR_OBS(J)%N_CELL = -1
       ENDIF
     ENDDO	
300  CONTINUE
!     IF(MSR) WRITE(200,*) j,X0,Y0,CUR_OBS(J)%N_CELL
     IF(CUR_OBS(J)%N_CELL <= 0) THEN
       IF(MSR) WRITE(IPT,*)'ERROR--CURRENT OBS SITE:',J,' OUT OF DOMAN',&
	       CUR_OBS(J)%N_CELL 
       CALL PSTOP
     ENDIF 	 
   ENDDO      

!--Gather AW0G,AWXG & AWYG use for interp grid to Obs station
   IF(.NOT. ALLOCATED(AW0G))   ALLOCATE(AW0G(0:NGL,3)) 
   IF(.NOT. ALLOCATED(AWXG))   ALLOCATE(AWXG(0:NGL,3))
   IF(.NOT. ALLOCATED(AWYG))   ALLOCATE(AWYG(0:NGL,3))

   IF(SERIAL)THEN
     AW0G = AW0
     AWXG = AWX
     AWYG = AWY
   END IF
# if defined(MULTIPROCESSOR)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AW0,AW0G)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWX,AWXG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWY,AWYG)
   IF(PAR)CALL MPI_BCAST(AW0G,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(AWXG,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(AWYG,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
# endif
!
!--Close Current Observation Global File---------------------------------------!
!
   CLOSE(1)   

!------------------------------------------------------------------------------!
!  Open Current Observation Files for Each Observation Point and Read Data     !
!------------------------------------------------------------------------------!
   ONAME = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)//"_cur.dat"
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_CUR
!----Count Number of Data Entries for Observation I----------------------------!

     READ(1,*,IOSTAT=IOS) Nsite_tmp,NCNT     
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_cur.dat at site number:',I
       CALL PSTOP
     END IF
     CUR_OBS(I)%N_TIMES = NCNT
     
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR CURRENT OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Current (UA,VA) and Time (TIME)-------------------!
     ALLOCATE(CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES))
     ALLOCATE(CUR_OBS(I)%UO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))
     ALLOCATE(CUR_OBS(I)%VO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))

!----Read in Current Data for Observation I------------------------------------!
     NLAY = CUR_OBS(I)%N_LAYERS
     DO J=1,CUR_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS) ODAYS,(CUR_OBS(I)%UO(J,K),CUR_OBS(I)%VO(J,K),K=1,NLAY)
       IF(IOS < 0)THEN         
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_cur.dat at site:',I,'Time No:',J
         CALL PSTOP
       END IF  
       CUR_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
     END DO

!----Convert Time to Seconds---------------------------------------------------!
!----Convert Current Data from cm/s to m/s-------------------------------------!
     CUR_OBS(I)%UO = CUR_OBS(I)%UO * .01_SP
     CUR_OBS(I)%VO = CUR_OBS(I)%VO * .01_SP
  END DO
  CLOSE(1)

!------------------------------------------------------------------------------!
!  Count Number of Points with Bad Data (UO = 0. + V0 = 0.)         
!------------------------------------------------------------------------------!
  NBD_CNT = 0
  DO I=1,N_ASSIM_CUR
    DO J=1,CUR_OBS(I)%N_TIMES
      DO K=1,CUR_OBS(I)%N_LAYERS
        IF(ABS(CUR_OBS(I)%UO(J,K)) + ABS(CUR_OBS(I)%VO(J,K)) < .0001) THEN
          NBD_CNT = NBD_CNT + 1
        END IF
      END DO
    END DO
  END DO
 
!------------------------------------------------------------------------------!
!  Compute Spatial Interpolation Weights for each Mooring Location 
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!  Compute Sigma Layer Weights for Vertical Interpolation                      
!------------------------------------------------------------------------------!
   ALLOCATE(ZZ_G(0:MGL,KB))
   IF(.NOT. ALLOCATED(HG)) ALLOCATE(HG(0:MGL))
   IF(SERIAL)THEN
     ZZ_G = ZZ
     HG = H
   END IF  
#  if defined (MULTIPROCESSOR)   
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,ZZ,ZZ_G) 
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,H,HG) 

   IF(PAR)CALL MPI_BCAST(ZZ_G,(MGL+1)*KB,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(HG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#  endif
   DO I=1,N_ASSIM_CUR
     NLAY = CUR_OBS(I)%N_LAYERS
     ALLOCATE(CUR_OBS(I)%S_INT(NLAY,2))
     ALLOCATE(CUR_OBS(I)%S_WEIGHT(NLAY,2))
     X0 = CUR_OBS(I)%X
     Y0 = CUR_OBS(I)%Y
# if defined (SPHERICAL)
     DO NN=1,NGL
       CALL ARC(X0,Y0,XCG(NN),YCG(NN),RDLIST(NN,1))
     END DO
# else
     RDLIST(1:NGL,1) = SQRT((XCG(1:NGL)-X0)**2 + (YCG(1:NGL)-Y0)**2)
# endif
     RDLAST = -1.0_SP
in:  DO WHILE(.TRUE.)
       LOCIJ = MINLOC(RDLIST,RDLIST>RDLAST)
       MIN_LOC = LOCIJ(1)
       IF(MIN_LOC == 0)THEN
         EXIT in
       END IF
       XTRI = XG(NVG(MIN_LOC,1:3))
       YTRI = YG(NVG(MIN_LOC,1:3))
       RDLAST = RDLIST(MIN_LOC,1)
       IF(ISINTRIANGLE1(XTRI,YTRI,X0,Y0))THEN
         JJ = MIN_LOC
	 EXIT IN
       END IF
       RDLAST = RDLIST(MIN_LOC,1)
     END DO IN  	 	      

     ND1=NVG(JJ,1)
     ND2=NVG(JJ,2)
     ND3=NVG(JJ,3)

# if defined (SPHERICAL)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),X0C)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),Y0C)
# else
     X0C = X0 - XCG(JJ)
     Y0C = Y0 - YCG(JJ)
# endif
     !----Linear Interpolation of Bathymetry------------------------------------------!
     H0 = AW0G(JJ,1)*HG(ND1)+AW0G(JJ,2)*HG(ND2)+AW0G(JJ,3)*HG(ND3)
     HX = AWXG(JJ,1)*HG(ND1)+AWXG(JJ,2)*HG(ND2)+AWXG(JJ,3)*HG(ND3)
     HY = AWYG(JJ,1)*HG(ND1)+AWYG(JJ,2)*HG(ND2)+AWYG(JJ,3)*HG(ND3)
     TEMP_H = H0 + HX*X0C + HY*Y0C !here is the interpolated observation depth

     DO K=1,KBM1
      !----Linear Interpolation of Bathymetry------------------------------------------!
      H0 = AW0G(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AW0G(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AW0G(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HX = AWXG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWXG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWXG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HY = AWYG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWYG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWYG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      ZZ_OB(K) = (H0 + HX*X0C + HY*Y0C)/TEMP_H
     END DO  

     DO J=1,NLAY
       SIGMA_C = -CUR_OBS(I)%ODEPTH(J)/TEMP_H
       DO K=2,KBM1
         IF(ZZ_OB(K) <= SIGMA_C .AND. ZZ_OB(K-1) > SIGMA_C)THEN 
           CUR_OBS(I)%S_INT(J,1) = K-1
           CUR_OBS(I)%S_INT(J,2) = K
           CUR_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ_OB(K))/(ZZ_OB(K-1)-ZZ_OB(K))
           CUR_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - CUR_OBS(I)%S_WEIGHT(J,1) 
         END IF  
       END DO
       IF(ZZ_OB(1) <= SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
         CUR_OBS(I)%S_INT(J,1) = 1
         CUR_OBS(I)%S_INT(J,2) = 1
         CUR_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         CUR_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       IF(ZZ_OB(KBM1) >= SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
         CUR_OBS(I)%S_INT(J,1) = KBM1
         CUR_OBS(I)%S_INT(J,2) = KBM1
         CUR_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         CUR_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       if (msr) print *, "obs",i,"lay", j,SIGMA_C,CUR_OBS(I)%S_INT(J,1),CUR_OBS(I)%S_INT(J,2), CUR_OBS(I)%ODEPTH(J),CUR_OBS(I)%DEPTH,'scoor',zz_ob(1:kbm1)
     END DO
   END DO
!#  if defined (MULTIPROCESSOR)   
   DEALLOCATE(ZZ_G)  
!#  endif
!------------------------------------------------------------------------------!
!  Report Number of Interpolation Points, Location and Number of Data 
!------------------------------------------------------------------------------!
   IF(.NOT. MSR)RETURN

   WRITE(IPT,*)
   WRITE(IPT,*)'!            CURRENT OBSERVATION DATA           '
# if defined (SPHERICAL)
   WRITE(IPT,*)" MOORING#    LON       LAT    #INTERP PTS  #DATA TIMES  NEAR_CELL  SITA"
# else
   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_CELL  SITA"
# endif
   DO I=1,N_ASSIM_CUR
!     MAXEL = MAXLOC(CUR_OBS(I)%X_WEIGHT,DIM=1)
!     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
!# if defined (SPHERICAL)
!     I,CUR_OBS(I)%X,CUR_OBS(I)%Y, &
!# else
!     I,(CUR_OBS(I)%X+VXMIN)/1000.,(CUR_OBS(I)%Y+VYMIN)/1000., &
!# endif
!       CUR_OBS(I)%N_INTPTS,CUR_OBS(I)%N_TIMES,CUR_OBS(I)%INTPTS(MAXEL),&
!       CUR_OBS(I)%SITA
   END DO
   WRITE(IPT,*)
   WRITE(IPT,*)'NUMBER OF BAD CURRENT DATA POINTS: ',NBD_CNT
   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
   DO I=1,N_ASSIM_CUR
    WRITE(IPT,*)I, CUR_OBS(I)%N_CELL,CUR_OBS(I)%innode,CUR_OBS(I)%s_weight
!   WRITE(IPT,*)I,CUR_OBS(I)%TIMES(1)/(24.*3600.),&
!       CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES)/(24.*3600.)
   WRITE(IPT,*)I,CUR_OBS(I)%TIMES(1)%MJD,&
       CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES)%MJD
   END DO
   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SET_CUR_ASSIM_DATA"

 
   RETURN
   END SUBROUTINE SET_ENKF_CUR_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
   SUBROUTINE SET_ENKF_T_ASSIM_DATA 
!
!!------------------------------------------------------------------------------!
!!  SET UP ASSIMILATION DATA FOR TEMP/SAL OBSERVATIONS                          |
!!------------------------------------------------------------------------------!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)   
   USE MOD_PAR  
# endif
   IMPLICIT NONE
   INTEGER I,J,K,NN,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   CHARACTER(LEN= 2 ) :: NAC   
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
   REAL(SP), PARAMETER :: ALF = 0.05_SP
   LOGICAL :: FEXIST
   INTEGER :: MAXEL,NBD_CNT
   REAL(SP) :: TS_RADIUS 
   INTEGER :: JMIN,INTMP(1) 
   REAL(SP):: LMIN

   REAL(SP), DIMENSION(1:NGL,1) :: RDLIST
   REAL(SP), DIMENSION(3) :: XTRI,YTRI
   REAL(SP) :: RDLAST
   INTEGER :: LOCIJ(2),MIN_LOC,JJ,IERR,Nsite_tmp
   INTEGER :: ND1,ND2,ND3
   REAL(SP) :: DELTA,COFA,COFB,COFC
   REAL(SP) ::S11,S22,S33,RTMP,RRTMP
   REAL(SP), DIMENSION(KB) :: ZZ_OB
!#  if defined (MULTIPROCESSOR)   
   REAL(SP), ALLOCATABLE :: ZZ_G(:,:)  
!# endif
   REAL(SP) :: X0C,Y0C,HX,HY,H0,TEMP_H
   REAL(DP) :: ODAYS
      
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_TS_ASSIM_DATA" 

!------------------------------------------------------------------------------!
!  Read Number of Scalar Observations and Coordinates of Each                  !
!------------------------------------------------------------------------------!
       
   FNAME = "enkf_t.xy"
!
!--Make Sure Temperature and Salinity Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'TEMP OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
     
!
!--Read Number of T/S Measurement Stations---------------------------------!
!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_T                  !nomber of TS Obs station
   ALLOCATE(T0_OBS(N_ASSIM_T))          !Type for TS_OBS

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!
   DO I=1,N_ASSIM_T  
     READ(1,*)ITMP,T0_OBS(I)%X,T0_OBS(I)%Y,T0_OBS(I)%DEPTH,NLAY,T0_OBS(I)%SITA
     T0_OBS(I)%N_LAYERS = NLAY
     ALLOCATE(T0_OBS(I)%ODEPTH(NLAY))
     DO J=1,NLAY
       READ(1,*)T0_OBS(I)%ODEPTH(J)
       IF(T0_OBS(I)%ODEPTH(J) > T0_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'OBSERVATION DEPTH',J,'OF TEMP OBS',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH=',T0_OBS(I)%ODEPTH(J),T0_OBS(I)%DEPTH
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO

!
!--Shift Coordinates-----------------------------------------------------------!
!
# if !defined (SPHERICAL) 
   T0_OBS(:)%X = T0_OBS(:)%X - VXMIN 
   T0_OBS(:)%Y = T0_OBS(:)%Y - VYMIN 
# endif
   
!# if defined (SPHERICAL) 
   IF(SERIAL)THEN
     XCG = XC
     YCG = YC
     XG = VX
     YG = VY
   END IF
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,XC,XCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,YC,YCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VX,XG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VY,YG)
   IF(PAR)CALL MPI_BCAST(XCG,NGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(YCG,NGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(XG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(YG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#  endif
!# endif
!  
!
!--find the cell number (TS_OBS(:)%N_CELL) of Obs station---------------------
!
   RRTMP = 100000.0     !100km
   DO J= 1,N_ASSIM_T 
     X0 = T0_OBS(J)%X
     Y0 = T0_OBS(J)%Y
     DO I=1,NGL
#    if !defined (SPHERICAL)
       Rtmp = SQRT((XCG(I)-X0)*(XCG(I)-X0)+(YCG(I)-Y0)*(YCG(I)-Y0))
#    else
       CALL ARC(X0,Y0,XCG(I),YCG(I),Rtmp)
#    endif
       IF(Rtmp < RRTMP) THEN
	 S11 = (XG(NVG(I,2))-X0)*(YG(NVG(I,3))-Y0)-&
	       (XG(NVG(I,3))-X0)*(YG(NVG(I,2))-Y0)
         S22 = (XG(NVG(I,3))-X0)*(YG(NVG(I,1))-Y0)-&
	       (XG(NVG(I,1))-X0)*(YG(NVG(I,3))-Y0)
         S33 = (XG(NVG(I,1))-X0)*(YG(NVG(I,2))-Y0)-&
	       (XG(NVG(I,2))-X0)*(YG(NVG(I,1))-Y0)
	 IF(S11 <= 0._SP .AND. S22 <= 0._SP .AND. S33 <= 0._SP) THEN
	   T0_OBS(J)%N_CELL = I
           S11=SQRT((XG(NVG(I,1))-X0)**2+(YG(NVG(I,1))-Y0)**2) ! CACLULATE DISTANCE
           S22=SQRT((XG(NVG(I,2))-X0)**2+(YG(NVG(I,2))-Y0)**2) 
           S33=SQRT((XG(NVG(I,3))-X0)**2+(YG(NVG(I,3))-Y0)**2)
           intmp=minloc([s11,s22,s33])
           T0_OBS(J)%innode=NVG(I,intmp(1))
	   GOTO 300
	 ELSE
	   T0_OBS(J)%N_CELL = 0
	 ENDIF
       ELSE
	 T0_OBS(J)%N_CELL = -1
       ENDIF
     ENDDO	
300  CONTINUE
     IF(T0_OBS(J)%N_CELL <= 0) THEN
       IF(MSR) WRITE(IPT,*)'ERROR--T/S OBS SITE:',J,' OUT OF DOMAN',&
	                    T0_OBS(J)%N_CELL 
       CALL PSTOP
     ENDIF 	 

   ENDDO             
!--Gather AW0G,AWXG & AWYG use for interp grid to Obs station
   IF(.NOT. ALLOCATED(AW0G)) ALLOCATE(AW0G(0:NGL,3))
   IF(.NOT. ALLOCATED(AWXG)) ALLOCATE(AWXG(0:NGL,3))
   IF(.NOT. ALLOCATED(AWYG)) ALLOCATE(AWYG(0:NGL,3))
   
   IF(SERIAL)THEN
     AW0G = AW0
     AWXG = AWX
     AWYG = AWY
   END IF
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN     
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AW0,AW0G)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWX,AWXG)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWY,AWYG)
     CALL MPI_BCAST(AW0G,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
     CALL MPI_BCAST(AWXG,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
     CALL MPI_BCAST(AWYG,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   END IF
#  endif     
!
!--Close Current Observation Global File---------------------------------------!
!
   CLOSE(1)
       
!------------------------------------------------------------------------------!
!  Open Temp/Sal Observation Files for Each Observation Point and Read Data    !
!------------------------------------------------------------------------------!

!----Make Sure Temperature/Salinity Observation File Exists--------------------!
   ONAME = 'enkf_t.dat'
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'TEMP/SALINITY OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
!----Open Temp/Salinity Observation File for Read------------------------------!
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_T
     READ(1,*,IOSTAT=IOS) nsite_tmp,NCNT
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_t.dat at site number:',I
       CALL PSTOP
     ENDIF	
     T0_OBS(I)%N_TIMES = NCNT
      
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR T OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Temp/Salinity (TEMP/SAL) and Time (TIME)----------!
     ALLOCATE(T0_OBS(I)%TIMES(T0_OBS(I)%N_TIMES))
     ALLOCATE(T0_OBS(I)%TEMP( T0_OBS(I)%N_TIMES , T0_OBS(I)%N_LAYERS ))

!----Read in Current Data for Observation I------------------------------------!
     NLAY = T0_OBS(I)%N_LAYERS
     DO J=1,T0_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS)ODAYS,(T0_OBS(I)%TEMP(J,K),K=1,NLAY)
       IF(IOS < 0) THEN         ! ios=0 if all goes ok.
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_t.dat at site:',I,  &
                     'Time No:',J
         CALL PSTOP
       ENDIF  
       T0_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
     END DO
  END DO
  CLOSE(1)

!------------------------------------------------------------------------------!
!  Count Number of Points with Bad Data (TEMP = 0. .OR. SAL = 0.)         
!------------------------------------------------------------------------------!
  NBD_CNT = 0
  DO I=1,N_ASSIM_T
    DO J=1,T0_OBS(I)%N_TIMES
      DO K=1,T0_OBS(I)%N_LAYERS
        IF(T0_OBS(I)%TEMP(J,K) < -90.) THEN
          NBD_CNT = NBD_CNT + 1
        END IF
      END DO
    END DO
  END DO
 
!------------------------------------------------------------------------------!
!  Compute Spatial Interpolation Weights for each Mooring Location 
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  Compute Sigma Layer Weights for Vertical Interpolation            
!------------------------------------------------------------------------------!
   ALLOCATE(ZZ_G(0:MGL,KB))  
   IF(.NOT. ALLOCATED(HG)) ALLOCATE(HG(0:MGL))
   IF(SERIAL)THEN
     ZZ_G = ZZ
     HG = H
   END IF      
#  if defined (MULTIPROCESSOR)   
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,ZZ,ZZ_G)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,H,HG) 
   IF(PAR)CALL MPI_BCAST(ZZ_G,(MGL+1)*KB,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(HG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#  endif

   DO I=1,N_ASSIM_T
     NLAY = T0_OBS(I)%N_LAYERS
     ALLOCATE(T0_OBS(I)%S_INT(NLAY,2))
     ALLOCATE(T0_OBS(I)%S_WEIGHT(NLAY,2))
     
     X0 = T0_OBS(I)%X
     Y0 = T0_OBS(I)%Y
# if defined (SPHERICAL)
     DO NN=1,NGL
       CALL ARC(X0,Y0,XCG(NN),YCG(NN),RDLIST(NN,1))
     END DO
# else
     RDLIST(1:NGL,1) = SQRT((XCG(1:NGL)-X0)**2 + (YCG(1:NGL)-Y0)**2)
# endif
     RDLAST = -1.0_SP
in:  DO WHILE(.TRUE.)
       LOCIJ = MINLOC(RDLIST,RDLIST>RDLAST)
       MIN_LOC = LOCIJ(1)
       IF(MIN_LOC == 0)THEN
         EXIT in
       END IF
       XTRI = XG(NVG(MIN_LOC,1:3))
       YTRI = YG(NVG(MIN_LOC,1:3))
       RDLAST = RDLIST(MIN_LOC,1)
       IF(ISINTRIANGLE1(XTRI,YTRI,X0,Y0))THEN
         JJ = MIN_LOC
	 EXIT IN
       END IF
       RDLAST = RDLIST(MIN_LOC,1)
     END DO IN  	 	      

     ND1=NVG(JJ,1)
     ND2=NVG(JJ,2)
     ND3=NVG(JJ,3)

# if defined (SPHERICAL)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),X0C)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),Y0C)
# else
     X0C = X0 - XCG(JJ)
     Y0C = Y0 - YCG(JJ)
# endif
     !----Linear Interpolation of Bathymetry------------------------------------------!
     H0 = AW0G(JJ,1)*HG(ND1)+AW0G(JJ,2)*HG(ND2)+AW0G(JJ,3)*HG(ND3)
     HX = AWXG(JJ,1)*HG(ND1)+AWXG(JJ,2)*HG(ND2)+AWXG(JJ,3)*HG(ND3)
     HY = AWYG(JJ,1)*HG(ND1)+AWYG(JJ,2)*HG(ND2)+AWYG(JJ,3)*HG(ND3)
     TEMP_H = H0 + HX*X0C + HY*Y0C !here is the interpolated observation depth

     DO K=1,KBM1
      !----Linear Interpolation of Bathymetry------------------------------------------!
      H0 = AW0G(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AW0G(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AW0G(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HX = AWXG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWXG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWXG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HY = AWYG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWYG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWYG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      ZZ_OB(K) = (H0 + HX*X0C + HY*Y0C)/TEMP_H
     END DO  

     DO J=1,NLAY
       SIGMA_C = -T0_OBS(I)%ODEPTH(J)/TEMP_H         !TS_OBS(I)%DEPTH
       DO K=2,KBM1
         IF(ZZ_OB(K) <= SIGMA_C .AND. ZZ_OB(K-1) > SIGMA_C)THEN 
           T0_OBS(I)%S_INT(J,1) = K-1
           T0_OBS(I)%S_INT(J,2) = K
           T0_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ_OB(K))/(ZZ_OB(K-1)-ZZ_OB(K))
           T0_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - T0_OBS(I)%S_WEIGHT(J,1) 
         END IF  
       END DO
       IF(ZZ_OB(1) <= SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
         T0_OBS(I)%S_INT(J,1) = 1
         T0_OBS(I)%S_INT(J,2) = 1
         T0_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         T0_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       IF(ZZ_OB(KBM1) > SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
         T0_OBS(I)%S_INT(J,1) = KBM1
         T0_OBS(I)%S_INT(J,2) = KBM1
         T0_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         T0_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF

     END DO
   END DO

   DEALLOCATE(ZZ_G)  

!------------------------------------------------------------------------------!
!  Report Number of Interpolation Points, Location and Number of Data 
!------------------------------------------------------------------------------!
   IF(.NOT. MSR)RETURN

   WRITE(IPT,*)
   WRITE(IPT,*)'!            TEMP/SALINITY OBSERVATION DATA           '
   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_NODE  SITA"
   DO I=1,N_ASSIM_T
!     MAXEL = MAXLOC(T0_OBS(I)%X_WEIGHT,DIM=1)
!     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
!     I,T0_OBS(I)%X/1000.,T0_OBS(I)%Y/1000., &
!       T0_OBS(I)%N_INTPTS,T0_OBS(I)%N_TIMES,T0_OBS(I)%INTPTS(MAXEL),&
!       T0_OBS(I)%SITA
   END DO
   WRITE(IPT,*)
   WRITE(IPT,*)'NUMBER OF BAD TS DATA POINTS: ',NBD_CNT
   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
   DO I=1,N_ASSIM_T
   WRITE(IPT,*)I,T0_OBS(I)%TIMES(1)/(24.*3600.),&
       T0_OBS(I)%TIMES(T0_OBS(I)%N_TIMES)/(24.*3600.)
   END DO

   RETURN
   END SUBROUTINE SET_ENKF_T_ASSIM_DATA
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
   SUBROUTINE SET_ENKF_S_ASSIM_DATA 
!
!!------------------------------------------------------------------------------!
!!  SET UP ASSIMILATION DATA FOR TEMP/SAL OBSERVATIONS                          |
!!------------------------------------------------------------------------------!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)   
   USE MOD_PAR  
# endif
   IMPLICIT NONE
   INTEGER I,J,K,NN,ECNT,ITMP,NCNT,IOS,NLAY
   CHARACTER(LEN=120) :: FNAME,ONAME
   CHARACTER(LEN= 2 ) :: NAC   
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
   REAL(SP), PARAMETER :: ALF = 0.05_SP
   LOGICAL :: FEXIST
   INTEGER :: MAXEL,NBD_CNT
   REAL(SP) :: TS_RADIUS 
   INTEGER :: JMIN,INTMP(1) 
   REAL(SP):: LMIN

   REAL(SP), DIMENSION(1:NGL,1) :: RDLIST
   REAL(SP), DIMENSION(3) :: XTRI,YTRI
   REAL(SP) :: RDLAST
   INTEGER :: LOCIJ(2),MIN_LOC,JJ,IERR,Nsite_tmp
   INTEGER :: ND1,ND2,ND3
   REAL(SP) :: DELTA,COFA,COFB,COFC
   REAL(SP) ::S11,S22,S33,RTMP,RRTMP
   REAL(SP), DIMENSION(KB) :: ZZ_OB
!#  if defined (MULTIPROCESSOR)   
   REAL(SP), ALLOCATABLE :: ZZ_G(:,:)  
!# endif
   REAL(SP) :: X0C,Y0C,HX,HY,H0,TEMP_H
   REAL(DP) :: ODAYS
      
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_TS_ASSIM_DATA" 

!------------------------------------------------------------------------------!
!  Read Number of Scalar Observations and Coordinates of Each                  !
!------------------------------------------------------------------------------!
       
   FNAME = "enkf_s.xy"
!
!--Make Sure Temperature and Salinity Assimilation Data File Exists-----------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'Salinity OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
     
!
!--Read Number of T/S Measurement Stations---------------------------------!
!
   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
   READ(1,*) N_ASSIM_S                  !nomber of TS Obs station
   ALLOCATE(S0_OBS(N_ASSIM_S))          !Type for TS_OBS

!
!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!
   DO I=1,N_ASSIM_S  
     READ(1,*)ITMP,S0_OBS(I)%X,S0_OBS(I)%Y,S0_OBS(I)%DEPTH,NLAY,S0_OBS(I)%SITA
     S0_OBS(I)%N_LAYERS = NLAY
     ALLOCATE(S0_OBS(I)%ODEPTH(NLAY))
     DO J=1,NLAY
       READ(1,*)S0_OBS(I)%ODEPTH(J)
       IF(S0_OBS(I)%ODEPTH(J) > S0_OBS(I)%DEPTH)THEN
         IF(MSR)WRITE(IPT,*)'OBSERVATION DEPTH',J,'OF salinity OBS',I
         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH=',S0_OBS(I)%ODEPTH(J),S0_OBS(I)%DEPTH
         IF(MSR)WRITE(IPT,*)'HALTING...........'
         CALL PSTOP
       END IF
     END DO
   END DO

!
!--Shift Coordinates-----------------------------------------------------------!
!
# if !defined (SPHERICAL) 
   S0_OBS(:)%X = S0_OBS(:)%X - VXMIN 
   S0_OBS(:)%Y = S0_OBS(:)%Y - VYMIN 
# endif
   
!# if defined (SPHERICAL) 
   IF(SERIAL)THEN
     XCG = XC
     YCG = YC
     XG = VX
     YG = VY
   END IF
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,XC,XCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,YC,YCG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VX,XG)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VY,YG)
   IF(PAR)CALL MPI_BCAST(XCG,NGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(YCG,NGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(XG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(YG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#  endif
!# endif
!  
!
!--find the cell number (TS_OBS(:)%N_CELL) of Obs station---------------------
!
   RRTMP = 100000.0     !100km
   DO J= 1,N_ASSIM_S 
     X0 = S0_OBS(J)%X
     Y0 = S0_OBS(J)%Y
     DO I=1,NGL
#    if !defined (SPHERICAL)
       Rtmp = SQRT((XCG(I)-X0)*(XCG(I)-X0)+(YCG(I)-Y0)*(YCG(I)-Y0))
#    else
       CALL ARC(X0,Y0,XCG(I),YCG(I),Rtmp)
#    endif
       IF(Rtmp < RRTMP) THEN
	 S11 = (XG(NVG(I,2))-X0)*(YG(NVG(I,3))-Y0)-&
	       (XG(NVG(I,3))-X0)*(YG(NVG(I,2))-Y0)
         S22 = (XG(NVG(I,3))-X0)*(YG(NVG(I,1))-Y0)-&
	       (XG(NVG(I,1))-X0)*(YG(NVG(I,3))-Y0)
         S33 = (XG(NVG(I,1))-X0)*(YG(NVG(I,2))-Y0)-&
	       (XG(NVG(I,2))-X0)*(YG(NVG(I,1))-Y0)
	 IF(S11 <= 0._SP .AND. S22 <= 0._SP .AND. S33 <= 0._SP) THEN
	   S0_OBS(J)%N_CELL = I
           S11=SQRT((XG(NVG(I,1))-X0)**2+(YG(NVG(I,1))-Y0)**2) ! CACLULATE DISTANCE
           S22=SQRT((XG(NVG(I,2))-X0)**2+(YG(NVG(I,2))-Y0)**2) 
           S33=SQRT((XG(NVG(I,3))-X0)**2+(YG(NVG(I,3))-Y0)**2)
           intmp=minloc([s11,s22,s33])
           S0_OBS(J)%innode=NVG(I,intmp(1))
	   GOTO 300
	 ELSE
	   S0_OBS(J)%N_CELL = 0
	 ENDIF
       ELSE
	 S0_OBS(J)%N_CELL = -1
       ENDIF
     ENDDO	
300  CONTINUE
     IF(S0_OBS(J)%N_CELL <= 0) THEN
       IF(MSR) WRITE(IPT,*)'ERROR--S OBS SITE:',J,' OUT OF DOMAN',&
	                    S0_OBS(J)%N_CELL 
       CALL PSTOP
     ENDIF 	 

   ENDDO             
!--Gather AW0G,AWXG & AWYG use for interp grid to Obs station
   IF(.NOT. ALLOCATED(AW0G)) ALLOCATE(AW0G(0:NGL,3))
   IF(.NOT. ALLOCATED(AWXG)) ALLOCATE(AWXG(0:NGL,3))
   IF(.NOT. ALLOCATED(AWYG)) ALLOCATE(AWYG(0:NGL,3))
   
   IF(SERIAL)THEN
     AW0G = AW0
     AWXG = AWX
     AWYG = AWY
   END IF
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN     
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AW0,AW0G)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWX,AWXG)
     CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,AWY,AWYG)
     CALL MPI_BCAST(AW0G,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
     CALL MPI_BCAST(AWXG,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
     CALL MPI_BCAST(AWYG,(NGL+1)*3,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   END IF
#  endif     
!
!--Close Current Observation Global File---------------------------------------!
!
   CLOSE(1)
       
!------------------------------------------------------------------------------!
!  Open Temp/Sal Observation Files for Each Observation Point and Read Data    !
!------------------------------------------------------------------------------!

!----Make Sure Temperature/Salinity Observation File Exists--------------------!
   ONAME = 'enkf_s.dat'
   INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'TEMP/SALINITY OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
!----Open Temp/Salinity Observation File for Read------------------------------!
   OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)

   DO I=1,N_ASSIM_S
     READ(1,*,IOSTAT=IOS) nsite_tmp,NCNT
     IF(IOS<0) then
       WRITE(IPT,*) 'ERROR in read ',trim(casename),'_s.dat at site number:',I
       CALL PSTOP
     ENDIF	
     S0_OBS(I)%N_TIMES = NCNT
      
     IF(NCNT == 0)THEN
       IF(MSR)WRITE(IPT,*)'NO DATA FOR S OBSERVATION',I
       CALL PSTOP
     END IF

!----Allocate Arrays to Hold Temp/Salinity (TEMP/SAL) and Time (TIME)----------!
     ALLOCATE(S0_OBS(I)%TIMES(S0_OBS(I)%N_TIMES))
     ALLOCATE(S0_OBS(I)%SAL( S0_OBS(I)%N_TIMES , S0_OBS(I)%N_LAYERS ))

!----Read in Current Data for Observation I------------------------------------!
     NLAY = S0_OBS(I)%N_LAYERS
     DO J=1,S0_OBS(I)%N_TIMES
       READ(1,*,IOSTAT=IOS)ODAYS,(S0_OBS(I)%SAL(J,K),K=1,NLAY)
       IF(IOS < 0) THEN         ! ios=0 if all goes ok.
         WRITE(IPT,*)'ERROR in read ',trim(casename),'_t.dat at site:',I,  &
                     'Time No:',J
         CALL PSTOP
       ENDIF  
       S0_OBS(I)%TIMES(J) = DAYS2TIME(ODAYS)
     END DO
  END DO
  CLOSE(1)

!------------------------------------------------------------------------------!
!  Count Number of Points with Bad Data (TEMP = 0. .OR. SAL = 0.)         
!------------------------------------------------------------------------------!
  NBD_CNT = 0
  DO I=1,N_ASSIM_S
    DO J=1,S0_OBS(I)%N_TIMES
      DO K=1,S0_OBS(I)%N_LAYERS
        IF(S0_OBS(I)%SAL(J,K) < -90.) THEN
          NBD_CNT = NBD_CNT + 1
        END IF
      END DO
    END DO
  END DO
 
!------------------------------------------------------------------------------!
!  Compute Spatial Interpolation Weights for each Mooring Location 
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  Compute Sigma Layer Weights for Vertical Interpolation            
!------------------------------------------------------------------------------!
   ALLOCATE(ZZ_G(0:MGL,KB))  
   IF(.NOT. ALLOCATED(HG)) ALLOCATE(HG(0:MGL))
   IF(SERIAL)THEN
     ZZ_G = ZZ
     HG = H
   END IF      
#  if defined (MULTIPROCESSOR)   
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,ZZ,ZZ_G)
   IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,H,HG) 
   IF(PAR)CALL MPI_BCAST(ZZ_G,(MGL+1)*KB,MPI_F,0,MPI_FVCOM_GROUP,IERR)
   IF(PAR)CALL MPI_BCAST(HG,MGL+1,MPI_F,0,MPI_FVCOM_GROUP,IERR)
#  endif

   DO I=1,N_ASSIM_S
     NLAY = S0_OBS(I)%N_LAYERS
     ALLOCATE(S0_OBS(I)%S_INT(NLAY,2))
     ALLOCATE(S0_OBS(I)%S_WEIGHT(NLAY,2))
     
     X0 = S0_OBS(I)%X
     Y0 = S0_OBS(I)%Y
# if defined (SPHERICAL)
     DO NN=1,NGL
       CALL ARC(X0,Y0,XCG(NN),YCG(NN),RDLIST(NN,1))
     END DO
# else
     RDLIST(1:NGL,1) = SQRT((XCG(1:NGL)-X0)**2 + (YCG(1:NGL)-Y0)**2)
# endif
     RDLAST = -1.0_SP
in:  DO WHILE(.TRUE.)
       LOCIJ = MINLOC(RDLIST,RDLIST>RDLAST)
       MIN_LOC = LOCIJ(1)
       IF(MIN_LOC == 0)THEN
         EXIT in
       END IF
       XTRI = XG(NVG(MIN_LOC,1:3))
       YTRI = YG(NVG(MIN_LOC,1:3))
       RDLAST = RDLIST(MIN_LOC,1)
       IF(ISINTRIANGLE1(XTRI,YTRI,X0,Y0))THEN
         JJ = MIN_LOC
	 EXIT IN
       END IF
       RDLAST = RDLIST(MIN_LOC,1)
     END DO IN  	 	      

     ND1=NVG(JJ,1)
     ND2=NVG(JJ,2)
     ND3=NVG(JJ,3)

# if defined (SPHERICAL)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),X0C)
     CALL ARCX(X0,Y0,XCG(JJ),YCG(JJ),Y0C)
# else
     X0C = X0 - XCG(JJ)
     Y0C = Y0 - YCG(JJ)
# endif
     !----Linear Interpolation of Bathymetry------------------------------------------!
     H0 = AW0G(JJ,1)*HG(ND1)+AW0G(JJ,2)*HG(ND2)+AW0G(JJ,3)*HG(ND3)
     HX = AWXG(JJ,1)*HG(ND1)+AWXG(JJ,2)*HG(ND2)+AWXG(JJ,3)*HG(ND3)
     HY = AWYG(JJ,1)*HG(ND1)+AWYG(JJ,2)*HG(ND2)+AWYG(JJ,3)*HG(ND3)
     TEMP_H = H0 + HX*X0C + HY*Y0C !here is the interpolated observation depth

     DO K=1,KBM1
      !----Linear Interpolation of Bathymetry------------------------------------------!
      H0 = AW0G(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AW0G(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AW0G(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HX = AWXG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWXG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWXG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      HY = AWYG(JJ,1)*ZZ_G(ND1,K)*HG(ND1)+AWYG(JJ,2)*ZZ_G(ND2,K)*HG(ND2)+AWYG(JJ,3)*ZZ_G(ND3,K)*HG(ND3)
      ZZ_OB(K) = (H0 + HX*X0C + HY*Y0C)/TEMP_H
     END DO  

     DO J=1,NLAY
       SIGMA_C = -S0_OBS(I)%ODEPTH(J)/TEMP_H         !TS_OBS(I)%DEPTH
       DO K=2,KBM1
         IF(ZZ_OB(K) <= SIGMA_C .AND. ZZ_OB(K-1) > SIGMA_C)THEN 
           S0_OBS(I)%S_INT(J,1) = K-1
           S0_OBS(I)%S_INT(J,2) = K
           S0_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ_OB(K))/(ZZ_OB(K-1)-ZZ_OB(K))
           S0_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - S0_OBS(I)%S_WEIGHT(J,1) 
         END IF  
       END DO
       IF(ZZ_OB(1) <= SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
         S0_OBS(I)%S_INT(J,1) = 1
         S0_OBS(I)%S_INT(J,2) = 1
         S0_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         S0_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF
       IF(ZZ_OB(KBM1) > SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
         S0_OBS(I)%S_INT(J,1) = KBM1
         S0_OBS(I)%S_INT(J,2) = KBM1
         S0_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
         S0_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
       END IF

     END DO
   END DO

   DEALLOCATE(ZZ_G)  

!------------------------------------------------------------------------------!
!  Report Number of Interpolation Points, Location and Number of Data 
!------------------------------------------------------------------------------!
   IF(.NOT. MSR)RETURN

   WRITE(IPT,*)
   WRITE(IPT,*)'!            TEMP/SALINITY OBSERVATION DATA           '
   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_NODE  SITA"
   DO I=1,N_ASSIM_S
!     MAXEL = MAXLOC(S0_OBS(I)%X_WEIGHT,DIM=1)
!     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
!     I,S0_OBS(I)%X/1000.,S0_OBS(I)%Y/1000., &
!       S0_OBS(I)%N_INTPTS,S0_OBS(I)%N_TIMES,S0_OBS(I)%INTPTS(MAXEL),&
!       S0_OBS(I)%SITA
   END DO
   WRITE(IPT,*)
   WRITE(IPT,*)'NUMBER OF BAD TS DATA POINTS: ',NBD_CNT
   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
   DO I=1,N_ASSIM_S
   WRITE(IPT,*)I,S0_OBS(I)%TIMES(1)/(24.*3600.),&
       S0_OBS(I)%TIMES(S0_OBS(I)%N_TIMES)/(24.*3600.)
   END DO

   RETURN
   END SUBROUTINE SET_ENKF_S_ASSIM_DATA
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|


!!=============================================================
!!==============================================================================|
   LOGICAL FUNCTION ISINTRIANGLE1(XT,YT,X0,Y0) 
!!==============================================================================|
!  determine if point (x0,y0) is in triangle defined by nodes (xt(3),yt(3))    |
!  using algorithm used for scene rendering in computer graphics               |
!  algorithm works well unless particle happens to lie in a line parallel      |
!  to the edge of a triangle.                                                  |
!  This can cause problems if you use a regular grid, say for idealized        |
!  modelling and you happen to see particles right on edges or parallel to     |
!  edges.                                                                      |

   USE MOD_PREC
   IMPLICIT NONE
   REAL(SP), INTENT(IN) :: X0,Y0
   REAL(SP), INTENT(IN) :: XT(3),YT(3)
   REAL(SP) :: F1,F2,F3
   REAL(SP) :: X1(2)
   REAL(SP) :: X2(2)
   REAL(SP) :: X3(2)
   REAL(SP) :: P(2)

!------------------------------------------------------------------------------|
   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: ISINTRIANGLE1" 

   ISINTRIANGLE1 = .FALSE.  

   IF(Y0 < MINVAL(YT) .OR. Y0 > MAXVAL(YT))THEN
     ISINTRIANGLE1 = .FALSE.
     RETURN
   END IF
   IF(X0 < MINVAL(XT) .OR. X0 > MAXVAL(XT))THEN
     ISINTRIANGLE1 = .FALSE.
     RETURN
   END IF

   F1 = (Y0-YT(1))*(XT(2)-XT(1)) - (X0-XT(1))*(YT(2)-YT(1))
   F2 = (Y0-YT(3))*(XT(1)-XT(3)) - (X0-XT(3))*(YT(1)-YT(3))
   F3 = (Y0-YT(2))*(XT(3)-XT(2)) - (X0-XT(2))*(YT(3)-YT(2))
  IF(F1*F3 >= 0.0_SP .AND. F3*F2 >= 0.0_SP) ISINTRIANGLE1 = .TRUE.

   IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: ISINTRIANGLE1" 
   RETURN
   END FUNCTION ISINTRIANGLE1
# endif
!==============================================================================|
END MODULE MOD_ENKF_OBS
