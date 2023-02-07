MODULE MOD_STATION_TIMESERIES

  USE MOD_TIME
  USE NETCDF
  USE VARS_WAVE
  IMPLICIT NONE
  SAVE
  
  LOGICAL OUT_STATION_TIMESERIES_ON
  CHARACTER(LEN=80) STATION_FILE
  CHARACTER(LEN=80) LOCATION_TYPE
  LOGICAL OUT_ELEVATION
  LOGICAL OUT_VELOCITY_3D
  LOGICAL OUT_VELOCITY_2D
  LOGICAL OUT_WIND_VELOCITY
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
  LOGICAL OUT_ATM_PRESSURE
# endif  
  LOGICAL OUT_SALT_TEMP
# if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
  LOGICAL OUT_SIG_WAVE_HEIGHT
  LOGICAL OUT_REL_PEAK_PERIOD
  LOGICAL OUT_WAVE_DIRECTION
  LOGICAL OUT_ENERGY_SPECTRUM
  LOGICAL OUT_WAVE_PARTITION
# endif  
  CHARACTER(LEN=80) OUT_INTERVAL
  
  NAMELIST /NML_STATION_TIMESERIES/        &
       & OUT_STATION_TIMESERIES_ON,        &
       & STATION_FILE,                     &
       & LOCATION_TYPE,                    &
       & OUT_ELEVATION,                    &
       & OUT_VELOCITY_3D,                  &
       & OUT_VELOCITY_2D,                  &
       & OUT_WIND_VELOCITY,                &
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
       & OUT_ATM_PRESSURE,                 &
# endif  
       & OUT_SALT_TEMP,                    &
# if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
       & OUT_SIG_WAVE_HEIGHT,              &
       & OUT_REL_PEAK_PERIOD,              &
       & OUT_WAVE_DIRECTION,               &
       & OUT_ENERGY_SPECTRUM,              &
       & OUT_WAVE_PARTITION,               &
# endif  
       & OUT_INTERVAL
       
  INTEGER NSTA       
  CHARACTER(LEN=20), ALLOCATABLE :: NAME_STA(:)
  REAL(SP), ALLOCATABLE :: LAT_STA(:),LON_STA(:),H_STA(:)
  INTEGER, ALLOCATABLE  :: NODE_STA(:),ELEMENT_STA(:),IDUMMY(:)
  INTEGER, ALLOCATABLE  :: NTVEGL(:),NBVEGL(:,:)
  
  TYPE(TIME) :: INTERVAL_TIME_SERIES, TIME_SERIES
  INTEGER    :: KD_START
  TYPE(TIME) :: KDD,KDD1
  
!--Control Variables----------------------------------------------!
   integer,private :: out_cnt            !!counts number of outputs
   integer,private :: stck_cnt           !!counts number of outputs in each file
   character(len=120),private :: cdfname !!netcdf file name
!--NetCDF IDs----------------------------------------------------!

   !--NetCDF File 
   integer,private :: nc_ofid

   !--Dimensions
   integer,private :: station_did,clen_did
   integer,private :: siglay_did,siglev_did
   integer,private :: time_did

   !--Grid Variables
   integer,private :: x_s_vid,y_s_vid,lat_s_vid,lon_s_vid
   integer,private :: siglay_vid,siglev_vid

   !--Flow Variables 
   integer,private :: time_s_vid
   integer,private :: iint_vid
   integer,private :: u_s_vid
   integer,private :: v_s_vid
   integer,private :: ww_s_vid
   integer,private :: s1_s_vid
   integer,private :: t1_s_vid
   integer,private :: el_s_vid
   integer,private :: h_s_vid
   integer,private :: ua_s_vid
   integer,private :: va_s_vid
   integer,private :: uuwind_s_vid
   integer,private :: vvwind_s_vid
   integer,private :: atmpres_s_vid
   integer,private :: name_s_vid
   
#  if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
   integer,private :: hs_s_vid
   integer,private :: pperiod_s_vid
   integer,private :: wdir_s_vid
   integer,private :: sdensity_s_vid
   integer,private :: msc_did
   integer,private :: hwind_s_vid
   integer,private :: dwind_s_vid
   integer,private :: twind_s_vid
   integer,private :: twindp_s_vid
   integer,private :: hswell_s_vid
   integer,private :: dswell_s_vid
   integer,private :: tswell_s_vid
   integer,private :: tswellp_s_vid
   integer,private :: imo_did   
#  endif   

   !--Info Variables
   character(len=120),public :: netcdf_timestring 


   INTERFACE PUTVAR
      MODULE PROCEDURE PUTVAR1D_INT
      MODULE PROCEDURE PUTVAR1D_REAL
      MODULE PROCEDURE PUTVAR2D_INT
      MODULE PROCEDURE PUTVAR2D_REAL
   END INTERFACE

  CONTAINS
!====================================================================================
   SUBROUTINE STATION_NAME_LIST_INITIALIZE
   
   IMPLICIT NONE
   
   OUT_STATION_TIMESERIES_ON = .False.
   STATION_FILE         = "'none'"
   LOCATION_TYPE             = "'node' or 'cell'"
   OUT_ELEVATION             = .False.
   OUT_VELOCITY_3D           = .False.
   OUT_VELOCITY_2D           = .False.
   OUT_WIND_VELOCITY         = .False.
#  if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
   OUT_ATM_PRESSURE          = .False.
#  endif   
   OUT_SALT_TEMP             = .False.
#  if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
   OUT_SIG_WAVE_HEIGHT       = .False.
   OUT_REL_PEAK_PERIOD       = .False.
   OUT_WAVE_DIRECTION        = .False.
   OUT_ENERGY_SPECTRUM       = .False.
   OUT_WAVE_PARTITION        = .False. 
#  endif  
   OUT_INTERVAL              = "A length of time: 'seconds= ','days= ', or 'cycles= '"
   
   
   RETURN
   END SUBROUTINE STATION_NAME_LIST_INITIALIZE
!----------------------------------------------------------------------------------
   SUBROUTINE STATION_NAME_LIST_PRINT
   USE CONTROL, ONLY : IPT
   
   IMPLICIT NONE
   
   WRITE(UNIT=IPT,NML=NML_STATION_TIMESERIES)
   
   RETURN
   END SUBROUTINE STATION_NAME_LIST_PRINT
!----------------------------------------------------------------------------------   
   SUBROUTINE STATION_NAME_LIST_READ
   USE CONTROL, ONLY : casename,NMLUNIT
   USE MOD_UTILS
   USE MOD_SET_TIME, ONLY : GET_OUTPUT_FILE_INTERVAL
   
   IMPLICIT NONE
   
   INTEGER :: IOS, I
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=160) :: PATHNFILE
   
   IF(DBG_SET(DBG_SBR)) &
         & WRITE(IPT,*) "Subroutine Begins: Read_Station_Name_List;"

   IOS = 0

   FNAME = "./"//trim(casename)//"_run.nml"

   CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

   !READ NAME LIST FILE
    REWIND(NMLUNIT)

   ! Read IO Information
   READ(UNIT=NMLUNIT, NML=NML_STATION_TIMESERIES,IOSTAT=ios)
   if(ios .NE. 0 ) Then
     if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_STATION_TIMESERIES)
     Call Fatal_error("Can Not Read NameList NML_STATION_TIMESERIES from file: "//trim(FNAME))
   end if
   CLOSE(NMLUNIT)

!   CALL GET_OUTPUT_FILE_INTERVAL(TRIM(OUT_INTERVAL),INTERVAL_TIME_SERIES)

   RETURN
   END SUBROUTINE STATION_NAME_LIST_READ
!----------------------------------------------------------------------------------   
   SUBROUTINE READ_STATION_FILE
   USE CONTROL, ONLY : MSR,input_dir,output_dir,casename,USE_REAL_WORLD_TIME,DATE_REFERENCE
   USE LIMS, ONLY : MGL,NGL,MYID,MSRID
   USE ALL_VARS, ONLY : SERIAL, PAR, NPROCS, H, ONE_THIRD, NVG,START_DATE
   USE MOD_UTILS, ONLY : Fatal_error
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR, ONLY : NMAP,EMAP,ACOLLECT
#  endif
   USE MOD_UTILS, ONLY : FOPEN
   USE MOD_OBCS, ONLY : GDAY1
   
   IMPLICIT NONE
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=3)   :: NAC
   INTEGER            :: IZAJ_MAX,IOS,IDUMMY1,I
   LOGICAL            :: FEXIST
   
   REAL(SP), ALLOCATABLE :: FTEMP(:),FTEMPC(:)
   INTEGER               :: IDD,IMM,IYY,ICC,IHH,IMI,ISS 
   INTEGER               :: IDD1,IMM1,IYY1,ICC1,IHH1,IMI1,ISS1 
   INTEGER               :: KD_REFERENCE  
   
   out_cnt = 0

   ALLOCATE(FTEMP(MGL))  ; FTEMP = 0.0_SP
   IF(SERIAL) FTEMP = H
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,H, FTEMP)
   ENDIF
#  endif

   IF(TRIM(LOCATION_TYPE) == "cell" .OR. TRIM(LOCATION_TYPE) == "CELL")THEN
     ALLOCATE(FTEMPC(NGL)) ;FTEMPC = 0.0_SP
     DO I = 1,NGL
       FTEMPC(I) = ONE_THIRD*(FTEMP(NVG(I,1))+FTEMP(NVG(I,2))+FTEMP(NVG(I,3)))
     END DO
   END IF    
     
   FNAME = trim(input_dir)//"/"//trim(STATION_FILE)
   CALL FOPEN(127,trim(FNAME),'cfr')

   IF(MSR)THEN
     IZAJ_MAX = 0
     READ(127,*,IOSTAT=IOS)
     DO WHILE(.TRUE.)
       READ(127,*,IOSTAT=IOS)idummy1
       IF(IOS < 0)EXIT
       IZAJ_MAX = IZAJ_MAX + 1
     END DO
     NSTA=IZAJ_MAX
     PRINT *,'NSTA=',NSTA
     ALLOCATE(NAME_STA(NSTA))
     ALLOCATE(LAT_STA(NSTA))
     ALLOCATE(LON_STA(NSTA))
     ALLOCATE(H_STA(NSTA))
     ALLOCATE(IDUMMY(NSTA))
     IF(TRIM(LOCATION_TYPE) == 'node' .OR. TRIM(LOCATION_TYPE) == 'NODE')THEN
       ALLOCATE(NODE_STA(NSTA))
     ELSE IF(TRIM(LOCATION_TYPE) == 'cell' .OR. TRIM(LOCATION_TYPE) == 'CELL')THEN  
       ALLOCATE(ELEMENT_STA(NSTA))
     ELSE
       CALL Fatal_error("LOCATION_TYPE should be either node or cell")
     END IF  
     REWIND(127)
 
     READ(127,*)
     DO I=1,NSTA
       IF(TRIM(LOCATION_TYPE) == "node" .OR. TRIM(LOCATION_TYPE) == "NODE")THEN
         READ(127,*)IDUMMY(I),LON_STA(I),LAT_STA(I),NODE_STA(I),H_STA(I),NAME_STA(I)
         WRITE(6,*) IDUMMY(I),LON_STA(I),LAT_STA(I),NODE_STA(I),H_STA(I),FTEMP(NODE_STA(I)),NAME_STA(I)
       ELSE IF(TRIM(LOCATION_TYPE) == "cell" .OR. TRIM(LOCATION_TYPE) == "CELL")THEN
         READ(127,*)IDUMMY(I),LON_STA(I),LAT_STA(I),ELEMENT_STA(I),H_STA(I),NAME_STA(I)
         WRITE(6,*) IDUMMY(I),LON_STA(I),LAT_STA(I),ELEMENT_STA(I),H_STA(I),FTEMPC(ELEMENT_STA(I)),NAME_STA(I)
       END IF	 
     ENDDO
     CLOSE(127)
     
   END IF

   DEALLOCATE(FTEMP)
   IF(TRIM(LOCATION_TYPE) == "cell" .OR. TRIM(LOCATION_TYPE) == "CELL")DEALLOCATE(FTEMPC)
   
   IF(USE_REAL_WORLD_TIME)THEN
     READ(START_DATE(1:2),*) ICC
     READ(START_DATE(3:4),*) IYY
     READ(START_DATE(6:7),*) IMM
     READ(START_DATE(9:10),*) IDD
     READ(START_DATE(12:13),*) IHH
     READ(START_DATE(15:16),*) IMI
     READ(START_DATE(18:19),*) ISS
     
     CALL GDAY1(IDD,IMM,IYY,ICC,KD_START)

     IF(DATE_REFERENCE /= 'default')THEN
       READ(DATE_REFERENCE(1:2),*) ICC1
       READ(DATE_REFERENCE(3:4),*) IYY1
       READ(DATE_REFERENCE(6:7),*) IMM1
       READ(DATE_REFERENCE(9:10),*) IDD1
       READ(DATE_REFERENCE(12:13),*) IHH1
       READ(DATE_REFERENCE(15:16),*) IMI1
       READ(DATE_REFERENCE(18:19),*) ISS1

       CALL GDAY1(IDD1,IMM1,IYY1,ICC1,KD_REFERENCE)
     ELSE
       CALL GDAY1(17,11,58,18,KD_REFERENCE)
     END IF  

     KD_START = KD_START - KD_REFERENCE
     KDD%MJD = KD_START      
     KDD%MuSOD = IHH*3600.+IMI*60.+ISS
   END IF
   
   RETURN
   END SUBROUTINE READ_STATION_FILE
!----------------------------------------------------------------------------------
   SUBROUTINE OUT_STATION_TIMESERIES 
   
   USE ALL_VARS
   USE NETCDF
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
   USE SWCOMM3, ONLY : MSC
   USE VARS_WAVE, ONLY : HSC1,TPEAK,DIRDEG1,SPEC_DENSITY
#  endif
   use mod_nctools,only : handle_ncerr
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP,T1TMP,S1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: UATMP,VATMP,UUWINDTMP,VVWINDTMP,ELTMP 
   REAL(SP) :: XTT,YTT
   INTEGER :: I1,I2,J 

   INTEGER  ::IZAJ_MAX,IZAJ_MIN,IZAJ,IZAJ_MAX_s,IZAJ_MIN_s,IZAJ_MAX_t,IZAJ_MIN_t
   INTEGER  ::kZAJ_MAX_s,kZAJ_MIN_s,kZAJ_MAX_t,kZAJ_MIN_t
   INTEGER  ::k,ierr
   REAL(SP) :: THOUR,THOUR1
   REAL(SP) :: STATMP(NSTA),STATMP1(NSTA),STATMPT(NSTA,KB),STATMPS(NSTA,KB)
   integer :: dims(1)
   real(sp), allocatable :: ftemp(:)
   integer :: VARID
   REAL(SP) :: KDD_TMP

!------------------------------------------------------------------------------!
!  WRITE TO FILES (SERIAL EXECUTION)                                           !
!------------------------------------------------------------------------------!
   IF(TIME_SERIES > IntTime) RETURN
   
   TIME_SERIES = IntTime + INTERVAL_TIME_SERIES
   THOUR = DTI*FLOAT(IINT-ISTART+1)/3600.0_SP
   THOUR1 = DTI*FLOAT(IINT)/3600.0_SP
!  WRITE ELEVATION, SALINITY, TEMPERATURE, VELOCITY AND WIND DATA
   out_cnt = out_cnt + 1
   stck_cnt = stck_cnt + 1 
   if(out_cnt == 1) call write_netcdf_setup

   dims(1) = stck_cnt
   
!--Open File
   if(msr)then
     ierr = nf90_open(cdfname,nf90_write,nc_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"file open error")
     end if

!--Dump Time/IINT to File
     ierr    = nf90_put_var(nc_ofid,iint_vid,iint,START=dims)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"error writing variable to netcdf")
     end if

     IF(USE_REAL_WORLD_TIME)THEN
       KDD1%MJD = KDD%MJD + INT((KDD%MuSOD/3600.+THOUR)/24.0)
       KDD1%MuSOD = KDD%MuSOD + THOUR * 3600 - INT((KDD%MuSOD/3600.+THOUR)/24.0) * 24 * 3600
       KDD_TMP = KDD1%MJD + KDD1%MuSOD/86400.0

       ierr    = nf90_put_var(nc_ofid,time_s_vid,kdd_tmp,START=dims)
       if(ierr /= nf90_noerr)then
         call handle_ncerr(ierr,"error writing variable to netcdf")
       end if

     ELSE
       ierr    = nf90_put_var(nc_ofid,time_s_vid,thour1*3600.,START=dims)
       if(ierr /= nf90_noerr)then
         call handle_ncerr(ierr,"error writing variable to netcdf")
       end if
     END IF
   end if

!--Write Variables to File
   if(msr) write(ipt,*)'dumping to netcdf file: ',trim(cdfname),stck_cnt

     IF(OUT_ELEVATION)THEN
       i1 = lbound(el,1) ; i2 = ubound(el,1)
       call putvar(i1,i2,m,mgl,1,1,"n",el,nc_ofid,el_s_vid,myid,nprocs&
            &,ipt, stck_cnt) 
     END IF

     IF(OUT_SALT_TEMP)THEN 
       i1 = lbound(t1,1) ; i2 = ubound(t1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",t1,nc_ofid,t1_s_vid,myid&
            &,nprocs,ipt, stck_cnt) 
       i1 = lbound(s1,1) ; i2 = ubound(s1,1)
       call putvar(i1,i2,m,mgl,kb,kb-1,"n",s1,nc_ofid,s1_s_vid,myid&
            &,nprocs,ipt, stck_cnt) 
     END IF

     IF(OUT_VELOCITY_3D)THEN 
       i1 = lbound(u,1) ; i2 = ubound(u,1)
       call putvar(i1,i2,n,ngl,kb,kb-1,"e",u,nc_ofid,u_s_vid,myid&
            &,nprocs,ipt, stck_cnt) 
       i1 = lbound(v,1) ; i2 = ubound(v,1)
       call putvar(i1,i2,n,ngl,kb,kb-1,"e",v,nc_ofid,v_s_vid,myid&
            &,nprocs,ipt, stck_cnt) 
       i1 = lbound(ww,1) ; i2 = ubound(ww,1)
       call putvar(i1,i2,n,ngl,kb,kb-1,"e",ww,nc_ofid,ww_s_vid,myid&
            &,nprocs,ipt, stck_cnt) 
     END IF 

     IF(OUT_VELOCITY_2D)THEN 
       allocate(ftemp(n))
       ftemp =ua(1:n)
       i1 = lbound(ftemp,1) ; i2 = ubound(ftemp,1)
       call putvar(i1,i2,n,ngl,1,1,"e",ftemp,nc_ofid,ua_s_vid,myid&
            &,nprocs,ipt, stck_cnt)
       deallocate(ftemp)
       allocate(ftemp(n))
       ftemp =va(1:n)
       i1 = lbound(ftemp,1) ; i2 = ubound(ftemp,1)
       call putvar(i1,i2,n,ngl,1,1,"e",ftemp,nc_ofid,va_s_vid,myid&
            &,nprocs,ipt, stck_cnt)
       deallocate(ftemp)
     END IF 

     IF(OUT_WIND_VELOCITY)THEN 
       allocate(ftemp(n))
       ftemp =uuwind(1:n)
       i1 = lbound(ftemp,1) ; i2 = ubound(ftemp,1)
       call putvar(i1,i2,n,ngl,1,1,"e",ftemp,nc_ofid,uuwind_s_vid,myid&
            &,nprocs,ipt, stck_cnt)
       deallocate(ftemp)
       allocate(ftemp(n))
       ftemp =vvwind(1:n)
       i1 = lbound(ftemp,1) ; i2 = ubound(ftemp,1)
       call putvar(i1,i2,n,ngl,1,1,"e",ftemp,nc_ofid,vvwind_s_vid,myid&
            &,nprocs,ipt, stck_cnt)
       deallocate(ftemp)
     END IF 

#    if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
     IF(OUT_ATM_PRESSURE)THEN 
       allocate(ftemp(m))
       ftemp =pa_air(1:m)
       i1 = lbound(ftemp,1) ; i2 = ubound(ftemp,1)
       call putvar(i1,i2,m,mgl,1,1,"n",ftemp,nc_ofid,atmpres_s_vid,myid&
            &,nprocs,ipt, stck_cnt)
       deallocate(ftemp)
     END IF 
#    endif

#    if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
     IF(OUT_SIG_WAVE_HEIGHT)THEN
       i1 = lbound(hsc1,1) ; i2 = ubound(hsc1,1)
       call putvar(i1,i2,m,mgl,1,1,"n",hsc1,nc_ofid,hs_s_vid,myid,nprocs&
            &,ipt, stck_cnt) 
     END IF

     IF(OUT_REL_PEAK_PERIOD)THEN
       i1 = lbound(tpeak,1) ; i2 = ubound(tpeak,1)
       call putvar(i1,i2,m,mgl,1,1,"n",tpeak,nc_ofid,pperiod_s_vid,myid,nprocs&
            &,ipt, stck_cnt) 
     END IF

     IF(OUT_WAVE_DIRECTION)THEN
       i1 = lbound(dirdeg1,1) ; i2 = ubound(dirdeg1,1)
       call putvar(i1,i2,m,mgl,1,1,"n",dirdeg1,nc_ofid,wdir_s_vid,myid,nprocs&
            &,ipt, stck_cnt) 
     END IF

     IF(OUT_ENERGY_SPECTRUM)THEN
       i1 = lbound(spec_density,1) ; i2 = ubound(spec_density,1)
       call putvar(i1,i2,m,mgl,msc,msc,"n",spec_density,nc_ofid,sdensity_s_vid,myid,nprocs&
            &,ipt, stck_cnt) 
     END IF

     IF(OUT_WAVE_PARTITION)THEN
     
       IF(TRIM(LOCATION_TYPE) == 'cell')THEN
         CALL Fatal_error("LOCATION_TYPE should be node for wave partition output") 
       END IF 

       i1 = lbound(hs_wind,1) ; i2 = ubound(hs_wind,1)
       call putvar(i1,i2,m,mgl,1,1,"n",hs_wind,nc_ofid,hwind_s_vid,myid,nprocs&
            &,ipt, stck_cnt)

       i1 = lbound(dirdeg_wind,1) ; i2 = ubound(dirdeg_wind,1)
       call putvar(i1,i2,m,mgl,1,1,"n",dirdeg_wind,nc_ofid,dwind_s_vid,myid,nprocs&
            &,ipt, stck_cnt)

       i1 = lbound(tpeak_wind,1) ; i2 = ubound(tpeak_wind,1)
       call putvar(i1,i2,m,mgl,1,1,"n",tpeak_wind,nc_ofid,twind_s_vid,myid,nprocs&
            &,ipt, stck_cnt)

       i1 = lbound(tpeak_wind_pos,1) ; i2 = ubound(tpeak_wind_pos,1)
       call putvar(i1,i2,m,mgl,1,1,"n",tpeak_wind_pos,nc_ofid,twindp_s_vid,myid,nprocs&
            &,ipt, stck_cnt)

       i1 = lbound(hs_swell_all,1) ; i2 = ubound(hs_swell_all,1)
       call putvar(i1,i2,m,mgl,50,50,"n",hs_swell_all,nc_ofid,hswell_s_vid,myid,nprocs&
            &,ipt, stck_cnt)

       i1 = lbound(dirdeg_swell_all,1) ; i2 = ubound(dirdeg_swell_all,1)
      call putvar(i1,i2,m,mgl,50,50,"n",dirdeg_swell_all,nc_ofid,dswell_s_vid,myid,nprocs&
            &,ipt, stck_cnt)

       i1 = lbound(tpeak_swell_all,1) ; i2 = ubound(tpeak_swell_all,1)
       call putvar(i1,i2,m,mgl,50,50,"n",tpeak_swell_all,nc_ofid,tswell_s_vid,myid,nprocs&
            &,ipt, stck_cnt)

       i1 = lbound(tpeak_swell_pos_all,1) ; i2 = ubound(tpeak_swell_pos_all,1)
       call putvar(i1,i2,m,mgl,50,50,"n",tpeak_swell_pos_all,nc_ofid,tswellp_s_vid,myid,nprocs&
            &,ipt, stck_cnt)
     END IF

#    endif

  IERR = NF90_CLOSE(NC_OFID)

   RETURN
   END SUBROUTINE out_station_timeseries
!==============================================================================|

!==============================================================================|
!  Write NetCDF Header and Static Variables                                    |
!==============================================================================|
   SUBROUTINE write_netcdf_setup 

   use all_vars
   use mod_clock
   use mod_nctools
#  if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
   use swcomm3, only : msc
#  endif

#  if defined (MULTIPROCESSOR)
   use mod_par 
#  endif
   use netcdf
   use mod_types
   use mod_utils
   implicit none
   integer, dimension(3) :: dynm2de_lay,dynm2dn_lay
   integer, dimension(2) :: dynm2ds
   integer, dimension(1) :: stat2ds
   integer, dimension(2) :: stat2ds_lev,stat2ds_lay 
   integer, dimension(2) :: stat2ds_char
   integer, dimension(1) :: dynmtime
   character(len=100)    :: netcdf_convention
   character(len=100)    :: timestamp ,temp
   integer               :: i,j,ierr,i1,i2
#  if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
   integer, dimension(3) :: dynm2dn_msc
   integer, dimension(3) :: dynm2dn_imo
#  endif

!==============================================================================|

!==============================================================================|
!  Set up Constants and Initialize Counters                                    |
!==============================================================================|
NETCDF_TIMESTRING = 'seconds after 00:00:00'
!--Initialize Stack Count
   stck_cnt = 1

!--NetCDF Convention String
   netcdf_convention = 'CF-1.0'

!--Time Stamp for History
   call get_timestamp(temp)
   timestamp = 'model started at: '//trim(temp)


!==============================================================================|
!  OPEN FILE AND DEFINE VARIABLES                                              |
!==============================================================================|
   if(msr)then

    cdfname = trim(OUTPUT_DIR)//trim(casename)//'_station_timeseries.nc'

!--Create File 
    ierr = nf90_create(path=cdfname,cmode=nf90_clobber,ncid=nc_ofid)
    if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"file creation error")
    end if

!--Description of File Contents
    ierr = nf90_put_att(nc_ofid,nf90_global,"title"      ,trim(case_title))
    ierr = nf90_put_att(nc_ofid,nf90_global,"institution",trim(institution))
    ierr = nf90_put_att(nc_ofid,nf90_global,"source"     ,trim(fvcom_version))
    ierr = nf90_put_att(nc_ofid,nf90_global,"history"    ,trim(timestamp))
    ierr = nf90_put_att(nc_ofid,nf90_global,"references" ,trim(fvcom_website))
    ierr = nf90_put_att(nc_ofid,nf90_global,"Conventions",trim(netcdf_convention))
#   if defined (SPHERICAL)
    ierr = nf90_put_att(nc_ofid,nf90_global,"CoordinateSystem","GeoReferenced")
#   endif

!--Define Fixed Model Dimensions 
    ierr = nf90_def_dim(nc_ofid,"siglay" ,kbm1   ,siglay_did )
    ierr = nf90_def_dim(nc_ofid,"siglev" ,kb     ,siglev_did )
    ierr = nf90_def_dim(nc_ofid,"station"   ,nsta      ,station_did   )
    ierr = nf90_def_dim(nc_ofid,"clen", 20, clen_did )
#   if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
    ierr = nf90_def_dim(nc_ofid,"msc",msc,msc_did)
    ierr = nf90_def_dim(nc_ofid,"IMO",50,imo_did)
#   endif    

!--Define Unlimited Model Dimension
    ierr = nf90_def_dim(nc_ofid,"time"   ,NF90_UNLIMITED,time_did)

!--Set Up Data Dimensioning - Static Vars
    stat2ds      = (/station_did/)            !!2d station vars
    stat2ds_lev  = (/station_did,siglev_did/)
    stat2ds_lay  = (/station_did,siglay_did/)
    stat2ds_char = (/clen_did,station_did/)

!--Set Up Data Dimensioning - Dynamic Vars 
    dynm2ds      = (/station_did,time_did/)            !!2d station vars
    dynm2de_lay  = (/station_did,siglay_did,time_did/) 
    dynm2dn_lay  = (/station_did,siglay_did,time_did/) 
    dynmtime     = (/time_did/)   
#   if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
    dynm2dn_msc  = (/station_did,msc_did,time_did/)
    dynm2dn_imo  = (/station_did,imo_did,time_did/)
#   endif    

!--Define Station Name Variables and Attributes

    !!====Station Name (NAME_STA)  ===================!
    ierr = nf90_def_var(nc_ofid,"name_station",nf90_char,stat2ds_char,name_s_vid)
    ierr = nf90_put_att(nc_ofid,name_s_vid,"long_name","Station Name")   

!--Define Coordinate Variables and Attributes
    !!====X Grid Coordinate at Nodes (VX) (Meters)===========!
    ierr = nf90_def_var(nc_ofid,"x",nf90_float,stat2ds,x_s_vid)
    ierr = nf90_put_att(nc_ofid,x_s_vid,"long_name","station x-coordinate")
    ierr = nf90_put_att(nc_ofid,x_s_vid,"units","meters")

    !!====Y Grid Coordinate at Nodes (VY) (Meters)===========!
    ierr = nf90_def_var(nc_ofid,"y",nf90_float,stat2ds,y_s_vid)
    ierr = nf90_put_att(nc_ofid,y_s_vid,"long_name","station y-coordinate")
    ierr = nf90_put_att(nc_ofid,y_s_vid,"units","meters")

    !!====Longitudinal Coordinate at Nodes (LON) (degrees)===!
    ierr = nf90_def_var(nc_ofid,"lon",nf90_float,stat2ds,lon_s_vid)
    ierr = nf90_put_att(nc_ofid,lon_s_vid,"long_name","Longitude")
    ierr = nf90_put_att(nc_ofid,lon_s_vid,"standard_name","longitude")
    ierr = nf90_put_att(nc_ofid,lon_s_vid,"units","degrees_east")

    !!====Latitudinal  Coordinate at Nodes (LAT) (degrees)===!
    ierr = nf90_def_var(nc_ofid,"lat",nf90_float,stat2ds,lat_s_vid)
    ierr = nf90_put_att(nc_ofid,lat_s_vid,"long_name","Latitude")
    ierr = nf90_put_att(nc_ofid,lat_s_vid,"standard_name","latitude")
    ierr = nf90_put_att(nc_ofid,lat_s_vid,"units","degrees_north")
    ierr = nf90_put_att(nc_ofid,lat_s_vid,"grid","Bathymetry_Mesh")

   !!====Sigma Coordinate for Sigma Layers (ZZ)  (-)========!
   ierr = nf90_def_var(nc_ofid,"siglay",nf90_float,stat2ds_lay,siglay_vid)
   ierr = nf90_put_att(nc_ofid,siglay_vid,"long_name","Sigma Layers")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"standard_name","ocean_sigma/general_coordinate")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"formula_terms","siglay:siglay eta:zeta depth:depth")

   !!====Sigma Coordinate for Sigma Levels (Z)   (-)========!
   ierr = nf90_def_var(nc_ofid,"siglev",nf90_float,stat2ds_lev,siglev_vid)
   ierr = nf90_put_att(nc_ofid,siglev_vid,"long_name","Sigma Levels")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"standard_name","ocean_sigma/general_coordinate")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"formula_terms","siglev:siglev eta:zeta depth:depth")

!--Define Mesh Relevant Variables and Attributes

    !!====Bathymetry at Nodes (H) (meters)===================!
    ierr = nf90_def_var(nc_ofid,"h",nf90_float,stat2ds,h_s_vid)
    ierr = nf90_put_att(nc_ofid,h_s_vid,"long_name","Bathymetry")   
    ierr = nf90_put_att(nc_ofid,h_s_vid,"units","meters")
    ierr = nf90_put_att(nc_ofid,h_s_vid,"positive","down")
    ierr = nf90_put_att(nc_ofid,h_s_vid,"standard_name","depth")

!--Define Model Time Variables and Attributes    
    IF(USE_REAL_WORLD_TIME)THEN
      ierr = nf90_def_var(nc_ofid,"time",nf90_float,dynmtime,time_s_vid)
      ierr = nf90_put_att(nc_ofid,time_s_vid,"long_name","time")
      if(DATE_REFERENCE == 'default')then
        ierr = nf90_put_att(nc_ofid,time_s_vid,"units",trim("days since 1858-11-17 00:00:00"))
        ierr = nf90_put_att(nc_ofid,time_s_vid,"format",trim("modified julian day (MJD)"))
      else
        ierr = nf90_put_att(nc_ofid,time_s_vid,"units","days since "//trim(DATE_REFERENCE))
        ierr = nf90_put_att(nc_ofid,time_s_vid,"format",trim("defined reference date"))
      end if
!JQI      ierr = nf90_put_att(nc_ofid,time_s_vid,"calendar","none")
      ierr = nf90_put_att(nc_ofid,time_s_vid,"time_zone","UTC")

    ELSE
      ierr = nf90_def_var(nc_ofid,"time",nf90_float,dynmtime,time_s_vid)
      ierr = nf90_put_att(nc_ofid,time_s_vid,"long_name","time")
      ierr = nf90_put_att(nc_ofid,time_s_vid,"units",trim(netcdf_timestring))
!JQI      ierr = nf90_put_att(nc_ofid,time_s_vid,"calendar","none")
      ierr = nf90_put_att(nc_ofid,time_s_vid,"time_zone","none")
    END IF  

    ierr = nf90_def_var(nc_ofid,"iint",nf90_int,dynmtime,iint_vid)
    ierr = nf90_put_att(nc_ofid,iint_vid,"long_name","internal mode iteration number")

!--Define Time Dependent Flow Variables (selected by user from input file)
    if(OUT_VELOCITY_3D)then
     ierr = nf90_def_var(nc_ofid,"u",nf90_float,dynm2de_lay,u_s_vid)
     ierr = nf90_put_att(nc_ofid,u_s_vid,"long_name","Eastward Water Velocity")
     ierr = nf90_put_att(nc_ofid,u_s_vid,"standard_name","eastward_sea_water_velocity")
     ierr = nf90_put_att(nc_ofid,u_s_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,u_s_vid,"type","data")
     ierr = nf90_put_att(nc_ofid,u_s_vid,"coordinates","time siglay station")
       
     ierr = nf90_def_var(nc_ofid,"v",nf90_float,dynm2de_lay,v_s_vid)
     ierr = nf90_put_att(nc_ofid,v_s_vid,"long_name","Northward Water Velocity")
     ierr = nf90_put_att(nc_ofid,u_s_vid,"standard_name","northward_sea_water_velocity")
     ierr = nf90_put_att(nc_ofid,v_s_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,v_s_vid,"type","data")
     ierr = nf90_put_att(nc_ofid,v_s_vid,"coordinates","time siglay station")

     ierr = nf90_def_var(nc_ofid,"ww",nf90_float,dynm2de_lay,ww_s_vid)
     ierr = nf90_put_att(nc_ofid,ww_s_vid,"long_name","Upward Water Velocity")
     ierr = nf90_put_att(nc_ofid,ww_s_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,ww_s_vid,"type","data")
    end if

    if(OUT_VELOCITY_2D)then
     ierr = nf90_def_var(nc_ofid,"ua",nf90_float,dynm2ds,ua_s_vid)
     ierr = nf90_put_att(nc_ofid,ua_s_vid,"long_name","Vertically Averaged x-velocity")
     ierr = nf90_put_att(nc_ofid,ua_s_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,ua_s_vid,"type","data")
       
     ierr = nf90_def_var(nc_ofid,"va",nf90_float,dynm2ds,va_s_vid)
     ierr = nf90_put_att(nc_ofid,va_s_vid,"long_name","Vertically Averaged y-velocity")
     ierr = nf90_put_att(nc_ofid,va_s_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,va_s_vid,"type","data")
    end if

    if(OUT_SALT_TEMP)then
     ierr = nf90_def_var(nc_ofid,"temp",nf90_float,dynm2dn_lay,t1_s_vid)
     ierr = nf90_put_att(nc_ofid,t1_s_vid,"long_name","temperature")
     ierr = nf90_put_att(nc_ofid,t1_s_vid,"standard_name","sea_water_temperature")
     ierr = nf90_put_att(nc_ofid,t1_s_vid,"units","degrees_C")
     ierr = nf90_put_att(nc_ofid,t1_s_vid,"type","data")
     ierr = nf90_put_att(nc_ofid,t1_s_vid,"coordinates","time siglay station")

     ierr = nf90_def_var(nc_ofid,"salinity",nf90_float,dynm2dn_lay,s1_s_vid)
     ierr = nf90_put_att(nc_ofid,s1_s_vid,"long_name","salinity")
     ierr = nf90_put_att(nc_ofid,s1_s_vid,"standard_name","sea_water_salinity")
     ierr = nf90_put_att(nc_ofid,s1_s_vid,"units","1e-3")
     ierr = nf90_put_att(nc_ofid,s1_s_vid,"type","data")
     ierr = nf90_put_att(nc_ofid,s1_s_vid,"coordinates","time siglay station")
    end if

    if(OUT_ELEVATION)then
     ierr = nf90_def_var(nc_ofid,"zeta",nf90_float,dynm2ds,el_s_vid)
     ierr = nf90_put_att(nc_ofid,el_s_vid,"long_name","Water Surface Elevation")
     ierr = nf90_put_att(nc_ofid,el_s_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,el_s_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,el_s_vid,"standard_name","sea_surface_height_above_geoid")
     ierr = nf90_put_att(nc_ofid,el_s_vid,"type","data")
     ierr = nf90_put_att(nc_ofid,el_s_vid,"coordinates","time station")
    end if

    if(OUT_WIND_VELOCITY)then
     ierr = nf90_def_var(nc_ofid,"uwind_speed",nf90_float,dynm2ds,uuwind_s_vid)
     ierr = nf90_put_att(nc_ofid,uuwind_s_vid,"long_name","Eastward wind velocity")
     ierr = nf90_put_att(nc_ofid,uuwind_s_vid,"units","(m/s)")
     ierr = nf90_put_att(nc_ofid,uuwind_s_vid,"standard_name","eastward wind")
     ierr = nf90_put_att(nc_ofid,uuwind_s_vid,"type","data")
     ierr = nf90_put_att(nc_ofid,uuwind_s_vid,"coordinates","time station")

     ierr = nf90_def_var(nc_ofid,"vwind_speed",nf90_float,dynm2ds,vvwind_s_vid)
     ierr = nf90_put_att(nc_ofid,vvwind_s_vid,"long_name","Northward wind velocity")
     ierr = nf90_put_att(nc_ofid,vvwind_s_vid,"units","(m/s)")
     ierr = nf90_put_att(nc_ofid,vvwind_s_vid,"standard_name","northward wind")
     ierr = nf90_put_att(nc_ofid,vvwind_s_vid,"type","data")
     ierr = nf90_put_att(nc_ofid,vvwind_s_vid,"coordinates","time station")
    end if

#   if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
    if(OUT_ATM_PRESSURE)then
     ierr = nf90_def_var(nc_ofid,"atm_pressure",nf90_float,dynm2ds,atmpres_s_vid)
     ierr = nf90_put_att(nc_ofid,atmpres_s_vid,"long_name","Sea level pressure")
     ierr = nf90_put_att(nc_ofid,atmpres_s_vid,"units","(Pa)")
     ierr = nf90_put_att(nc_ofid,atmpres_s_vid,"standard_name","atmospheric pressure")
     ierr = nf90_put_att(nc_ofid,atmpres_s_vid,"type","data")
    end if
#   endif    

#   if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
    if(OUT_SIG_WAVE_HEIGHT)then
     ierr = nf90_def_var(nc_ofid,"wave_height",nf90_float,dynm2ds,hs_s_vid)
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"long_name","Signaficant Wave Height")
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"standard_name","signaficant_wave_height")
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"type","data")
    end if

    if(OUT_REL_PEAK_PERIOD)then
     ierr = nf90_def_var(nc_ofid,"wave_period",nf90_float,dynm2ds,pperiod_s_vid)
     ierr = nf90_put_att(nc_ofid,pperiod_s_vid,"long_name","Relative Peak Period")
     ierr = nf90_put_att(nc_ofid,pperiod_s_vid,"units","seconds")
     ierr = nf90_put_att(nc_ofid,pperiod_s_vid,"standard_name","relative_peak_period")
     ierr = nf90_put_att(nc_ofid,pperiod_s_vid,"type","data")
    end if

    if(OUT_WAVE_DIRECTION)then
     ierr = nf90_def_var(nc_ofid,"wave_dir",nf90_float,dynm2ds,wdir_s_vid)
     ierr = nf90_put_att(nc_ofid,wdir_s_vid,"long_name","Wave Direction")
     ierr = nf90_put_att(nc_ofid,wdir_s_vid,"units","degrees")
     ierr = nf90_put_att(nc_ofid,wdir_s_vid,"standard_name","wave_direction")
     ierr = nf90_put_att(nc_ofid,wdir_s_vid,"type","data")
    end if

    if(OUT_ENERGY_SPECTRUM)then
     ierr = nf90_def_var(nc_ofid,"energy_spectrum",nf90_float,dynm2dn_msc,sdensity_s_vid)
     ierr = nf90_put_att(nc_ofid,sdensity_s_vid,"long_name","Energy Spectrum")
     ierr = nf90_put_att(nc_ofid,sdensity_s_vid,"units","m^2/hz")
     ierr = nf90_put_att(nc_ofid,sdensity_s_vid,"standard_name","energy_dpectrum")
     ierr = nf90_put_att(nc_ofid,sdensity_s_vid,"type","data")
    end if

    if(OUT_WAVE_PARTITION)then
     ierr = nf90_def_var(nc_ofid,"hs_wind",nf90_float,dynm2ds,hwind_s_vid)
     ierr = nf90_put_att(nc_ofid,hwind_s_vid,"long_name","Wind Wave HS")
     ierr = nf90_put_att(nc_ofid,hwind_s_vid,"units","m")
     ierr = nf90_put_att(nc_ofid,hwind_s_vid,"standard_name","hs_windwave")
     ierr = nf90_put_att(nc_ofid,hwind_s_vid,"type","data") 

     ierr = nf90_def_var(nc_ofid,"dirdeg_wind",nf90_float,dynm2ds,dwind_s_vid)
     ierr = nf90_put_att(nc_ofid,dwind_s_vid,"long_name","Wind Wave DIR")
     ierr = nf90_put_att(nc_ofid,dwind_s_vid,"units","degree")
     ierr = nf90_put_att(nc_ofid,dwind_s_vid,"standard_name","dirdeg_windwave")
     ierr = nf90_put_att(nc_ofid,dwind_s_vid,"type","data")

     ierr = nf90_def_var(nc_ofid,"tpeak_wind",nf90_float,dynm2ds,twind_s_vid)
     ierr = nf90_put_att(nc_ofid,twind_s_vid,"long_name","Wind Wave HS")
     ierr = nf90_put_att(nc_ofid,twind_s_vid,"units","s")
     ierr = nf90_put_att(nc_ofid,twind_s_vid,"standard_name","tpeak_windwave")
     ierr = nf90_put_att(nc_ofid,twind_s_vid,"type","data")

     ierr = nf90_def_var(nc_ofid,"dirdeg_pos_wind",nf90_float,dynm2ds,twindp_s_vid)
     ierr = nf90_put_att(nc_ofid,twindp_s_vid,"long_name","Wind Wave DIR POS")
     ierr = nf90_put_att(nc_ofid,twindp_s_vid,"units","non")
     ierr = nf90_put_att(nc_ofid,twindp_s_vid,"standard_name","dirdeg_pos_windwave")
     ierr = nf90_put_att(nc_ofid,twindp_s_vid,"type","data")

     ierr = nf90_def_var(nc_ofid,"hs_swell",nf90_float,dynm2dn_imo,hswell_s_vid)
     ierr = nf90_put_att(nc_ofid,hswell_s_vid,"long_name","Swell HS")
     ierr = nf90_put_att(nc_ofid,hswell_s_vid,"units","m")
     ierr = nf90_put_att(nc_ofid,hswell_s_vid,"standard_name","hs_swell")
     ierr = nf90_put_att(nc_ofid,hswell_s_vid,"type","data")

     ierr = nf90_def_var(nc_ofid,"dirdeg_swell",nf90_float,dynm2dn_imo,dswell_s_vid)
     ierr = nf90_put_att(nc_ofid,dswell_s_vid,"long_name","Swell DIR")
     ierr = nf90_put_att(nc_ofid,dswell_s_vid,"units","degree")
     ierr = nf90_put_att(nc_ofid,dswell_s_vid,"standard_name","dirdeg_swell")
     ierr = nf90_put_att(nc_ofid,dswell_s_vid,"type","data")

     ierr = nf90_def_var(nc_ofid,"tpeak_swell",nf90_float,dynm2dn_imo,tswell_s_vid)
     ierr = nf90_put_att(nc_ofid,tswell_s_vid,"long_name","Swell HS")
     ierr = nf90_put_att(nc_ofid,tswell_s_vid,"units","s")
     ierr = nf90_put_att(nc_ofid,tswell_s_vid,"standard_name","tpeak_swell")
     ierr = nf90_put_att(nc_ofid,tswell_s_vid,"type","data")

     ierr = nf90_def_var(nc_ofid,"dirdeg_pos_swell",nf90_float,dynm2dn_imo,tswellp_s_vid)
     ierr = nf90_put_att(nc_ofid,tswellp_s_vid,"long_name","Swell DIR POS")
     ierr = nf90_put_att(nc_ofid,tswellp_s_vid,"units","non")
     ierr = nf90_put_att(nc_ofid,tswellp_s_vid,"standard_name","dirdeg_pos_swell")
     ierr = nf90_put_att(nc_ofid,tswellp_s_vid,"type","data")
    end if

#   endif
!--Exit Define Mode
    ierr = nf90_enddef(nc_ofid)
    ierr = nf90_close(nc_ofid)

   end if !(msr)

!==============================================================================|
!  WRITE VARIABLES TO FILE                                                     |
!==============================================================================|
   if(msr)then
     ierr = nf90_open(cdfname,nf90_write,nc_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"file open error")
     end if
   end if
   
   !!====Longitude at Nodes (LON) ==========================!
   i1 = lbound(lon,1) ; i2 = ubound(lon,1)
   call putvar(i1,i2,m,mgl,1,1,"n",lon,nc_ofid,lon_s_vid,myid&
        &,nprocs,ipt, stck_cnt)

   !!====Latitude  at Nodes (LAT) ==========================!
   i1 = lbound(lat,1) ; i2 = ubound(lat,1)
   call putvar(i1,i2,m,mgl,1,1,"n",lat,nc_ofid,lat_s_vid,myid&
        &,nprocs,ipt, stck_cnt) 

   !!====X Grid Coordinate at Nodes (VX)====================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar(i1,i2,m,mgl,1,1,"n",vx+vxmin,nc_ofid,x_s_vid,myid,nprocs&
        &,ipt, stck_cnt) 

   !!====Y Grid Coordinate at Nodes (VY)====================!
   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar(i1,i2,m,mgl,1,1,"n",vy+vymin,nc_ofid,y_s_vid,myid,nprocs&
        &,ipt, stck_cnt) 

   !!====Bathymetry at Nodes (H)============================!
   i1 = lbound(h,1) ; i2 = ubound(h,1)
   call putvar(i1,i2,m,mgl,1,1,"n",h,nc_ofid,h_s_vid,myid,nprocs,ipt,&
        & stck_cnt) 

   !!====Sigma Layers (zz)==================================!
   i1 = lbound(zz,1) ; i2 = ubound(zz,1)
   call putvar(i1,i2,m,mgl,kb,kb-1,"n",zz,nc_ofid,siglay_vid,myid&
        &,nprocs,ipt, stck_cnt) 

   !!====Sigma Levels (z)==================================!
   i1 = lbound(z,1) ; i2 = ubound(z,1)
   call putvar(i1,i2,m,mgl,kb,kb,"n",z,nc_ofid,siglev_vid,myid,nprocs&
        &,ipt, stck_cnt) 

   !!====Station Name (NAME_STA)============================!
   if(msr)then
    i1 = lbound(name_sta,1) ; i2 = ubound(name_sta,1)
    call putvar_char(i1,i2,name_sta,nc_ofid,name_s_vid,myid,nprocs,ipt,&
        & stck_cnt) 
   end if	

!==============================================================================|
!  close the file                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(nc_ofid)

   return
   end subroutine write_netcdf_setup
!==============================================================================|
 !==============================================================================|
 !  Collect Data to Global Array and Write to Netcdf File                       |
 !==============================================================================|
 
 SUBROUTINE PUTVAR1D_REAL(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid&
      &,myid,nprocs,ipt,stk)
   
   !------------------------------------------------------------------------------|
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs&
        &,ipt,stk
   character(len=*),intent(in)   :: map_type
   real(sp), dimension(i1:i2) :: var
   
   real(sp), allocatable, dimension(:,:) :: temp

   allocate(temp(i1:i2,kt))
   temp(i1:i2,1)=var

   CALL PUTVAR2D_REAL(i1,i2,n1,n1gl,kt,k1,map_type,temp,nc_fid,vid&
        &,myid,nprocs,ipt,stk)
   
   deallocate(temp)

 END SUBROUTINE PUTVAR1D_REAL
   

 subroutine PUTVAR2D_REAL(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid&
      &,myid,nprocs,ipt,stk)
!------------------------------------------------------------------------------|

#  if defined (MULTIPROCESSOR)
   use mod_par,only : nmap,emap,acollect
   use lims, only : msrid
#  endif
   use all_vars, only : nvg
   use mod_nctools,only : handle_ncerr
   use mod_types
   use mod_utils, only : pstop
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs&
        &,ipt, stk
   character(len=*),intent(in)   :: map_type
   real(sp), allocatable :: var(:,:)

   real(sp), allocatable, dimension(:,:) :: temp,gtemp
   integer :: ierr,k1m1,i,j
   integer, allocatable :: dims(:)
   

   k1m1 = k1 
   if(k1m1 == 1)then
     allocate(dims(2))
     dims(1) = 1 
     dims(2) = stk
   else
     allocate(dims(3))
     dims(1) = 1 
     dims(2) = 1 
     dims(3) = stk      
   end if
     

   if(map_type(1:1) /= "e" .and. map_type(1:1) /= "n")then
     write(ipt,*)'map_type input to putvar should be "e" OR "n"'
     call pstop
   end if

   if(nprocs==1)then
     allocate(temp(nsta,k1m1))  

     if(map_type(1:1) == 'n' .and. TRIM(LOCATION_TYPE) == 'node')then
      do i=1,nsta
       temp(i,1:k1m1) = var(node_sta(i),1:k1m1)
      end do  
     else if(map_type(1:1) == 'e' .and. TRIM(LOCATION_TYPE) == 'cell')then
      do i=1,nsta
       temp(i,1:k1m1) = var(element_sta(i),1:k1m1)
      end do  
     else if(map_type(1:1) == 'n' .and. TRIM(LOCATION_TYPE) == 'cell')then
      do i=1,nsta
       temp(i,:) = 0.0_SP
       do j =1,3
        temp(i,1:k1m1) = temp(i,1:k1m1) + var(nvg(element_sta(i),j),1:k1m1)
       end do
       temp(i,:) = temp(i,:)/3.0_SP
      end do   
     else if(map_type(1:1) == 'e' .and. TRIM(LOCATION_TYPE) == 'node')then
      do i=1,nsta
       temp(i,:) = 0.0_SP
       do j =1, ntvegl(node_sta(i))
        temp(i,1:k1m1) = temp(i,1:k1m1) + var(nbvegl(node_sta(i),j),1:k1m1)
       end do
       temp(i,:) = temp(i,:)/ntvegl(node_sta(i))
      end do   
     end if
   end if

#  if defined (MULTIPROCESSOR)
   if(nprocs > 1)then
     allocate(gtemp(n1gl,kt))
     if(map_type(1:1) == "e")then
       call acollect(myid,msrid,nprocs,emap,var,gtemp)
     else 
        call acollect(myid,msrid,nprocs,nmap,var,gtemp)
    end if
     allocate(temp(nsta,k1m1))
     if(map_type(1:1) == 'n' .and. TRIM(LOCATION_TYPE) == 'node')then
      do i=1,nsta
       temp(i,1:k1m1) = gtemp(node_sta(i),1:k1m1)
      end do  
     else if(map_type(1:1) == 'e' .and. TRIM(LOCATION_TYPE) == 'cell')then
      do i=1,nsta
       temp(i,1:k1m1) = gtemp(element_sta(i),1:k1m1)
      end do  
     else if(map_type(1:1) == 'n' .and. TRIM(LOCATION_TYPE) == 'cell')then
      do i=1,nsta
       temp(i,:) = 0.0_SP
       do j =1,3
        temp(i,1:k1m1) = temp(i,1:k1m1) + gtemp(nvg(element_sta(i),j),1:k1m1)
       end do
       temp(i,:) = temp(i,:)/3.0_SP
      end do   
     else if(map_type(1:1) == 'e' .and. TRIM(LOCATION_TYPE) == 'node')then
      do i=1,nsta
       temp(i,:) = 0.0_SP
       do j =1, ntvegl(node_sta(i))
        temp(i,1:k1m1) = temp(i,1:k1m1) + gtemp(nbvegl(node_sta(i),j),1:k1m1)
       end do
       temp(i,:) = temp(i,:)/ntvegl(node_sta(i))
      end do   
     end if
     deallocate(gtemp)
   end if
#  endif

!   if(myid /= 1) return
   if(myid == 1) then
     ierr = nf90_put_var(nc_fid,vid,temp,START=dims)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"error writing variable to netcdf")
     end if
   end if  
   deallocate(dims,temp)

   return
 end subroutine PUTVAR2D_REAL
!==============================================================================|

 SUBROUTINE PUTVAR1D_INT(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid&
      &,myid,nprocs,ipt,stk )
   
   !------------------------------------------------------------------------------|
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs&
        &,ipt, stk
   character(len=*),intent(in)   :: map_type
   INTEGER, dimension(i1:i2) :: var
   
   INTEGER, allocatable, dimension(:,:) :: temp
   
   allocate(temp(i1:i2,kt))
   temp(i1:i2,kt)= var
   
   call PUTVAR2D_INT(i1,i2,n1,n1gl,kt,k1,map_type,temp,nc_fid,vid&
        &,myid,nprocs,ipt, stk)
   
   deallocate(temp)
   
 END SUBROUTINE PUTVAR1D_INT
 
 subroutine PUTVAR2D_INT(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid&
      &,myid,nprocs,ipt, stk)
   
   !------------------------------------------------------------------------------|
   
#  if defined (MULTIPROCESSOR)
   use mod_par,only : nmap,emap,acollect
   use lims,only : msrid
#  endif
   use all_vars, only : nvg
   use mod_nctools,only : handle_ncerr
   use mod_utils, only : pstop
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs&
        &,ipt,stk
   character(len=*),intent(in)   :: map_type
   INTEGER, allocatable :: var(:,:)

   INTEGER, allocatable, dimension(:,:) :: temp,gtemp
   integer :: ierr,k1m1,i,j
   integer, allocatable :: dims(:)
   

   k1m1 = k1 
   if(k1m1 == 1)then
     allocate(dims(2))
     dims(1) = 1 
     dims(2) = stk
   else
     allocate(dims(3))
     dims(1) = 1 
     dims(2) = 1 
     dims(3) = stk 
   end if
     

   if(map_type(1:1) /= "e" .and. map_type(1:1) /= "n")then
     write(ipt,*)'map_type input to putvar should be "e" OR "n"'
     call pstop
   end if

   if(nprocs==1)then
     allocate(temp(nsta,k1m1))

     if(map_type(1:1) == 'n' .and. TRIM(LOCATION_TYPE) == 'node')then
      do i=1,nsta
       temp(i,1:k1m1) = var(node_sta(i),1:k1m1)
      end do  
     else if(map_type(1:1) == 'e' .and. TRIM(LOCATION_TYPE) == 'cell')then
      do i=1,nsta
       temp(i,1:k1m1) = var(element_sta(i),1:k1m1)
      end do  
     else if(map_type(1:1) == 'n' .and. TRIM(LOCATION_TYPE) == 'cell')then
      do i=1,nsta
       temp(i,:) = 0.0_SP
       do j =1,3
        temp(i,1:k1m1) = temp(i,1:k1m1) + var(nvg(element_sta(i),j),1:k1m1)
       end do
       temp(i,:) = temp(i,:)/3.0_SP
      end do   
     else if(map_type(1:1) == 'e' .and. TRIM(LOCATION_TYPE) == 'node')then
      do i=1,nsta
       temp(i,:) = 0.0_SP
       do j =1, ntvegl(node_sta(i))
        temp(i,1:k1m1) = temp(i,1:k1m1) + var(nbvegl(node_sta(i),j),1:k1m1)
       end do
       temp(i,:) = temp(i,:)/ntvegl(node_sta(i))
      end do   
     end if
   end if

#  if defined (MULTIPROCESSOR)
   if(nprocs > 1)then
     allocate(gtemp(n1gl,kt))
     if(map_type(1:1) == "e")then
       call acollect(myid,msrid,nprocs,emap,var,gtemp)
     else 
       call acollect(myid,msrid,nprocs,nmap,var,gtemp)
     end if
     allocate(temp(nsta,k1m1))
     if(map_type(1:1) == 'n' .and. TRIM(LOCATION_TYPE) == 'node')then
      do i=1,nsta
       temp(i,1:k1m1) = gtemp(node_sta(i),1:k1m1)
      end do  
     else if(map_type(1:1) == 'e' .and. TRIM(LOCATION_TYPE) == 'cell')then
      do i=1,nsta
       temp(i,1:k1m1) = gtemp(element_sta(i),1:k1m1)
      end do  
     else if(map_type(1:1) == 'n' .and. TRIM(LOCATION_TYPE) == 'cell')then
      do i=1,nsta
       temp(i,:) = 0.0_SP
       do j =1,3
        temp(i,1:k1m1) = temp(i,1:k1m1) + gtemp(nvg(element_sta(i),j),1:k1m1)
       end do
       temp(i,:) = temp(i,:)/3.0_SP
      end do   
     else if(map_type(1:1) == 'e' .and. TRIM(LOCATION_TYPE) == 'node')then
      do i=1,nsta
       temp(i,:) = 0.0_SP
       do j =1, ntvegl(node_sta(i))
        temp(i,1:k1m1) = temp(i,1:k1m1) + gtemp(nbvegl(node_sta(i),j),1:k1m1)
       end do
       temp(i,:) = temp(i,:)/ntvegl(node_sta(i))
      end do   
     end if
     deallocate(gtemp)
   end if
#  endif

!   if(myid /= 1) return
   if(myid == 1) then
     ierr = nf90_put_var(nc_fid,vid,temp,START=dims)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"error writing variable to netcdf")
     end if
   end if  
   deallocate(dims,temp)

   return
 end subroutine PUTVAR2D_INT
!==============================================================================|
!==============================================================================|

 SUBROUTINE PUTVAR_CHAR(i1,i2,var,nc_fid,vid,myid,nprocs,ipt,stk )
   
   !------------------------------------------------------------------------------|
   use mod_nctools,only : handle_ncerr
   implicit none
   integer, intent(in) :: i1,i2,nc_fid,vid,myid,nprocs,ipt, stk
   CHARACTER(LEN=20), dimension(i1:i2) :: var
   CHARACTER(LEN=20), allocatable,dimension(:,:) :: var_tmp
   integer :: ierr
   integer, allocatable :: dims(:)

   allocate(dims(2))
   dims(1) = 1 
   dims(2) = stk
   
   allocate(var_tmp(i1:i2,1))
   var_tmp(i1:i2,1) = var(i1:i2)

   ierr = nf90_put_var(nc_fid,vid,var_tmp,START=dims)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing variable to netcdf")
   end if

   deallocate(dims)

   return
 end subroutine PUTVAR_CHAR
!==============================================================================|
!==============================================================================|
   SUBROUTINE TRIANGLE_GRID_EDGE_GL

!==============================================================================|
!             DEFINE NTVEGL, NBVEGL                                          !
!                                                                              !
! ntvegl(1:mgl):           total number of the surrounding triangles               !
!                      connected to the given node                             !
! nbvegl(1:mgl, 1:ntvegl+1): the identification number of surrounding                !
!                      triangles with a common node (counted clockwise)        !
!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   INTEGER I,J,NCNT,MX_NBR_ELEM_GL

!
!----DETERMINE MAX NUMBER OF SURROUNDING ELEMENTS------------------------------!
!
   MX_NBR_ELEM_GL = 0
   DO I=1,MGL
     NCNT = 0
     DO J=1,NGL
       IF( FLOAT(NVG(J,1)-I)*FLOAT(NVG(J,2)-I)*FLOAT(NVG(J,3)-I) == 0.0_SP) &
         NCNT = NCNT + 1
     END DO
     MX_NBR_ELEM_GL = MAX(MX_NBR_ELEM_GL,NCNT)
   END DO

!
!----ALLOCATE ARRAYS BASED ON MX_NBR_ELEM--------------------------------------!
! 
   ALLOCATE(NTVEGL(MGL));                  NTVEGL = 0
   ALLOCATE(NBVEGL(MGL,MX_NBR_ELEM_GL+1)); NBVEGL = 0
!
!--DETERMINE NUMBER OF SURROUNDING ELEMENTS FOR NODE I = NTVEGL(I)---------------!
!--DETERMINE NBVEGL - INDICES OF NEIGHBORING ELEMENTS OF NODE I------------------!
!
       
   DO I=1,MGL
     NCNT=0
     DO J=1,NGL
       IF (FLOAT(NVG(J,1)-I)*FLOAT(NVG(J,2)-I)*FLOAT(NVG(J,3)-I) == 0.0_SP)THEN
         NCNT = NCNT+1
         NBVEGL(I,NCNT)=J
       END IF
     ENDDO
     NTVEGL(I)=NCNT
   ENDDO

   RETURN
   END SUBROUTINE TRIANGLE_GRID_EDGE_GL
!==============================================================================!
END MODULE MOD_STATION_TIMESERIES



