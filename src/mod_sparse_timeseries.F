MODULE MOD_SPARSE_TIMESERIES
# if defined (WAVE_CURRENT_INTERACTION)

  USE MOD_TIME
  USE NETCDF
  USE VARS_WAVE
  USE MOD_STATION_TIMESERIES, ONLY : PUTVAR,NSTA,NODE_STA
  IMPLICIT NONE
  SAVE
  
  LOGICAL  OUT_WAVE_SPARSE_TIMESERIES_ON
  REAL(SP) SPARSE_DISTANCE
  LOGICAL  OUT_WIND_VELOCITY_SPARSE
  LOGICAL  OUT_SIG_WAVE_HEIGHT_SPARSE
  LOGICAL  OUT_REL_PEAK_PERIOD_SPARSE
  LOGICAL  OUT_WAVE_DIRECTION_SPARSE
  LOGICAL  OUT_ENERGY_SPECTRUM_SPARSE
  LOGICAL  OUT_WAVE_PARTITION_SPARSE
  CHARACTER(LEN=80) OUT_INTERVAL_SPARSE
  
  NAMELIST /NML_WAVE_SPARSE_TIMESERIES/    &
       & OUT_WAVE_SPARSE_TIMESERIES_ON,    &
       & SPARSE_DISTANCE,                  &
       & OUT_WIND_VELOCITY_SPARSE,                &
       & OUT_SIG_WAVE_HEIGHT_SPARSE,              &
       & OUT_REL_PEAK_PERIOD_SPARSE,              &
       & OUT_WAVE_DIRECTION_SPARSE,               &
       & OUT_ENERGY_SPECTRUM_SPARSE,              &
       & OUT_WAVE_PARTITION_SPARSE,               &
       & OUT_INTERVAL_SPARSE
       
!  INTEGER NSTA       
!  INTEGER, ALLOCATABLE  :: NODE_STA(:)
  TYPE(TIME) :: WAVE_INTERVAL_TIME_SERIES, WAVE_TIME_SERIES
  TYPE(TIME) :: KDD,KDD1
  
!--Control Variables----------------------------------------------!
   integer,private :: out_cnt            !!counts number of outputs
   integer,private :: stck_cnt           !!counts number of outputs in each file
   character(len=120),private :: cdfname !!netcdf file name
!--NetCDF IDs----------------------------------------------------!

   !--NetCDF File 
   integer,private :: nc_ofid

   !--Dimensions
   integer,private :: station_did    !,clen_did
   integer,private :: time_did

   !--Grid Variables
   integer,private :: x_s_vid,y_s_vid,lat_s_vid,lon_s_vid

   !--Flow Variables 
   integer,private :: time_s_vid
   integer,private :: iint_vid
   integer,private :: h_s_vid
   integer,private :: uuwind_s_vid
   integer,private :: vvwind_s_vid
   
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

   !--Info Variables
   character(len=120),public :: netcdf_timestring 

  CONTAINS
!====================================================================================
   SUBROUTINE SPARSE_NAME_LIST_INITIALIZE
   
   IMPLICIT NONE
   
   OUT_WAVE_SPARSE_TIMESERIES_ON    = .False.
   SPARSE_DISTANCE                  = 1000.  
   OUT_WIND_VELOCITY_SPARSE         = .False.
   OUT_SIG_WAVE_HEIGHT_SPARSE       = .False. 
   OUT_REL_PEAK_PERIOD_SPARSE       = .False. 
   OUT_WAVE_DIRECTION_SPARSE        = .False.  
   OUT_ENERGY_SPECTRUM_SPARSE       = .False.  
   OUT_WAVE_PARTITION_SPARSE        = .False. 
   OUT_INTERVAL_SPARSE              = "A length of time: 'seconds= ','days= ', or 'cycles= '"
   
   RETURN
   END SUBROUTINE SPARSE_NAME_LIST_INITIALIZE
!----------------------------------------------------------------------------------
   SUBROUTINE SPARSE_NAME_LIST_PRINT
   USE CONTROL, ONLY : IPT
   
   IMPLICIT NONE
   
   WRITE(UNIT=IPT,NML=NML_WAVE_SPARSE_TIMESERIES)
   RETURN
   END SUBROUTINE SPARSE_NAME_LIST_PRINT
!----------------------------------------------------------------------------------   
   SUBROUTINE SPARSE_NAME_LIST_READ
   USE CONTROL, ONLY : casename,NMLUNIT
   USE MOD_UTILS
   USE MOD_SET_TIME, ONLY : GET_OUTPUT_FILE_INTERVAL
   
   IMPLICIT NONE
   
   INTEGER :: IOS, I
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=160) :: PATHNFILE
   
   IF(DBG_SET(DBG_SBR)) &
         & WRITE(IPT,*) "Subroutine Begins: Read_Sparse_Name_List;"

   IOS = 0

   FNAME = "./"//trim(casename)//"_run.nml"

   CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

   !READ NAME LIST FILE
    REWIND(NMLUNIT)

   ! Read IO Information
   READ(UNIT=NMLUNIT, NML=NML_WAVE_SPARSE_TIMESERIES,IOSTAT=ios)
   if(ios .NE. 0 ) Then
     if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_WAVE_SPARSE_TIMESERIES)
     Call Fatal_error("Can Not Read NameList NML_WAVE_SPARSE_TIMESERIES from file: "//trim(FNAME))
   end if
   CLOSE(NMLUNIT)

!   CALL GET_OUTPUT_FILE_INTERVAL(TRIM(OUT_INTERVAL_SPARSE),WAVE_INTERVAL_TIME_SERIES)

   RETURN
   END SUBROUTINE SPARSE_NAME_LIST_READ
!----------------------------------------------------------------------------------   
   SUBROUTINE SPARSE_STATION
   USE LIMS, ONLY : MGL
   IMPLICIT NONE
   INTEGER , ALLOCATABLE :: NODE_STA_TMP(:)   
   INTEGER               :: I
   
   ALLOCATE(NODE_STA(0:MGL))   ;NODE_STA=0 
   CALL OUTPUT_SPARSE_STATION(SPARSE_DISTANCE,NODE_STA)
   DO I = 1,MGL
     IF(NODE_STA(I)==0)THEN
       NSTA=I-1
       EXIT
     END IF
   END DO

   ALLOCATE(NODE_STA_TMP(0:NSTA)) 
   NODE_STA_TMP(0:NSTA) = NODE_STA(0:NSTA)
   DEALLOCATE(NODE_STA)

   ALLOCATE(NODE_STA(0:NSTA))
   NODE_STA = NODE_STA_TMP
   DEALLOCATE(NODE_STA_TMP)

   RETURN
   END SUBROUTINE SPARSE_STATION
!----------------------------------------------------------------------------------
   SUBROUTINE OUT_WAVE_SPARSE_TIMESERIES 
   
   USE NETCDF
   USE SWCOMM3, ONLY : MSC,MDC
   USE MOD_NCTOOLS, ONLY : HANDLE_NCERR
   IMPLICIT NONE
   INTEGER  :: I1,I2,ierr
   REAL(SP) :: THOUR,THOUR1
   INTEGER  :: dims(1)
   REAL(SP) :: KDD_TMP
   REAL(SP), ALLOCATABLE :: ftemp(:)

!------------------------------------------------------------------------------!
!  WRITE TO FILES (SERIAL EXECUTION)                                           !
!------------------------------------------------------------------------------!
   IF(WAVE_TIME_SERIES > IntTime) RETURN
   
   WAVE_TIME_SERIES = IntTime + WAVE_INTERVAL_TIME_SERIES
   THOUR = DTI*FLOAT(IINT-ISTART+1)/3600.0_SP
   THOUR1 = DTI*FLOAT(IINT)/3600.0_SP

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

   IF(OUT_WIND_VELOCITY_SPARSE)THEN 
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

   IF(OUT_SIG_WAVE_HEIGHT_SPARSE)THEN
     i1 = lbound(hsc1,1) ; i2 = ubound(hsc1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",hsc1,nc_ofid,hs_s_vid,myid,nprocs&
          &,ipt, stck_cnt) 
   END IF

   IF(OUT_REL_PEAK_PERIOD_SPARSE)THEN
     i1 = lbound(tpeak,1) ; i2 = ubound(tpeak,1)
     call putvar(i1,i2,m,mgl,1,1,"n",tpeak,nc_ofid,pperiod_s_vid,myid,nprocs&
          &,ipt, stck_cnt) 
   END IF

   IF(OUT_WAVE_DIRECTION_SPARSE)THEN
     i1 = lbound(dirdeg1,1) ; i2 = ubound(dirdeg1,1)
     call putvar(i1,i2,m,mgl,1,1,"n",dirdeg1,nc_ofid,wdir_s_vid,myid,nprocs&
          &,ipt, stck_cnt) 
   END IF

   IF(OUT_ENERGY_SPECTRUM_SPARSE)THEN
     i1 = lbound(spec_density,1) ; i2 = ubound(spec_density,1)
     call putvar(i1,i2,m,mgl,msc,msc,"n",spec_density,nc_ofid,sdensity_s_vid,myid,nprocs&
          &,ipt, stck_cnt) 
   END IF
    
   IF(OUT_WAVE_PARTITION_SPARSE)THEN
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

   IERR = NF90_CLOSE(NC_OFID)
   RETURN
   END SUBROUTINE OUT_WAVE_SPARSE_TIMESERIES
!==============================================================================|

!==============================================================================|
!  Write NetCDF Header and Static Variables                                    |
!==============================================================================|
   SUBROUTINE write_netcdf_setup 

   use mod_clock, only : get_timestamp
   use mod_nctools, only : handle_ncerr
   use swcomm3, only : msc

   use netcdf
   use mod_utils
   implicit none
   integer, dimension(2) :: dynm2ds
   integer, dimension(1) :: stat2ds
   integer, dimension(1) :: dynmtime
   character(len=100)    :: netcdf_convention
   character(len=100)    :: timestamp ,temp
   integer               :: ierr,i1,i2
   integer, dimension(3) :: dynm2dn_msc
   integer, dimension(3) :: dynm2dn_imo

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

    cdfname = trim(OUTPUT_DIR)//trim(casename)//'_sparse_timeseries.nc'

!--Create File 
# if defined (USE_NETCDF4)
    ierr = nf90_create(path=cdfname,cmode=nf90_netcdf4,ncid=nc_ofid)
# else
    ierr = nf90_create(path=cdfname,cmode=nf90_clobber,ncid=nc_ofid)
# endif
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
    ierr = nf90_def_dim(nc_ofid,"station"   ,nsta      ,station_did   )
!    ierr = nf90_def_dim(nc_ofid,"clen", 20, clen_did )
    ierr = nf90_def_dim(nc_ofid,"msc",msc,msc_did)
    ierr = nf90_def_dim(nc_ofid,"IMO",50,imo_did)

!--Define Unlimited Model Dimension
    ierr = nf90_def_dim(nc_ofid,"time"   ,NF90_UNLIMITED,time_did)

!--Set Up Data Dimensioning - Static Vars
    stat2ds      = (/station_did/)            !!2d station vars
!    stat2ds_char = (/clen_did,station_did/)

!--Set Up Data Dimensioning - Dynamic Vars 
    dynm2ds      = (/station_did,time_did/)            !!2d station vars
    dynmtime     = (/time_did/)   
    dynm2dn_msc  = (/station_did,msc_did,time_did/)
    dynm2dn_imo  = (/station_did,imo_did,time_did/)

!--Define Station Name Variables and Attributes

    !!====Station Name (NAME_STA)  ===================!
!    ierr = nf90_def_var(nc_ofid,"name_station",nf90_char,stat2ds_char,name_s_vid)
!    ierr = nf90_put_att(nc_ofid,name_s_vid,"long_name","Station Name")   

!--Define Coordinate Variables and Attributes
    !!====X Grid Coordinate at Nodes (VX) (Meters)===========!
    ierr = nf90_def_var(nc_ofid,"x",nf90_float,stat2ds,x_s_vid)
    ierr = nf90_put_att(nc_ofid,x_s_vid,"long_name","station x-coordinate")
    ierr = nf90_put_att(nc_ofid,x_s_vid,"units","meters")
#   if defined (USE_NETCDF4) && (USE_COMPRESSION)
    ierr = nf90_def_var_deflate(nc_ofid, x_s_vid, shuffle=1, deflate=1, deflate_level=9)
#   endif

    !!====Y Grid Coordinate at Nodes (VY) (Meters)===========!
    ierr = nf90_def_var(nc_ofid,"y",nf90_float,stat2ds,y_s_vid)
    ierr = nf90_put_att(nc_ofid,y_s_vid,"long_name","station y-coordinate")
    ierr = nf90_put_att(nc_ofid,y_s_vid,"units","meters")
#   if defined (USE_NETCDF4) && (USE_COMPRESSION)
    ierr = nf90_def_var_deflate(nc_ofid, y_s_vid, shuffle=1, deflate=1, deflate_level=9)
#   endif

    !!====Longitudinal Coordinate at Nodes (LON) (degrees)===!
    ierr = nf90_def_var(nc_ofid,"lon",nf90_float,stat2ds,lon_s_vid)
    ierr = nf90_put_att(nc_ofid,lon_s_vid,"long_name","Longitude")
    ierr = nf90_put_att(nc_ofid,lon_s_vid,"standard_name","longitude")
    ierr = nf90_put_att(nc_ofid,lon_s_vid,"units","degrees_east")
#   if defined (USE_NETCDF4) && (USE_COMPRESSION)
    ierr = nf90_def_var_deflate(nc_ofid, lon_s_vid, shuffle=1, deflate=1, deflate_level=9)
#   endif

    !!====Latitudinal  Coordinate at Nodes (LAT) (degrees)===!
    ierr = nf90_def_var(nc_ofid,"lat",nf90_float,stat2ds,lat_s_vid)
    ierr = nf90_put_att(nc_ofid,lat_s_vid,"long_name","Latitude")
    ierr = nf90_put_att(nc_ofid,lat_s_vid,"standard_name","latitude")
    ierr = nf90_put_att(nc_ofid,lat_s_vid,"units","degrees_north")
    ierr = nf90_put_att(nc_ofid,lat_s_vid,"grid","Bathymetry_Mesh")
#   if defined (USE_NETCDF4) && (USE_COMPRESSION)
    ierr = nf90_def_var_deflate(nc_ofid, lat_s_vid, shuffle=1, deflate=1, deflate_level=9)
#   endif

!--Define Mesh Relevant Variables and Attributes

    !!====Bathymetry at Nodes (H) (meters)===================!
    ierr = nf90_def_var(nc_ofid,"h",nf90_float,stat2ds,h_s_vid)
    ierr = nf90_put_att(nc_ofid,h_s_vid,"long_name","Bathymetry")   
    ierr = nf90_put_att(nc_ofid,h_s_vid,"units","meters")
    ierr = nf90_put_att(nc_ofid,h_s_vid,"positive","down")
    ierr = nf90_put_att(nc_ofid,h_s_vid,"standard_name","depth")
#   if defined (USE_NETCDF4) && (USE_COMPRESSION)
    ierr = nf90_def_var_deflate(nc_ofid, h_s_vid, shuffle=1, deflate=1, deflate_level=9)
#   endif

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
#   if defined (USE_NETCDF4) && (USE_COMPRESSION)
    ierr = nf90_def_var_deflate(nc_ofid, time_s_vid, shuffle=1, deflate=1, deflate_level=9)
#   endif

    ierr = nf90_def_var(nc_ofid,"iint",nf90_int,dynmtime,iint_vid)
    ierr = nf90_put_att(nc_ofid,iint_vid,"long_name","internal mode iteration number")
#   if defined (USE_NETCDF4) && (USE_COMPRESSION)
    ierr = nf90_def_var_deflate(nc_ofid, iint_vid, shuffle=1, deflate=1, deflate_level=9)
#   endif

!--Define Time Dependent Flow Variables (selected by user from input file)

    if(OUT_WIND_VELOCITY_SPARSE)then
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

#    if defined (USE_NETCDF4) && (USE_COMPRESSION)
     ierr = nf90_def_var_deflate(nc_ofid, uuwind_s_vid, shuffle=1, deflate=1, deflate_level=9)
     ierr = nf90_def_var_deflate(nc_ofid, vvwind_s_vid, shuffle=1, deflate=1, deflate_level=9)
#    endif
    end if


    if(OUT_SIG_WAVE_HEIGHT_SPARSE)then
     ierr = nf90_def_var(nc_ofid,"wave_height",nf90_float,dynm2ds,hs_s_vid)
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"long_name","Signaficant Wave Height")
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"standard_name","signaficant_wave_height")
     ierr = nf90_put_att(nc_ofid,hs_s_vid,"type","data")
#    if defined (USE_NETCDF4) && (USE_COMPRESSION)
     ierr = nf90_def_var_deflate(nc_ofid, hs_s_vid, shuffle=1, deflate=1, deflate_level=9)
#    endif
    end if

    if(OUT_REL_PEAK_PERIOD_SPARSE)then
     ierr = nf90_def_var(nc_ofid,"wave_period",nf90_float,dynm2ds,pperiod_s_vid)
     ierr = nf90_put_att(nc_ofid,pperiod_s_vid,"long_name","Relative Peak Period")
     ierr = nf90_put_att(nc_ofid,pperiod_s_vid,"units","seconds")
     ierr = nf90_put_att(nc_ofid,pperiod_s_vid,"standard_name","relative_peak_period")
     ierr = nf90_put_att(nc_ofid,pperiod_s_vid,"type","data")
#    if defined (USE_NETCDF4) && (USE_COMPRESSION)
     ierr = nf90_def_var_deflate(nc_ofid, pperiod_s_vid, shuffle=1, deflate=1, deflate_level=9)
#    endif
    end if

    if(OUT_WAVE_DIRECTION_SPARSE)then
     ierr = nf90_def_var(nc_ofid,"wave_dir",nf90_float,dynm2ds,wdir_s_vid)
     ierr = nf90_put_att(nc_ofid,wdir_s_vid,"long_name","Wave Direction")
     ierr = nf90_put_att(nc_ofid,wdir_s_vid,"units","degrees")
     ierr = nf90_put_att(nc_ofid,wdir_s_vid,"standard_name","wave_direction")
     ierr = nf90_put_att(nc_ofid,wdir_s_vid,"type","data")
#    if defined (USE_NETCDF4) && (USE_COMPRESSION)
     ierr = nf90_def_var_deflate(nc_ofid, wdir_s_vid, shuffle=1, deflate=1, deflate_level=9)
#    endif
    end if

    if(OUT_ENERGY_SPECTRUM_SPARSE)then
     ierr = nf90_def_var(nc_ofid,"energy_spectrum",nf90_float,dynm2dn_msc,sdensity_s_vid)
     ierr = nf90_put_att(nc_ofid,sdensity_s_vid,"long_name","Energy Spectrum")
     ierr = nf90_put_att(nc_ofid,sdensity_s_vid,"units","m^2/hz")
     ierr = nf90_put_att(nc_ofid,sdensity_s_vid,"standard_name","energy_dpectrum")
     ierr = nf90_put_att(nc_ofid,sdensity_s_vid,"type","data")
#    if defined (USE_NETCDF4) && (USE_COMPRESSION)
     ierr = nf90_def_var_deflate(nc_ofid, sdensity_s_vid, shuffle=1, deflate=1, deflate_level=9)
#    endif
    end if

    if(OUT_WAVE_PARTITION_SPARSE)then
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

#    if defined (USE_NETCDF4) && (USE_COMPRESSION)
     ierr = nf90_def_var_deflate(nc_ofid, hwind_s_vid, shuffle=1, deflate=1, deflate_level=9)
     ierr = nf90_def_var_deflate(nc_ofid, dwind_s_vid, shuffle=1, deflate=1, deflate_level=9)
     ierr = nf90_def_var_deflate(nc_ofid, twind_s_vid, shuffle=1, deflate=1, deflate_level=9)
     ierr = nf90_def_var_deflate(nc_ofid, twindp_s_vid, shuffle=1, deflate=1, deflate_level=9)
     ierr = nf90_def_var_deflate(nc_ofid, hswell_s_vid, shuffle=1, deflate=1, deflate_level=9)
     ierr = nf90_def_var_deflate(nc_ofid, dswell_s_vid, shuffle=1, deflate=1, deflate_level=9)
     ierr = nf90_def_var_deflate(nc_ofid, tswell_s_vid, shuffle=1, deflate=1, deflate_level=9)
     ierr = nf90_def_var_deflate(nc_ofid, tswellp_s_vid, shuffle=1, deflate=1, deflate_level=9)
#    endif
    end if

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
   call putvar(i1,i2,m,mgl,1,1,"n",lon,nc_ofid,lon_s_vid,myid,nprocs,ipt, stck_cnt)

   !!====Latitude  at Nodes (LAT) ==========================!
   i1 = lbound(lat,1) ; i2 = ubound(lat,1)
   call putvar(i1,i2,m,mgl,1,1,"n",lat,nc_ofid,lat_s_vid,myid,nprocs,ipt, stck_cnt) 

   !!====X Grid Coordinate at Nodes (VX)====================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar(i1,i2,m,mgl,1,1,"n",vx+vxmin,nc_ofid,x_s_vid,myid,nprocs,ipt, stck_cnt) 

   !!====Y Grid Coordinate at Nodes (VY)====================!
   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar(i1,i2,m,mgl,1,1,"n",vy+vymin,nc_ofid,y_s_vid,myid,nprocs,ipt, stck_cnt) 

   !!====Bathymetry at Nodes (H)============================!
   i1 = lbound(h,1) ; i2 = ubound(h,1)
   call putvar(i1,i2,m,mgl,1,1,"n",h,nc_ofid,h_s_vid,myid,nprocs,ipt,stck_cnt) 

!==============================================================================|
!  close the file                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(nc_ofid)
   return
   end subroutine write_netcdf_setup
 
!==============================================================================|
!==============================================================================|

   SUBROUTINE OUTPUT_SPARSE_STATION(OUT_DIS,OUT_NODES)

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR , ONLY : NMAP,ACOLLECT
#  endif
   IMPLICIT NONE
   INTEGER  :: I,OUT_NODES(0:MGL),IERR
   REAL(SP) :: OUT_DIS
   REAL(SP), ALLOCATABLE :: VX_ALL(:),VY_ALL(:) 

!==============================================================================|

   IF(SERIAL)THEN
     OUT_NODES = 0
     CALL SPARSH_GRID(VX,VY,OUT_DIS,OUT_NODES)
   END IF 

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     ALLOCATE(VX_ALL(0:MGL))
     ALLOCATE(VY_ALL(0:MGL))

     CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VX,VX_ALL)
     CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,VY,VY_ALL)

     OUT_NODES = 0
     CALL SPARSH_GRID(VX_ALL,VY_ALL,OUT_DIS,OUT_NODES)
     CALL MPI_BCAST(OUT_NODES,MGL+1,MPI_INTEGER,0,MPI_FVCOM_GROUP,IERR)

     DEALLOCATE(VX_ALL)
     DEALLOCATE(VY_ALL)
   END IF
#  endif

   RETURN
   END SUBROUTINE OUTPUT_SPARSE_STATION
!==============================================================================|

!==============================================================================|
   SUBROUTINE SPARSH_GRID(VX_ALL,VY_ALL,OUT_DIS,OUT_NODES)
   USE VARS_WAVE
   IMPLICIT NONE
   REAL(SP) :: XCC(0:MGL),YCC(0:MGL),XX,YY,OUT_DIS
   INTEGER  :: RPLOT,RMIN,NP(0:MGL),OUT_NODES(0:MGL)
   INTEGER  :: I,II,IP

   REAL(SP) :: VX_ALL(0:MGL),VY_ALL(0:MGL)
   REAL(SP) :: R
   INTEGER,PARAMETER :: DIR = 1 ! LL--->XY
   INTEGER,PARAMETER :: ALPHAA =1 ! NORTH POLE

#  if defined(SPHERICAL)
   DO I =0,MGL
     CALL POLAR_STEREO(VX_ALL(I),VY_ALL(I),XCC(I),YCC(I),ALPHAA,DIR)
   END DO
#  else
   DO I =0,MGL
     XCC(I)=VX_ALL(I)
     YCC(I)=VY_ALL(I)
   END DO
#  endif
   RPLOT = OUT_DIS
   IP = 0
   DO I = 1,MGL
     XX = XCC(I)
     YY = YCC(I)
     RMIN = RPLOT+1
     DO II = 1,IP
       R = SQRT((XX-XCC(NP(II)))**2+(YY-YCC(NP(II)))**2)
       IF(R < RMIN) RMIN = R
     END DO
     IF(RMIN > RPLOT)THEN
       OUT_NODES(IP) = I
       IP = IP+1
       NP(IP) = I
     END IF
   END DO

   RETURN
   END SUBROUTINE SPARSH_GRID


!=============================================================
   SUBROUTINE POLAR_STEREO(LONG,LAT,X,Y,ALPHAA,DIR)

   IMPLICIT NONE
   REAL(SP), PARAMETER :: PI        = 3.1415926
   REAL(SP), PARAMETER :: R_EARTH   = 6371.E+3
   REAL(SP), PARAMETER :: DEG_2_RAD = .017453292519943
   REAL(SP), PARAMETER :: RAD_2_DEG = 57.2957795130823

   REAL(SP) :: LONG,LAT,X,Y ,RHO
   REAL(SP) :: MAP_FACTOR, RLAMP
   INTEGER  :: DIR,ALPHAA

   IF(DIR ==1)THEN

     MAP_FACTOR = 2.0/(1+ALPHAA*SIN(LAT*DEG_2_RAD))

!    Convert from geodetic latitude and longitude to polar stereographic      
     RHO = R_EARTH * MAP_FACTOR * COS(LAT*DEG_2_RAD)

     RLAMP=(LONG)*DEG_2_RAD

     X = RHO * COS(RLAMP)
     Y = RHO * SIN(RLAMP)

   ELSE

     LAT =ALPHAA* asin((4*R_EARTH*R_EARTH - X*X - Y*Y)   &
          /(4*R_EARTH*R_EARTH + X*X + Y*Y))
     LONG =ATAN2(Y,X)

     LAT = LAT*RAD_2_DEG
     LONG= LONG*RAD_2_DEG

   ENDIF

   RETURN

   END SUBROUTINE POLAR_STEREO
!=============================================================

# endif  
END MODULE MOD_SPARSE_TIMESERIES



