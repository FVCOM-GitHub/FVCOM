!===============================Make Executable============================
!  Make executable:
!
!   ifort wrf_to_fvcom.f90 -L/hosts/salmon01/data00/medm/netcdf/3.6.2/em64t/lib -lnetcdf -lm -I/hosts/salmon01/data00/medm/netcdf/3.6.2/em64t/include -FR -o wrf2fvcom
!
!
!==========================================================================
!! This program read a WRF netcdf file and create forcing data for fvcom 
!! version 0.13             2007/06/15
!! updated using  coare26sn
!! version 0.12             2007/06/14
!  add flag -global           to check if we want to save global attributes or not
!! version 0.11             2007/06/12
!! temporarily for lin's case, output forcing and save file including SLP
!! This is for ocean use, the sea level pressure is approximately 
!! set as (P+PB)/100 mb  at the lowest layer.
!
!===============================Run Program================================
!  Run program:
!
!   wrf2fvcom -i wrf_netcdf_input_file -o fvcom_forcing_file \
!   -forecast -s surface_meterological_variables_file         
!
!===============================Options====================================
!
!   -help     : Print this information
!   -h        : Print this information
!   -i        : WRF input file (netcdf format)
!   -o        : forcing output name - default is ARWout
!   -s        : Extract and save necessary fields for backup
!   -debug    : Print some debug information
!   -forecast : flag for forecast data
!   -hindcast : flag for hindcast data
!   -noglobal : don't save global attributes   
!
!==========================================================================
!
!  June 2007  VERSION 0.13
!  Song Hu 
!
!==========================================================================
!
!  MODIFICATIONS:
!
!   Changed name of the time dimension to 'Time' - the same as
! original wrf files: David Stuebe 6/27/07
!
!   Changed the definition of the time dimension to 'NF90_UNLIMITED'
! instead of the fixed length, 73: David Stuebe 6/27/07
!
!   Changed the code to get the length of the record dimension from
! the input file: Song Hu 6/27/07
!
!   Changed the code to get the correct dimension for -s option: 
!   Song Hu  7/02/07


	program wrf2fvcom 

	use netcdf
	implicit none
	include 'netcdf.inc'

!--------------Variable for reading input netcdf file--------------------------------------------------------------------------------------
	integer                                             :: cdfid                       ! id of input file
	real                                                :: rcall                       ! status code for netcdf
	character (len=200)                                 :: attname                     ! name of attributes
	integer                                             :: ndims, nvars, natts         ! number of dims, var and att in file
	integer                                             :: south_north,west_east       ! dimension of horizontal grids
	integer                                             :: bottom_n                    ! dimension of vertical grid
	integer                                             :: DateStr_id                  ! string length of time
        integer                                             :: unlimdimid                  ! unlimited dims ID - not used
	character (len=31), allocatable, dimension(:)       :: dnam                        ! name of netcdf DIMS
	integer,            allocatable, dimension(:)       :: dval                        ! value of netcdf DIMS
        integer                                             :: iweg, isng, ibtg            ! staggered dims in input_file
	character (len=1)                                   :: gridtype                    ! global att from netcdf file
	real                                                :: stand_lon                   ! global att from netcdf file
	character (len=40)                                  :: title                       ! global att from netcdf file
	real                                                :: truelat1, truelat2          ! global att from netcdf file
	real                                                :: dx, dy                      ! global att from netcdf file
	real                                                :: cen_lon, cen_lat            ! global att from netcdf file
	integer                                             :: map_proj                    ! global att from netcdf file
	character (len=80)                                  :: varnam                      ! var from netcdf file
	integer                                             :: idvar                       ! varnam id number
	integer                                             :: ivtype                      ! varnam type
	integer                                             :: idm                         ! varnam - number of dims
	integer                                             :: natt                        ! varnam attributes
	character (len=80)                                  :: description                 ! varnam description
	character (len=80)                                  :: units                       ! varnam units
	integer                                             :: len_unt, len_title          ! length of character strings
	integer                                             :: len_string, len_des         ! length of character strings
	integer                                             :: len_var                     ! length of character strings
	integer                                             :: dims(4), ishape(10)         ! netcdf array dims and shapes

        integer                                             :: jjj,j1,j2,j3,j4             ! tmp variables for dimension
	integer                                             :: i,j,k                         ! tmp variables for loop
!----------------------Input and output forcing file---------------------------------------------------------------------------------------------------------
        integer                                             :: nc_ofid                     ! id of -o and -s  net cdf files

        character (len=200)                                 :: input_file                  ! netcdf input file (ARW WRF raw output) or product from -s save_file
        character (len=200)                                 :: save_file                   ! output meteorological for backup (only keep necessary fields to save space)
        character (len=200)                                 :: history,source,institute,gti! output global attributes


!----------------------WRF Meteorological Variables---------------------------------------------------------------------------------------
        integer                                             :: south_north_id              ! id of output horizontal dimension
        integer                                             :: west_east_id                ! id of output horizontal dimension
        integer                                             :: bottom_up_id                ! id of output vertical   dimension
        integer                                             :: time_id                     ! id of output time
        integer, dimension(2)                               :: stat0m                      ! dimension for 2D variable (Time String)
        integer, dimension(3)                               :: stat2m                      ! dimension for surface 2D variable (except P and PB)
        integer, dimension(4)                               :: stat3m                      ! dimensino for 3D varialbe (P,PB)

!---------id for meteorological variables---------------block
        integer                                             :: t2_id,q2_id,sst_id          ! id of air temperature, qv, sst
	integer                                             :: u10_id,v10_id,tsk_id        ! id of u10,v10 and skin temperature
        integer                                             :: rainc_id,rainnc_id          ! id of rainc and rain-nonc
        integer                                             :: glw_id,swdown_id            ! id of longwave(atmosphere) and shortwave
        integer                                             :: p_id,pb_id                  ! id of P and PB
        integer                                             :: xlat_id,xlong_id            ! id of Latitude and Longitude
        integer                                             :: temp_id                     ! temp id
!--------array for meteorological variables-------------block
	character,          allocatable, dimension(:,:)     :: times                       ! Times array from netcdf file
	real,               allocatable, dimension(:,:,:,:) :: xlat                        ! Latitude array from netcdf file
	real,               allocatable, dimension(:,:,:,:) :: xlong                       ! Longitud array from netcdf file  
	real,               allocatable, dimension(:,:,:,:) :: p                           ! Pressure pertubation from netcdf (4 dimensions)
	real,               allocatable, dimension(:,:,:,:) :: pb                          ! Pressure base from netcdf        (4 dimensions)
	real,               allocatable, dimension(:,:,:,:) :: ptemp                       ! Pressure pertubation from netcdf (3 dimensions)
	real,               allocatable, dimension(:,:,:,:) :: pbtemp                      ! Pressure base from netcdf        (3 dimensions)
	real,               allocatable, dimension(:,:,:,:) :: q2                          ! QV from netcdf 
	real,               allocatable, dimension(:,:,:,:) :: sst                         ! Sea Surface Temperature
	real,               allocatable, dimension(:,:,:,:) :: swdown                      ! DOWNWARD SHORTWAVE
	real,               allocatable, dimension(:,:,:,:) :: glw                         ! DOWNWARD LONGWAVE (meteorological)
	real,               allocatable, dimension(:,:,:,:) :: t2                          ! T at 2 meter
	real,               allocatable, dimension(:,:,:,:) :: u10                         ! U at 10 meter
	real,               allocatable, dimension(:,:,:,:) :: v10                         ! V at 10 meter
	real,               allocatable, dimension(:,:,:,:) :: rainc                       ! rainc of WRF 
	real,               allocatable, dimension(:,:,:,:) :: rainnc                      ! rainnc of WRF
	real,               allocatable, dimension(:,:,:,:) :: tsk                         ! Skin surface temperauter( we need this because SST only availalbe on ocean)

!-------------for time string output-----------------------------------------------------------------------------------------------------------
        integer                                             :: ntimes                      ! number of times in input file
        integer, dimension(2)                               :: vstart                      ! start index of time string array
        integer, dimension(2)                               :: vcount                      ! end   index of time string array

!------------for read arguments part------------------------------------------------------------------------------------------------------------
	logical                                             :: noglobal                    ! for global attributes
	logical                                             :: debug                       ! for debug
	logical                                             :: small                       ! for save small meteorological data
	logical                                             :: output                      ! for output of forcing

!----------------------forcing---------coare2.6 part-----------------------------------------------------------------------------
        character (len=80)                                  :: case                        ! forcing file name
        character (len=80)                                  :: output_M                       ! forcing file name
        character (len=200), allocatable, dimension(:)      :: M_vars                      ! variables to forcing
	integer                                             :: time_index                  ! time dimension
	integer,dimension(3)                                :: force3m                     ! forcing file dimension
	integer,dimension(2)                                :: force2m                     ! forcing file dimension
	integer,dimension(2)                                :: force1m                     ! forcing file dimension

!----------------------------id of forcing variables----------------------
        integer                                             :: stress_u_id,stress_v_id
        integer                                             :: netheat_id,shortwave_id,longwave_id,sensible_id,latent_id,evap_id,prec_id
        integer                                             :: pressure_slp_id, rh_gen_id
!---------------------------array of forcing variables-------------------
        real,               allocatable, dimension(:,:,:)   :: stress_u                    ! stress u
        real,               allocatable, dimension(:,:,:)   :: stress_v                    ! stress v
        real,               allocatable, dimension(:,:,:)   :: netheat                     ! Net Heat Flux
        real,               allocatable, dimension(:,:,:)   :: longwave                    ! Long wave
        real,               allocatable, dimension(:,:,:)   :: sensible                    ! Sensible Heat flux
        real,               allocatable, dimension(:,:,:)   :: latent                      ! Latent Heat flux
        real,               allocatable, dimension(:,:,:)   :: evaporation                 ! Evaporation
        real,               allocatable, dimension(:,:,:)   :: precipitation               ! Precipitation
        real,               allocatable, dimension(:,:,:)   :: pressure_slp                ! pressure_slp
	real,               allocatable, dimension(:,:,:)   :: rh_gen                      ! Relative humidity generated by QV

!--------------------for coare 2.6 bulk-------------------------------------------------------------------------------------------
	real                                                :: urxx,taxx,paxx,tsxx,zuxx,ztxx,zqxx,rh_valxx
	real                                                :: dswxx, dlwxx, tauxx,hsbxx, hlbxx
	real                                                :: theta, AG, gb_lat,gb_zi
!-----------------END OF VARIABLE DEFINITION---------------------------------------------------------------------------------------



!-------------Generate source information-----------------------------------------------------------
	call date_and_time(gti)
        source ="wrf2fvcom version 0.13 (2007-06-27) (Bulk method: COARE 2.6SN)"
	institute = "School of Marine Science and Technology, UMASSD, at time of "//trim(gti) 
        write(*,*) "  "
        write(*,*) "================================================ "
        write(*,*) trim(source)




! GET INPUT FILE AND OUTPUT CASE NAME FOR COMMAND LINE
	if (debug) write(*,*) "READ arguments from command line "
	call read_args(input_file,case,debug,small,output,save_file,history,noglobal)

!! OPEN THE INPUT FILE
	if (debug) write(*,*) "OPENING netcdf input file "//trim(input_file)
	rcall = nf_open (input_file, NF_NOWRITE, cdfid)
	call handle_err ("Opening input netCDF file", rcall)
	rcall = nf_get_att_text(cdfid, nf_global, 'GRIDTYPE', gridtype)
	if ( trim(gridtype) .ne. 'C' ) then
	write(*,*) 'Warning: the input file is not the ARW raw output'
	write(*,*) 'It probably is the backup meteorological data'
!      STOP
	endif



!! GET BASIC INFORMTION ABOUT THE FILE
	if (debug) write(*,*) "READ global attributes from netcdf file"
	stand_lon = -99999.99   !!! backward compatibility
	rcall = nf_get_att_text(cdfid, nf_global, 'TITLE',     title)
	rcall = nf_get_att_real(cdfid, nf_global, 'DX',        dx)
	rcall = nf_get_att_real(cdfid, nf_global, 'DY',        dy)
	rcall = nf_get_att_real(cdfid, nf_global, 'CEN_LAT',   cen_lat)
	rcall = nf_get_att_real(cdfid, nf_global, 'CEN_LON',   cen_lon)
	rcall = nf_get_att_real(cdfid, nf_global, 'STAND_LON', stand_lon)
	rcall = nf_get_att_real(cdfid, nf_global, 'TRUELAT1',  truelat1)
	rcall = nf_get_att_real(cdfid, nf_global, 'TRUELAT2',  truelat2)
	rcall = nf_get_att_int (cdfid, nf_global, 'MAP_PROJ',  map_proj)
	if (stand_lon == -99999.99 ) stand_lon = cen_lat    !!! dealing with an old WRF file

	if (debug) write(*,*) "READ dims from netcdf file"
	rcall = nf_inq(cdfid, nDims, nVars, nAtts, unlimDimID)
	allocate (M_vars(nVars))
	allocate (dval(nDims))
	allocate (dnam(nDims))
	iweg = 1
	isng = 1
	ibtg = 1
	do i = 1,nDims
	rcall = nf_inq_dim(cdfid, i, dnam(i), dval(i))


	if( dnam(i)(1:12) == 'south_north ')  south_north=dval(i)
        if( dnam(i)(1:12) == 'west_east   ')  west_east  =dval(i)
	if( dnam(i)(1:12) == 'bottom_top  ')  bottom_n   =i


	if ( dnam(i)(1:9)  == 'west_east'          ) iweg = max(iweg,dval(i))
	if ( dnam(i)(1:11) == 'south_north'        ) isng = max(isng,dval(i))
	if ( dnam(i)(1:10) == 'bottom_top'         ) ibtg = max(ibtg,dval(i))
	if ( trim(dnam(i)) == 'num_metgrid_levels' ) ibtg = max(ibtg,dval(i))
	if ( dnam(i)(1:4)  == 'land'               ) ibtg = max(ibtg,dval(i))
	if ( dnam(i)(1:4)  == 'soil'               ) ibtg = max(ibtg,dval(i))
	if ( dnam(i)(1:5)  == 'month'              ) ibtg = max(ibtg,dval(i))
	if ( dnam(i)(1:5)  == 'z-dim'              ) ibtg = max(ibtg,dval(i))
	enddo




        DO idvar = 1,nVars        !! LOOP THROUGH ALL VARIABLES IN FILE
        rcall = nf_inq_var(cdfid, idvar, varnam, ivtype, idm, ishape, natt)


        call handle_err ("ERROR reading variable", rcall)
        if (debug) write(*,*)
        if (debug) write(*,*) "Variable:",idvar,"out of",nVars
        if (debug) write(*,*) "DEALING with variable: ", trim(varnam)
	do i=1,idm
	enddo
        if (.not. debug ) write(*,'("   Variable ",i2," :  ",A10,$)') idvar,varnam
        if (.not. debug .and. idm .gt. 2 ) write(*,*) " (*)"
        if (.not. debug .and. idm .le. 2 ) write(*,*) "    "

!! GET THE DIMS FOR INPUT AND OUTPUT FROM THE SHAPE
        dims  = 1
        do i = 1,idm
          dims(i)  = dval(ishape(i))
        enddo
        if (debug) write(*,*) " DIMS: ", dims
!! GET THE UNITS AND DESCRIPTION OF THE FIELD
        units = "                                                                                "
        description = "                                                                                "
        rcall = NF_GET_ATT_TEXT(cdfid, idvar, "units", units )
        rcall = NF_GET_ATT_TEXT(cdfid, idvar, "description", description )
        len_var = len_trim(varnam)
        len_des = len_trim(description)
        len_unt = len_trim(units)
        if (debug) write(*,*) " DESCRIPTION: ", description(1:len_des)
        if (debug) write(*,*) " UNITS: ", units(1:len_unt)

!! GET THE NUMBER OF TIMES 

        IF (varnam == 'Times' ) THEN
           allocate (times(dims(1), dims(2)))
           rcall = nf_get_var_text(cdfid, idvar, times)
           ntimes = dims(2)
	write(*,*)dims(1),dims(2)

        ENDIF

        IF (varnam == 'XLAT'  ) THEN
           allocate (xlat(dims(1), dims(2), dims(3), dims(4)))
	write(*,*) trim(varnam)
	write(*,*) dims(1),dims(2),dims(3),dims(4)
           rcall = nf_get_var_real(cdfid, idvar, xlat)
        ENDIF


        IF (varnam == 'XLONG'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)

           allocate (xlong(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, xlong)
        ENDIF

        IF (varnam == 'P'  ) THEN          !! save the lowest P to ptemp
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)

           allocate (p(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, p)

	do jjj=1,4
	    if(ishape(jjj)==bottom_n)  dims(jjj)=1
	enddo
	   allocate (ptemp(dims(1), dims(2), dims(3), dims(4)))
           do j4=1,dims(4)
	   do j3=1,dims(3)
           do j2=1,dims(2)
	   do j1=1,dims(1)
           ptemp(j1,j2,j3,j4)=p(j1,j2,j3,j4)
	   enddo
	   enddo
           enddo
           enddo


        ENDIF






        IF (varnam == 'PB'  ) THEN  !! save the lowest PB to pbtemp
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (pb(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, pb)

        do jjj=1,4
            if(ishape(jjj)==bottom_n)  dims(jjj)=1
        enddo
           allocate (pbtemp(dims(1), dims(2), dims(3), dims(4)))
           do j4=1,dims(4)
           do j3=1,dims(3)
           do j2=1,dims(2)
           do j1=1,dims(1)
           pbtemp(j1,j2,j3,j4)=pb(j1,j2,j3,j4)
           enddo
           enddo
           enddo
           enddo
        ENDIF


        IF (varnam == 'Q2'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (q2(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, q2)
        ENDIF

        IF (varnam == 'SST'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (sst(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, sst)
        ENDIF

        IF (varnam == 'TSK'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (tsk(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, tsk)
        ENDIF




        IF (varnam == 'SWDOWN'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (swdown(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, swdown)
        ENDIF


        IF (varnam == 'GLW'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (glw(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, glw)
        ENDIF



        IF (varnam == 'T2'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (t2(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, t2)
        ENDIF


        IF (varnam == 'U10'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (u10(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, u10)
        ENDIF


        IF (varnam == 'V10'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (v10(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, v10)
        ENDIF

        IF (varnam == 'RAINC'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (rainc(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, rainc)
        ENDIF



        IF (varnam == 'RAINNC'  ) THEN
        write(*,*) trim(varnam)
        write(*,*) dims(1),dims(2),dims(3),dims(4)
           allocate (rainnc(dims(1), dims(2), dims(3), dims(4)))
           rcall = nf_get_var_real(cdfid, idvar, rainnc)
        ENDIF



	enddo             !! LOOP of Variables


!!     Extract neccessary meteorological field and save in a smaller file
	if(small)then               !! create netcdf file to store backup
	rcall = nf90_create(path=trim(save_file),cmode=nf90_clobber,ncid=nc_ofid)
	call handle_err ("Create small netCDF file", rcall)

!--Define Fixed Model Dimensions
        rcall = nf90_def_dim(nc_ofid,"south_north",south_north,south_north_id)
        rcall = nf90_def_dim(nc_ofid,"west_east",west_east,west_east_id)
        rcall = nf90_def_dim(nc_ofid,"bottom_up",1,bottom_up_id)
        rcall = nf90_def_dim(nc_ofid,"DateStrLen",19,DateStr_id)
!        rcall = nf90_def_dim(nc_ofid,"Times",73,time_id)

!        rcall = nf90_def_dim(nc_ofid,"Time",73,time_id)

        !USE UNLIMITED DIMENSION FOR TIME
      rcall = nf90_def_dim(nc_ofid,"Time",nf90_unlimited,time_id)

!--Set up static dimensions
       
       stat3m=(/west_east_id,south_north_id,bottom_up_id,time_id/)  ! for P and PB
       stat2m=(/west_east_id,south_north_id,time_id/)
       stat0m=(/DateStr_id,time_id/)                             ! for time
!--Define Variables and Attributes

	rcall = nf90_def_var(nc_ofid,"XLAT",nf90_float,stat2m,xlat_id)
        rcall =nf90_inq_varid(cdfid,"XLAT",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,xlat_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,xlat_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,xlat_id)


        rcall = nf90_def_var(nc_ofid,"XLONG",nf90_float,stat2m,xlong_id)
        rcall =nf90_inq_varid(cdfid,"XLONG",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,xlong_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,xlong_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,xlong_id)


        rcall=nf90_def_var(nc_ofid,"T2",nf90_float,stat2m,t2_id)
        rcall =nf90_inq_varid(cdfid,"T2",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,t2_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,t2_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,t2_id)


        rcall = nf90_def_var(nc_ofid,"Q2",nf90_float,stat2m,q2_id)
        rcall =nf90_inq_varid(cdfid,"Q2",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,q2_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,q2_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,q2_id)


        rcall = nf90_def_var(nc_ofid,"SST",nf90_float,stat2m,sst_id)
        rcall =nf90_inq_varid(cdfid,"SST",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,sst_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,sst_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,sst_id)


        rcall = nf90_def_var(nc_ofid,"TSK",nf90_float,stat2m,tsk_id)
        rcall =nf90_inq_varid(cdfid,"TSK",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,tsk_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,tsk_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,tsk_id)



        rcall = nf90_def_var(nc_ofid,"U10",nf90_float,stat2m,u10_id)
        rcall =nf90_inq_varid(cdfid,"U10",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,u10_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,u10_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,u10_id)


        rcall = nf90_def_var(nc_ofid,"V10",nf90_float,stat2m,v10_id)
        rcall =nf90_inq_varid(cdfid,"V10",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,v10_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,v10_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,v10_id)

        rcall = nf90_def_var(nc_ofid,"RAINC",nf90_float,stat2m,rainc_id)
        rcall =nf90_inq_varid(cdfid,"RAINC",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,rainc_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,rainc_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,rainc_id)

        rcall = nf90_def_var(nc_ofid,"RAINNC",nf90_float,stat2m,rainnc_id)
        rcall =nf90_inq_varid(cdfid,"RAINNC",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,rainnc_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,rainnc_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,rainnc_id)


        rcall = nf90_def_var(nc_ofid,"GLW",nf90_float,stat2m,glw_id)
        rcall =nf90_inq_varid(cdfid,"GLW",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,glw_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,glw_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,glw_id)

        rcall = nf90_def_var(nc_ofid,"SWDOWN",nf90_float,stat2m,swdown_id)
        rcall =nf90_inq_varid(cdfid,"SWDOWN",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,swdown_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,swdown_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,swdown_id)

        rcall = nf90_def_var(nc_ofid,"P ",nf90_float,stat3m,p_id)
        rcall =nf90_inq_varid(cdfid,"P ",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,p_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,p_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,p_id)

        rcall = nf90_def_var(nc_ofid,"PB",nf90_float,stat3m,pb_id)
        rcall =nf90_inq_varid(cdfid,"PB",temp_id)
        rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,pb_id)
        rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,pb_id)
        rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,pb_id)

        rcall=nf90_def_var(nc_ofid,"Times",nf90_char,stat0m,time_id)
        rcall=nf90_put_att(nc_ofid,time_id,"description","GMT time") 
 

!--------Write global attributes

	if(.NOT. noglobal)then
	do k=1,nAtts
	rcall=nf90_inq_attname(cdfid,nf_global,k,attname)
         call handle_err ("Error reading Global Attributes", rcall)
	rcall=nf90_copy_att(cdfid,nf_global,trim(attname),nc_ofid,nf90_global)
         call handle_err ("Error writing Global Attributes", rcall)
	enddo
	endif
	rcall=nf90_put_att(nc_ofid,nf90_global,"History",history)
        rcall=nf90_put_att(nc_ofid,nf90_global,"Source",source)
	rcall=nf90_put_att(nc_ofid,nf90_global,"Institute",institute)

!--Exit Define Mode
	rcall=nf90_enddef(nc_ofid)



!-- Write Data--------------------------------

      rcall =nf90_put_var(nc_ofid,xlat_id,xlat(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,xlong_id,xlong(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,t2_id,t2(:,:,:,1))
      call handle_err ("Error writing t2", rcall)

      rcall =nf90_put_var(nc_ofid,q2_id,q2(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,sst_id,sst(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,tsk_id,tsk(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,u10_id,u10(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,v10_id,v10(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,rainc_id,rainc(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,rainnc_id,rainnc(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,glw_id,glw(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,swdown_id,swdown(:,:,:,1))

      rcall =nf90_put_var(nc_ofid,p_id,ptemp(:,:,:,:))

      rcall =nf90_put_var(nc_ofid,pb_id,pbtemp(:,:,:,:))     
  
      vstart(1)=1
      vstart(2)=1
      vcount(1)=19
      vcount(2)=ntimes
      rcall =nf90_put_var(nc_ofid,time_id,times,vstart,vcount)  ! write time string

      rcall=nf90_close(nc_ofid)                                 ! close netcdf file
        endif                                                   ! end  of small
    


!------------------------------------------------------------------------------
	if(output)then

	time_index=ntimes

	allocate(stress_u(west_east,south_north,time_index))
        allocate(stress_v(west_east,south_north,time_index))
        allocate(netheat(west_east,south_north,time_index))
        allocate(longwave(west_east,south_north,time_index))
        allocate(sensible(west_east,south_north,time_index))
        allocate(latent(west_east,south_north,time_index))
        allocate(evaporation(west_east,south_north,time_index))
        allocate(precipitation(west_east,south_north,time_index))
	allocate(pressure_slp(west_east,south_north,time_index))
	allocate(rh_gen(west_east,south_north,time_index))

!---------Use TSK to fill the SST data for the land (where SST is set as -1)
    	do k=1,time_index
	do j=1,south_north
	do i=1,west_east
	if(sst(i,j,k,1).le.100)then       ! the SST units here is K
	sst(i,j,k,1)=tsk(i,j,k,1)
	endif
	enddo
	enddo 
	enddo

 

	do k=1,time_index
	do j=1,south_north
	do i=1,west_east
	
	urxx=sqrt(u10(i,j,k,1)**2+v10(i,j,k,1)**2)	! wind speed	
        taxx=t2(i,j,k,1)-273.16                         ! air temperature
        paxx=(pbtemp(i,j,k,1)+ptemp(i,j,k,1))/100.      ! pressure(mb)
	pressure_slp(i,j,k)=paxx
        tsxx=sst(i,j,k,1)-273.16                        ! sst
        zuxx=10.                                        ! wind height (10m)
	ztxx=2.                                         ! air t height (2m)
	zqxx=2.                                         ! Humidity height (2m)

	rh_valxx=t2(i,j,k,1)-273.16
	rh_valxx=6.112*exp(17.67*rh_valxx/(rh_valxx+243.5))       
        rh_valxx=0.622*rh_valxx/(paxx-rh_valxx)
        rh_valxx=q2(i,j,k,1)*100./rh_valxx              ! Relative Humidity
        if(rh_valxx.gt.100) rh_valxx=100.

	rh_gen(i,j,k)=rh_valxx

	dswxx=swdown(i,j,k,1)
	dlwxx=glw(i,j,k,1)


! -------------calculate wind stress and heat flux
	gb_lat=42.          ! for latitude of Georges Bank regions
	gb_zi=600.          ! default for PBL height

	call coare26sn(urxx,zuxx,taxx,ztxx,rh_valxx,zqxx,paxx,tsxx,dswxx,dlwxx,gb_lat,gb_zi,tauxx,hsbxx,hlbxx)


	theta=180.+atan2(u10(i,j,k,1),v10(i,j,k,1))*180./3.1415926


      IF(theta.GT.270.)THEN
      AG=3.14159/180.*(theta-270.)
      stress_u(i,j,k)=tauxx*COS(AG)
      stress_v(i,j,k)=-tauxx*SIN(AG)
      GO TO 55
      ENDIF

      IF(theta.GT.0.AND.theta.LT.90.)THEN
      AG=3.14159/180.*(90.-theta)
      stress_u(i,j,k)=-tauxx*COS(AG)
      stress_v(i,j,k)=-tauxx*SIN(AG)
      GO TO 55
      ENDIF

      IF(theta.GT.90.AND.theta.LT.180)THEN
      AG=3.14159/180.*(theta-90.)
      stress_u(i,j,k)=-tauxx*COS(AG)
      stress_v(i,j,k)=tauxx*SIN(AG)
      GO TO 55
      ENDIF

      IF(theta.EQ.90.)THEN
      stress_u(i,j,k)=-tauxx
      stress_v(i,j,k)=0.0
      GO TO 55
      ENDIF

      IF(theta.EQ.180.)THEN
      stress_u(i,j,k)=0.0
      stress_v(i,j,k)=tauxx
      GO TO 55
      ENDIF

      IF(theta.EQ.270.)THEN
      stress_u(i,j,k)=tauxx
      stress_v(i,j,k)=0.0
      GO TO 55
      ENDIF

      IF(theta.EQ.360.OR.theta.EQ.0)THEN
      stress_u(i,j,k)=0.0
      stress_v(i,j,k)=-tauxx
      GO TO 55
      ENDIF

      AG=3.14159/180.*(theta-180.)
      stress_u(i,j,k)=tauxx*SIN(AG)
      stress_v(i,j,k)=tauxx*COS(AG)
55    continue 


!---------longwave, sensible and latent heat flux-------------------------------------

       longwave(i,j,k)=glw(i,j,k,1)-0.98*5.6697*((sst(i,j,k,1)*0.01)**4)
       sensible(i,j,k)=hsbxx
       latent(i,j,k)=hlbxx



!--------------Net heat flux-----------------------------------------------------------
       netheat(i,j,k)=swdown(i,j,k,1)+longwave(i,j,k)+sensible(i,j,k)+latent(i,j,k)

!----------precipitation and evaportion (Units: m/s)
       if(k==1)then
       precipitation(i,j,k)=0.
       else
       precipitation(i,j,k)=(rainc(i,j,k,1)+rainnc(i,j,k,1)-rainc(i,j,k-1,1)-rainnc(i,j,k-1,1))*1000./3600.            ! units (m/s)
       endif

       evaporation(i,j,k)=hlbxx/ ((2.501-0.00237*(sst(i,j,k,1)-273.16))/(1.0E-9))  ! units(m/s)

	enddo
	enddo
	enddo






!-----------Create   forcing data file----------------------------------------

	output_M=trim(case)	
	rcall = nf90_create(path=output_M,cmode=nf90_clobber,ncid=nc_ofid)
	   call handle_err ("Create forcing netCDF file", rcall)

!--Define Fixed Model Dimensions
      rcall = nf90_def_dim(nc_ofid,"south_north",south_north,south_north_id)
      rcall = nf90_def_dim(nc_ofid,"west_east",west_east,west_east_id)
      rcall = nf90_def_dim(nc_ofid,"DateStrLen",19,DateStr_id)
!      rcall = nf90_def_dim(nc_ofid,"Times",73,time_id)
!      rcall = nf90_def_dim(nc_ofid,"Time",73,time_id)

!--Define Unlimited Model Dimensin
      rcall = nf90_def_dim(nc_ofid,"Time",nf90_unlimited,time_id)
!--Set up static dimensions

       force3m=(/west_east_id,south_north_id,time_id/)  ! for forcing
       force2m=(/west_east_id,south_north_id/)          ! for latitude, longitude
       force1m=(/DateStr_id,time_id/)                   ! for time



       rcall = nf90_def_var(nc_ofid,"XLAT",nf90_float,force2m,xlat_id)
       rcall =nf90_inq_varid(cdfid,"XLAT",temp_id)
       rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,xlat_id)
       rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,xlat_id)
       rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,xlat_id)

  
       rcall = nf90_def_var(nc_ofid,"XLONG",nf90_float,force2m,xlong_id)
       rcall =nf90_inq_varid(cdfid,"XLONG",temp_id)
       rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,xlong_id)
       rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,xlong_id)
       rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,xlong_id)
  
       rcall = nf90_def_var(nc_ofid,"U10",nf90_float,force3m,u10_id)
       rcall =nf90_inq_varid(cdfid,"U10",temp_id)
       rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,u10_id)
       rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,u10_id)
       rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,u10_id)
  
  
       rcall = nf90_def_var(nc_ofid,"V10",nf90_float,force3m,v10_id)
       rcall =nf90_inq_varid(cdfid,"V10",temp_id)
       rcall=nf90_copy_att(cdfid,temp_id,'description',nc_ofid,v10_id)
       rcall=nf90_copy_att(cdfid,temp_id,'units',nc_ofid,v10_id)
       rcall=nf90_copy_att(cdfid,temp_id,'coordinates',nc_ofid,v10_id)
  
       rcall = nf90_def_var(nc_ofid,"Stress_U",nf90_float,force3m,stress_u_id)
       rcall=nf90_put_att(nc_ofid,stress_u_id,"description","U Wind stress at sea surface, westward is negative")
       rcall=nf90_put_att(nc_ofid,stress_u_id,"units","Pa")
       rcall=nf90_put_att(nc_ofid,stress_u_id,"coordinates","XLONG XLAT")
  
  
       rcall = nf90_def_var(nc_ofid,"Stress_V",nf90_float,force3m,stress_v_id)
       rcall=nf90_put_att(nc_ofid,stress_v_id,"description","V Wind stress at sea surface, southward is negative")
       rcall=nf90_put_att(nc_ofid,stress_v_id,"units","Pa")
       rcall=nf90_put_att(nc_ofid,stress_v_id,"coordinates","XLONG XLAT")
  
  
       rcall = nf90_def_var(nc_ofid,"Net_Heat",nf90_float,force3m,netheat_id)
       rcall=nf90_put_att(nc_ofid,netheat_id,"description","Sum of shortwave, longwave, sensible and latent heat fluxes, ocean lose heat is negative")
       rcall=nf90_put_att(nc_ofid,netheat_id,"units","W m-2")
       rcall=nf90_put_att(nc_ofid,netheat_id,"coordinates","XLONG XLAT")
  
       rcall = nf90_def_var(nc_ofid,"Shortwave",nf90_float,force3m,shortwave_id)
       rcall=nf90_put_att(nc_ofid,shortwave_id,"description","Shortwave, upward is negative")
       rcall=nf90_put_att(nc_ofid,shortwave_id,"units","W m-2")
       rcall=nf90_put_att(nc_ofid,shortwave_id,"coordinates","XLONG XLAT")

       rcall = nf90_def_var(nc_ofid,"Longwave",nf90_float,force3m,longwave_id)
       rcall=nf90_put_att(nc_ofid,longwave_id,"description","Longwave, upward is negative")
       rcall=nf90_put_att(nc_ofid,longwave_id,"units","W m-2")
       rcall=nf90_put_att(nc_ofid,longwave_id,"coordinates","XLONG XLAT")

       rcall = nf90_def_var(nc_ofid,"Sensible",nf90_float,force3m,sensible_id)
       rcall=nf90_put_att(nc_ofid,sensible_id,"description","Sensible Heat flux, upward is negative")
       rcall=nf90_put_att(nc_ofid,sensible_id,"units","W m-2")
       rcall=nf90_put_att(nc_ofid,sensible_id,"coordinates","XLONG XLAT")

       rcall = nf90_def_var(nc_ofid,"Latent",nf90_float,force3m,latent_id)
       rcall=nf90_put_att(nc_ofid,latent_id,"description","Latent Heat flux, upward is negative")
       rcall=nf90_put_att(nc_ofid,latent_id,"units","W m-2")
       rcall=nf90_put_att(nc_ofid,latent_id,"coordinates","XLONG XLAT")


       rcall = nf90_def_var(nc_ofid,"Evaporation",nf90_float,force3m,evap_id)
       rcall=nf90_put_att(nc_ofid,evap_id,"description","Evaporation, ocean lose water is negative")
       rcall=nf90_put_att(nc_ofid,evap_id,"units","m s-1")
       rcall=nf90_put_att(nc_ofid,evap_id,"coordinates","XLONG XLAT")

       rcall = nf90_def_var(nc_ofid,"Precipitation",nf90_float,force3m,prec_id)
       rcall=nf90_put_att(nc_ofid,prec_id,"description","Precipitation, ocean lose water is negative")
       rcall=nf90_put_att(nc_ofid,prec_id,"units","m s-1")
       rcall=nf90_put_att(nc_ofid,prec_id,"coordinates","XLONG XLAT")

       rcall = nf90_def_var(nc_ofid,"SLP",nf90_float,force3m,pressure_slp_id)
       rcall=nf90_put_att(nc_ofid,pressure_slp_id,"description","Sea level pressure only for ocean")
       rcall=nf90_put_att(nc_ofid,pressure_slp_id,"units","mb")
       rcall=nf90_put_att(nc_ofid,pressure_slp_id,"coordinates","XLONG XLAT")

       rcall = nf90_def_var(nc_ofid,"RH",nf90_float,force3m,rh_gen_id)
       rcall=nf90_put_att(nc_ofid,rh_gen_id,"description","Relative Humidity (generated from Qv)")
       rcall=nf90_put_att(nc_ofid,rh_gen_id,"units","percent")
       rcall=nf90_put_att(nc_ofid,rh_gen_id,"coordinates","XLONG XLAT")



       rcall=nf90_def_var(nc_ofid,"Times",nf90_char,force1m,time_id)
       rcall=nf90_put_att(nc_ofid,time_id,"description","GMT time")

!--------Write global attributes

	if(.NOT. noglobal) then
        do k=1,nAtts
        rcall=nf90_inq_attname(cdfid,nf_global,k,attname)
        call handle_err ("Error reading Global Attributes", rcall)
        rcall=nf90_copy_att(cdfid,nf_global,trim(attname),nc_ofid,nf90_global)
        call handle_err ("Error writing Global Attributes", rcall)
        enddo
	endif
        rcall=nf90_put_att(nc_ofid,nf90_global,"History",history)
        rcall=nf90_put_att(nc_ofid,nf90_global,"Source",source)
        rcall=nf90_put_att(nc_ofid,nf90_global,"Institute",institute)


!--Exit Define Mode
        rcall=nf90_enddef(nc_ofid)


!!---------write data for forcing file---------------------------------------


      rcall =nf90_put_var(nc_ofid,xlat_id,xlat(:,:,1,1))
      rcall =nf90_put_var(nc_ofid,xlong_id,xlong(:,:,1,1))
      rcall =nf90_put_var(nc_ofid,u10_id,u10(:,:,:,1))
      call handle_err ("Error writing u10", rcall)
      rcall =nf90_put_var(nc_ofid,v10_id,v10(:,:,:,1))
      rcall =nf90_put_var(nc_ofid,stress_u_id,stress_u)
      rcall =nf90_put_var(nc_ofid,stress_v_id,stress_v)
      rcall =nf90_put_var(nc_ofid,netheat_id,netheat)
      rcall =nf90_put_var(nc_ofid,shortwave_id,swdown(:,:,:,1))
      rcall =nf90_put_var(nc_ofid,longwave_id,longwave)
      rcall =nf90_put_var(nc_ofid,sensible_id,sensible)
      rcall =nf90_put_var(nc_ofid,latent_id,latent)
      rcall =nf90_put_var(nc_ofid,evap_id,evaporation)
      rcall =nf90_put_var(nc_ofid,prec_id,precipitation)
      rcall =nf90_put_var(nc_ofid,pressure_slp_id,pressure_slp)     
      rcall =nf90_put_var(nc_ofid,rh_gen_id,rh_gen)
      vstart(1)=1
      vstart(2)=1
      vcount(1)=19
      vcount(2)=ntimes
      rcall =nf90_put_var(nc_ofid,time_id,times,vstart,vcount)  !! writing time

	rcall = nf90_close(nc_ofid)                             !! close forcing netcdf file

	endif

      rcall = nf_close(cdfid)

      write(*,*) "Successful "
      write(*,*) "================================================ "

	end program wrf2fvcom



!------------------------------------------------------------------------------

  subroutine help_info

  print*," "
  print*," wrf2fvcom -i wrf_netcdf_input_file -o fvcom_forcing_file \"
  print*," -forecast -s surface_meterological_variables_file         "
  print*," "
  print*," Current options available are:"
  print*," -help     : Print this information"
  print*," -h        : Print this information"
  print*," -i        : WRF input file (netcdf format)"
  print*," -o        : forcing output name - default is ARWout"
  print*," -s        : Extract and save necessary fields for backup"
  print*," -debug    : Print some debug information"
  print*," -forecast : flag for forecast data"
  print*," -hindcast : flag for hindcast data"
  print*," -noglobal : don't save global attributes"
  STOP

  end subroutine help_info



!------------------------------------------------------------------------------
  subroutine read_args(input_file,case,debug,small,output,save_file,history,noglobal)

  implicit none
  character (len=200)   :: input_file,save_file,history
  character (len=80)    :: case
  logical               :: debug
  logical               :: small
  logical               :: output
  logical               :: noglobal
  integer, external     :: iargc
  integer               :: numarg, i
  character (len=80)    :: dummy


! set up some defaults first
  
  input_file = "NaN"
  case       = "NaN"
  numarg     = iargc()
  i          = 1
  debug = .FALSE.
  save_file  = "NaN"
  history    = "Output from WRF2.2 (Unknown forecast/hindcast, please check global attributes)"
  noglobal   = .FALSE.

  if (numarg == 0) call help_info

  do while (i <= numarg)
    call getarg(i,dummy)

    if (dummy(1:1) == "-") then    ! We have an option, else it is the filename

      SELECTCASE (trim(dummy))
          CASE ("-help")
               call help_info
          CASE ("-h")
               call help_info
          CASE ("-i")
               i = i+1
               call getarg(i,dummy)         ! read input_file name
               input_file = dummy
          CASE ("-o")
	       output =.TRUE.
               i = i+1
               call getarg(i,dummy)         ! read case name
               case = dummy
          CASE ("-s")
               small = .TRUE. 
               i = i+1
               call getarg(i,dummy)         ! read save_file name
               save_file = dummy
          CASE ("-debug")
               debug = .TRUE.
	  CASE ("-forecast")
		history = "Output from WRF 2.2 Forecast"
	  CASE ("-hindcast")
		history = "Output from WRF 2.2 Hindcast"
          CASE ("-noglobal")
                noglobal=.TRUE.
          CASE DEFAULT
               call help_info
      END SELECT
    else
      input_file = dummy                     ! if option -i was not used for input file
    endif

      i = i+1

  enddo

  if (input_file == " ") call help_info


  end subroutine read_args

!------------------------------------------------------------------------------
      subroutine handle_err(message,nf_status)
      include "netcdf.inc"
      integer                 :: nf_status
      character (len=80)      :: message
      if (nf_status .ne. nf_noerr) then
         write(*,*)  'ERROR: ', trim(message)
         STOP
      endif
      end subroutine handle_err
!------------------------------------------------------------------------------





	subroutine coare26sn(u,zu,t,zt,rh,zq,p,ts,rs,rl,lat,zi,tau,hsb,hlb)
	implicit none

!    version coare26sn
! Scalar version of COARE3 code (Fairall et al, 2003) with cool-skin option
! retained but warm layer and surface wave options removed (hence the name
! COARE2.6). Assumes all input variables are single- valued scalars, with
! default values assigned for air pressure, radiation, PBL height and latitude
! if these data are not available. Input NaNs to indicate missing data.
! Defaults should be set to representative regional values if possible.
!
! Input:
!
!     u = relative wind speed (m/s) at height zu(m)
!     t = bulk air temperature (degC) at height zt(m)
!    rh = relative humidity (%) at height zq(m)
!     P = surface air pressure (mb) (default = 1015)
!    ts = water temperature (degC)
!    Rs = downward shortwave radiation (W/m^2) (default = 150)
!    Rl = downward longwave radiation (W/m^2) (default = 370)
!   lat = latitude (default = +45 N)
!    zi = PBL height (m) (default = 600m)
!
! Output:  A = [usr,tau,qsen,qlat,Cd,Ch,Ce,L,zet,dter,tkt] where
!
!   usr = friction velocity (m/s)
!   tau = wind stress (N/m^2)
!  qsen = sensible heat flux into ocean (W/m^2)
!  qlat = latent heat flux into ocean (W/m^2)
!    Cd = wind stress transfer coefficient at height zu
!    Ch = sensible heat transfer coefficient at height zu
!    Ce = latent heat transfer coefficient at height zu
!     L = Obukhov length scale (m)
!   zet = Monin-Obukhov stability parameter zu/L
!  dter = cool-skin temperature depression (degC)
!   tkt = cool-skin thickness (m)



! Notes: 1) u is the relative wind speed, i.e., the magnitude of the
!           difference between the wind (at zu) and ocean surface current
!           vectors.
!        2) Set jcool=0 in code if ts is true surface skin temperature,
!           otherwise ts is assumed the bulk temperature and jcool=1.
!        3) Set P=NaN to assign default value if no air pressure data
!           available.
!        4) Set Rs=NaN, Rl=NaN if no radiation data available.  This assigns
!           default values to Rs, Rl so that cool skin option can be applied.
!        5) Set lat=NaN and/or zi=NaN to assign default values if latitude
!           and/or PBL height not given.
!        6) The code to compute the heat flux caused by precipitation is
!           included if rain data is available (default is no rain).
!        7) Code updates the cool-skin temperature depression dter and thickness
!           tkt during iteration loop for consistency.
!        8) Number of iterations set to nits = 6.


! Reference:
!
!  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
!  Bulk parameterization of air sea fluxes: updates and verification for the
!  COARE algorithm, J. Climate, 16, 571-590.

! Code history:
!
! 1. 12/14/05 - coare26sn.m created based on J. Edson's 12/2004 Matlab code and
!    B. Weller's UOP 10/3/2003 Fortran code with additional input from J. Edson %    and C. Fairall. coare26sn differs from these earlier codes in the following %    ways: a) the number of computation loop iterations nits increased from 3 to %    6 for better convergence and b) the net longwave radiation flux Rln is
!    updated during the iteration loop to improve the cool-skin depression
!    temperature.
! 2. 12/21/05 - sign error in psiu_26 corrected.



	real jcool,lat,zi,us,u,qs,q,rain,beta,von,fdg,tdk,grav,rgas
	real le,cpa,cpv,rhoa,visa,al,be,cpw,rhow,visw,tcw,bigc,wetc
	real rns,rnl
	real qsat26sea,qsat26air,grv,ts,p,t,rs,rl,zu,zt,rh,zq,tau,qsen,qlat
	real DU,DT,DQ,TA,UG,DTER,UT,U10,USR,ZO10,CD10,CH10,CT10,ZOT10
	real CD,CT,CC,RIBCU,RIBU,ZETU,L10,NITS,PSIU_26,TSR,PSIT_26,QSR,TKT,CHARN
	real ZET, ZO, RR, L, ZOQ, ZOT, BF, HSB, HLB, QOUT, DELS, QCOL, ALQ, XLAMX, DQER
	real ch,ce
	integer i

!------------------------------------------------------------------------------

! set jcool=1 if Ts is bulk ocean temperature (default),
!     jcool=0 if Ts is true ocean skin temperature.

      jcool=1






! set local variables to default values if input is NaN
!     P=1015      ! pressure
!     Rs=370      !incident shortwave radiation
!     Rl=150      !incident longwave radiation
!     lat=42      !latitude
!     zi=600      !PBL height


! input variable u is assumed relative wind speed (magnitude of difference
! between wind and surface current vectors). to follow orginal Fairall code, set
! surface current speed us=0. if us data are available, construct u prior to
! using this code.
     us = 0*u

! convert rh to specific humidity
     Qs = qsat26sea(ts,P)/1000    ! surface water specific humidity (g/kg)
     Q  = qsat26air(t,P,rh)/1000  ! specific humidity of air (g/kg)
! set rain to zero
     rain = 0*u ! rain rate (mm/hr) - keep as option

     !***********  set constants **********************************************
     Beta = 1.2
     von  = 0.4
     fdg  = 1.00
     tdk  = 273.16
     grav = grv(lat)

     !***********  air constants **********************************************
     Rgas = 287.1
     Le   = (2.501-.00237*ts)*1e6
     cpa  = 1004.67
     cpv  = cpa*(1+0.84*Q)
     rhoa = P*100./(Rgas*(t+tdk)*(1+0.61*Q))
     visa = 1.326e-5*(1+6.542e-3*t+8.301e-6*t**2-4.84e-9*t**3)

     !***********  cool skin constants  ***************************************
     Al   = 2.1e-5*(ts+3.2)**0.79
     be   = 0.026
     cpw  = 4000
     rhow = 1022
     visw = 1e-6
     tcw  = 0.6
     bigc = 16*grav*cpw*(rhow*visw)**3./(tcw**2*rhoa**2)
     wetc = 0.622*Le*Qs/(Rgas*(ts+tdk)**2)


     !***********  net radiation fluxes ***************************************
     Rns = 0.945*Rs ! albedo correction
     Rnl = 0.97*(5.67e-8*(ts-0.3*jcool+tdk)**4-Rl) ! initial value



!****************  begin bulk loop ********************************************

     !***********  first guess ************************************************
     du = u-us
     dt = ts-t-.0098*zt
     dq = Qs-Q
     ta = t+tdk
     ug = 0.5
     dter  = 0.3
     ut    = sqrt(du**2+ug**2)
     u10   = ut*log(10/1e-4)/log(zu/1e-4)
     usr   = .035*u10
     zo10  = 0.011*usr**2/grav+0.11*visa/usr
     Cd10  = (von/log(10./zo10))**2
     Ch10  = 0.00115
     Ct10  = Ch10/sqrt(Cd10)
     zot10 = 10./exp(von/Ct10)
     Cd    = (von/log(zu/zo10))**2
     Ct    = von/log(zt/zot10)
     CC    = von*Ct/Cd
     Ribcu = -zu/(zi*0.004*Beta**3)
     Ribu  = -grav*zu/ta*((dt-dter*jcool)+0.61*ta*dq)/ut**2
     if (Ribu<0) then
        zetu = CC*Ribu/(1.+Ribu/Ribcu) ! unstable
     else
        zetu = CC*Ribu*(1.+27./9.*Ribu/CC) ! stable
     endif
     L10 = zu/zetu
     nits=6
     if (zetu>50) then ! stable with very thin M-O length relative to zu
        nits=1
     endif
     usr = ut*von/(log(zu/zo10)-psiu_26(zu/L10))
     tsr = -(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_26(zt/L10))
     qsr = -(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-psit_26(zq/L10))
     tkt = 0.001  ! cool skin thickness (m)
     charn = 0.011
     if (ut>10)then
        charn = 0.011+(ut-10)/(18-10)*(0.018-0.011)
     endif
     if (ut>18)then
        charn = .018
     endif









   !**************  bulk loop *****************************************


   do i=1,nits
     zet=von*grav*zu/ta*(tsr +0.61*ta*qsr)/(usr**2)
     zo=charn*usr**2/grav+0.11*visa/usr ! surface roughness
     rr=zo*usr/visa
     L=zu/zet
     zoq=min(1.15e-4,5.5e-5/rr**0.6) ! mositure roughness
     zot=zoq                        ! temperature roughness
     usr=ut*von/(log(zu/zo)-psiu_26(zu/L))
     tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot)-psit_26(zt/L))
     qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zoq)-psit_26(zq/L))
     Bf=-grav/ta*usr*(tsr+0.61*ta*qsr)



     if (Bf>0) then
        ug=Beta*(Bf*zi)**0.333
     else
        ug=0.2
     endif
     ut=sqrt(du**2+ug**2)
     hsb=-rhoa*cpa*usr*tsr ! flux out
     hlb=-rhoa*Le*usr*qsr  ! flux out
     qout=Rnl+hsb+hlb
     dels=Rns*(0.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4)))
     qcol=qout-dels
     alq=Al*qcol+be*hlb*cpw/Le
     if (alq>0) then
        xlamx=6./(1.+(bigc*alq/usr**4)**0.75)**0.333
        tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)
     else
        xlamx=6.0
        tkt=min(0.01, xlamx*visw/(sqrt(rhoa/rhow)*usr))
     endif
     dter=qcol*tkt/tcw
     dqer=wetc*dter
     Rnl=0.97*(5.67e-8*(ts-dter*jcool+tdk)**4-Rl) ! update dter



   enddo



!*****************  compute fluxes  *******************************************
     tau=rhoa*usr*usr*du/ut   ! wind stress
     hsb=rhoa*cpa*usr*tsr     ! sensible heat flux
     hlb=rhoa*Le*usr*qsr      ! latent heat flux

!********  compute transfer coeffs relative to ut @ meas. ht ******************
     Cd=tau/rhoa/ut/max(.1,du)
     Ch=-usr*tsr/ut/(dt-dter*jcool)
     Ce=-usr*qsr/(dq-dqer*jcool)/ut






















	end subroutine coare26sn








	real function psit_26(zet)
	implicit none
	real zet,dzet,x,psik,psic,f	
	if(zet>=0)then           ! stable
	dzet=min(50.,0.35*zet)
        psit_26=-((1.+0.6667*zet)**1.5+0.6667*(zet-14.28)*exp(-dzet)+8.525)
	else
	x=(1.-15.*zet)**0.5
	psik=2*log((1.+x)/2.)
	x=(1.-34.15*zet)**0.3333
	psic=1.5*log((1.+x+x**2)/3.)-sqrt(3.)*atan(1.+2.*x)/sqrt(3.)+4.*atan(1.)/sqrt(3.)
	f=zet**2/(1+zet**2)
	psit_26=(1-f)*psik+f*psic
	endif
	end function

	real function psiu_26(zet)
	implicit none
	real zet, dzet, x,psik,psic,f
	if(zet>=0) then
	dzet=min(50.,0.35*zet)
	psiu_26=-((1.+zet)+0.6667*(zet-14.28)*exp(-dzet)+8.525)
	else
	x=(1.-15.*zet)**0.25
	psik=2.*log((1.+x)/2.)+log((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.)
	x=(1.-10.15*zet)**0.3333
	psic=1.5*log((1.+x+x*x)/3)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.)
	f=zet**2/(1.+zet**2)
	psiu_26=(1-f)*psik+f*psic
	endif
	end function	

	real function bucksat(T,P)
	implicit none
	real T,P
	bucksat=6.1121*exp(17.502*T/(T+240.97))*(1.0007+3.46E-6*P)
	end function

	real function qsat26sea(T,P)
	implicit none
	real T,P,es,ex,bucksat
	ex=bucksat(T,P)
	es=0.98*ex
	qsat26sea=622*es/(P-0.378*es)
	end function

	real function qsat26air(T,P,rh)
	implicit none
	real T,P,es,em,bucksat,rh
	es=bucksat(T,P)
	em=0.01*rh*es
	qsat26air=622*em/(P-0.378*em)
	end function 

	real function grv(lat)
	implicit none
	real lat,pi,gamma,c1,c2,c3,c4,phi,x
	pi=3.1415926
        gamma=9.7803267715
        c1=0.0052790414
        c2=0.0000232718
        c3=0.0000001262
        c4=0.0000000007
        phi=lat*pi/180
        x=sin(phi)
        grv=gamma*( 1+c1*(x**2)+c2*(x**4)+c3*(x**6)+c4*(x**8) )
	end function
 

