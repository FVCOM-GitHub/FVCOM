PROGRAM wrf2fvcom
  USE MOD_WRF2GRD
  USE MOD_FORCE
  USE MOD_NCDIO
  ! SEE PARAMETER SETTINGS IN MOD_SST2GRD.F90
  IMPLICIT NONE
 
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  INTEGER :: STATUS
  CHARACTER(len=80) :: FNAME
  TYPE(TIME) :: NOW
  TYPE(NCFTIME), POINTER ::FTM

  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("WRF2FVCOM")

# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
# endif

  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ========================="
     WRITE(IPT,*) "!   BEGIN WRF 2 FVCOM"
     WRITE(IPT,*) "! ========================="
  END IF
  
  CALL GET_RUN_FILE

# if defined (PROJ)  
  USE_PROJ = HAVE_PROJ(PROJECTION_REFERENCE)
  IF (USE_PROJ) THEN
     if(dbg_set(dbg_log)) write(IPT,*) "! PROJ 4 CARTOGRAPHIC PROJECTION ENABLED!"
  ELSE
     CALL FATAL_ERROR("YOU MUST BUILD WITH THE CARTOGRAPHIC PROJECTION&
          & LIBRARY ENABLED TO USE WRF2FVCOM")
  END IF
# else
  USE_PROJ = .FALSE.
# endif


  CALL SET_TIME

  CALL GET_FVCOM_GRID


  CALL LOAD_WRF_IN


  ALLOCATE(WUSURF(0:NT),WVSURF(0:NT))

  ALLOCATE(QPREC2(0:MT),QEVAP2(0:MT))
  
  ALLOCATE(SWRAD_WATTS(0:MT), WTSURF_WATTS(0:MT))


  CALL MAKE_NC_OUTFILE

  FTM => NC_OUT%FTIME
  now = StartTime
  IINT = 1
  DO WHILE(now <= ENDTIME)
     
     IF (DBG_SET(DBG_LOG)) CALL PRINT_REAL_TIME(NOW,IPT,"TIME IS:")

     CALL UPDATE_WINDSTRESS(NOW,WUSURF,WVSURF)
     CALL UPDATE_PRECIPITATION(NOW,Qprec2,Qevap2)
     CALL UPDATE_HEAT(NOW,SWRAD_WATTS,WTSURF_WATTS)

     SELECT CASE(OUT_FILE_TYPE)
     CASE('binary')

        CALL DUMP_BINARY_FORCE(NOW)

     CASE('netcdf')
        
        FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
        CALL UPDATE_IODATA(NOW)
        CALL NC_WRITE_FILE(NC_OUT)

     CASE('both')

        
        CALL DUMP_BINARY_FORCE(NOW)
        
        FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
        CALL UPDATE_IODATA(NOW)
        CALL NC_WRITE_FILE(NC_OUT)
        

     CASE DEFAULT
        
        CALL FATAL_ERROR("INVALID OUT_FILE_TYPE: select netcdf, binary&
             & or both!")
     END SELECT



     NOW = NOW + INTERVAL
     IINT = IINT + 1

  END DO

  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN

!!$
!!$
!!$  ! CREATE THE OUTPUT FILE TREE
!!$  CALL MY_OUTFILE
!!$
!!$  CALL PRINT_FILE(NC_OUT)
!!$  CALL NC_WRITE_FILE(NC_OUT)
!!$
!!$
!!$
!!$
!!$  FTM => NC_OUT%FTIME
!!$  FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
!!$  CALL NC_WRITE_FILE(NC_OUT)
!!$  

  
END PROGRAM wrf2fvcom

