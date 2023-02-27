PROGRAM sst2grd
  USE MOD_NCDIO, only : update_iodata
  USE MOD_SST2GRD
  USE LIMS
  USE CONTROL, only : IPT
  USE MOD_UTILS
  USE MOD_INTERP
  USE MOD_INPUT

  ! SEE PARAMETER SETTINGS IN MOD_SST2GRD.F90
  IMPLICIT NONE
  
  INTEGER :: I

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  WRITE(IPT,*) "BEGIN SST INTERP:"

  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("SST2GRID")


# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
# endif

  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  CALL NAME_LIST_INITIALIZE

  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ========================="
     WRITE(IPT,*) "!   BEGIN SST INTERP"
     WRITE(IPT,*) "! ========================="
  END IF
  
  
  CALL READ_NAMELIST

  CALL CNTRL_PRMTRS

  CALL CHECK_IO_DIRS

  CALL GET_FVCOM_GRID

  CALL SET_TIME


  ! CREATE THE OUTPUT FILE TREE
  CALL MY_OUTFILE

  CALL PRINT_FILE(NC_OUT)
  CALL NC_WRITE_FILE(NC_OUT)

  IF (.not. CONSTANT) THEN

     CALL ALLOCATE_SPACE
     
     CALL CREATE_INTERP


     IDAY = 0
     DO   !365 or 366
        IDAY = IDAY + 1
        WRITE(*,*) "PROCESSING IDAY: ",IDAY
        
        CALL UPDATE_SST
        
        CALL UPDATE_IODATA(NC_OUT,NOW)
        write(*,*) "TIME IS: "// WRITE_DATETIME(NOW,6,TIMEZONE)
        PRINT*,"sst: min/max/mean", minval(sstout),maxval(sstout),sum(sstout,1)/size(sstout,1)
        
        
        FTM => NC_OUT%FTIME
        FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
        CALL NC_WRITE_FILE(NC_OUT)
        
        NOW = NOW + SST_INTERVAL
        
     ENDDO   !
     
  ELSE ! IF IT SHOULD JUST BE A CONSTANT VALUE...

     SSTout = cval

     DO I = 1,365
     
        call update_ioData(NC_OUT,NOW)
        
        FTM => NC_OUT%FTIME
        FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
        CALL NC_WRITE_FILE(NC_OUT)
        
        NOW = NOW + SST_INTERVAL
        
     END DO
  END IF

     
END PROGRAM sst2grd

