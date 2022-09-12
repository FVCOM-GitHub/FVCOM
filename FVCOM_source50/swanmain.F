# if defined (WAVE_CURRENT_INTERACTION)
!***********************************************************************
!                                                                      *
!JQI   PROGRAM SWAN
!***********************************************************************
!
!     Main program
!
!***********************************************************************
!
!JQI   USE VARS_WAVE
!   USE ALL_VARS
!JQI#  if defined (NETCDF_IO)
!JQI   USE MOD_NCDIO
!JQI#  endif      
!JQI#  if defined (MULTIPROCESSOR)
!JQI   USE MOD_PAR
!JQI#  endif      

!JQI   IMPLICIT NONE

!
!==============================================================================!
!  UGSWAN VERSION                                                               !
!==============================================================================!

!JQI   UGSWAN_VERSION     = 'UGSWAN_1.0'
!JQI#  if defined (NETCDF_IO)
!JQI   INSTITUTION       = 'School for Marine Science and Technology'
!JQI   NETCDF_TIMESTRING = 'seconds after 00:00:00'
!JQI#  endif

!==============================================================================!
!   SETUP PARALLEL ENVIRONMENT                                                 !
!==============================================================================!

!JQI   SERIAL = .TRUE. 
!JQI   PAR    = .FALSE. 
!JQI   MSR    = .TRUE.
!JQI   MYID   = 1
!JQI   NPROCS = 1
!JQI#  if defined (MULTIPROCESSOR)
!JQI   CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR)
!JQI#  endif
!==============================================================================!
!   IMPORT CASENAME FROM COMMAND LINE                                          !
!==============================================================================!

!JQI   CALL GET_CASENAME(CASENAME)

!  --- start SWAN run

!JQI   CALL SWMAIN                                                         

!JQI   END
 
!***********************************************************************
!                                                                      *
   SUBROUTINE SWMAIN_SETUP                                                   
!                                                                      *
!***********************************************************************
!
!     SWMAIN subroutine, calling SWINIT, SWREAD, SWCOMP 
!
!***********************************************************************
!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif      
   USE TIMECOMM                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OUTP_DATA                                                       
   USE M_GENARR                                                        
#  if defined (EXPLICIT)
   USE MOD_ACTION_EX
#  else
   USE MOD_PETSC
   USE MOD_ACTION_IM
#  endif   
!   USE MOD_USGRID 
   USE VARS_WAVE
#  if defined (SPHERICAL)   
   USE MOD_SPHERICAL
#  endif   
#  if defined (NETCDF_IO)
   USE MOD_NCDIO
#  endif
#  if defined (WAVE_SETUP)
   USE MOD_WAVESETUP                                                     
#  endif   
!
   IMPLICIT NONE

!   INTEGER :: IUNIT, IOSTAT, IT0, ITW, SAVITE, ILEN, INERR, IERR
   INTEGER :: IUNIT, IOSTAT,  SAVITE, ILEN, INERR, IERR
   INTEGER :: ISTAT, IF1, IL1, IDC,ISC,IP   
!   CHARACTER :: PTYPE, PNAME *8, COMPUT *4, DTTIWR*18                     
   CHARACTER :: PTYPE, PNAME *8, DTTIWR*18                     
   CHARACTER*20 :: NUMSTR, CHARS(1)                                       
   CHARACTER*80 :: MSGSTR                                                 
   LOGICAL :: LOPEN     

!   INTEGER, ALLOCATABLE :: CROSS(:,:)                                  
!   INTEGER, ALLOCATABLE :: BGRIDP(:)                                   
!   REAL   , ALLOCATABLE :: BSPECS(:,:,:,:)                             
!   REAL   , ALLOCATABLE :: AC1(:,:,:), AC2_TMP(:,:,:),AC2_TMP1(:)                     
   REAL   , ALLOCATABLE :: AC2_TMP(:,:,:)   !,AC2_TMP1(:)                     
!
!   REAL   , ALLOCATABLE :: BLKND(:), BLKNDC(:), OURQT(:)               

   REAL, ALLOCATABLE :: FTEMP(:)
   INTEGER :: I
   CHARACTER(LEN=100) :: NCFILE
!     --- initialize various data
   LEVERR=0                                                            
   MAXERR=1                                                            
   ITRACE=0                                                            
   INERR =0                                                            

   CALL SWINIT (INERR)                                                 
   IF(INERR > 0) RETURN                                              
   
#  if defined (NETCDF_IO)
  ! CALL SET_NCD_IO
#  endif      
!
   COMPUT = '    '

!  --- read and process user commands

   CALL SWREAD !(COMPUT)   

!JQI#  if defined (SPHERICAL)
!JQI   CALL ALLOC_SPHERE_VARS
!JQI#  endif     

!JQI   CALL TRIANGLE_GRID_EDGE
!JQI   CALL CELL_AREA
     
!  SET ISBCE AND ISONB CORRECTLY IN HALO CELLS/NODES
!JQI#  if defined (MULTIPROCESSOR)
!JQI   ALLOCATE(FTEMP(0:NT)) ; FTEMP = ISBCE
!JQI   IF(PAR)CALL EXCHANGE(EC,NT,1,MYID,NPROCS,FTEMP)
!JQI   ISBCE = FTEMP
!JQI   DEALLOCATE(FTEMP)
!JQI   ALLOCATE(FTEMP(0:MT)) ; FTEMP = ISONB
!JQI   IF(PAR)CALL EXCHANGE(NC,MT,1,MYID,NPROCS,FTEMP)
!JQI   ISONB = FTEMP
!JQI   DEALLOCATE(FTEMP)
!JQI#  endif

!#  if !defined (EXPLICIT)
!   CALL PETSc_SET   
!#  endif     

!
!  EXCHANGE SHAPE FACTOR INFORMATION
!
!JQI#  if defined (MULTIPROCESSOR)
!JQI   IF(PAR)CALL EXCHANGE(EC,NT,1,MYID,NPROCS,ART)
!JQI#  endif

!  --- allocate some arrays meant for computation                    
   IF(NUMOBS > 0)THEN
     IF(.NOT.ALLOCATED(CROSS)) ALLOCATE(CROSS(2,M))            
   ELSE
     IF(.NOT.ALLOCATED(CROSS)) ALLOCATE(CROSS(0,0))                
   ENDIF                                                             
   IF(.NOT.ALLOCATED(BSPECS)) ALLOCATE(BSPECS(MDC,MSC,NBSPEC,2))    
   IF(.NOT.ALLOCATED(BGRIDP)) ALLOCATE(BGRIDP(6*NBGRPT))            
!  --- do some preparations before computation                       
   CALL SWPREP ( BSPECS, BGRIDP, CROSS , SPCDIR, SPCSIG )                            

!  --- check all possible flags and if necessary change
!      if option is not correct
     if(msr)print*,"IWIND = 1 ",IWIND
   CALL ERRCHK                                                       
     if(msr)print*,"IWIND = 2 ",IWIND
!  --- initialisation of necessary grids for depth,
!      current, wind and friction

   IF(.NOT.ALLOCATED(COMPDA))ALLOCATE(COMPDA(MT,MCMVAR),STAT=ISTAT)
   IF(ISTAT /= 0)THEN                                            
     MSGSTR = 'Allocation problem: array COMPDA'                                        
     WRITE(*,*) MSGSTR
     CALL PSTOP
   END IF                                                            

   IF(.NOT. ALLOCATED(UWWIND)) ALLOCATE(UWWIND(0:NT)) 
   IF(.NOT. ALLOCATED(VWWIND)) ALLOCATE(VWWIND(0:NT)) 
   UWWIND = 0.0_SP     !ZERO
   VWWIND = 0.0_SP     !ZERO      

# if !defined (WAVE_ONLY)
   CALL SWRBC 
# else
   IF(ICEIN_ON)THEN
     CALL SWRBC_ICE
   ELSE
     CALL SWRBC 
   END IF   
# endif
   IF(.NOT.ALLOCATED(AC2)) ALLOCATE(AC2(MDC,MSC,0:MT))     
!  --- allocate AC1 in case of non-stationary situation or in case   
!      of using the S&L scheme                                       

   IF(NSTATM == 1 .AND. MXITNS > 1 .OR. PROPSC == 3)THEN        
     IF(.NOT.ALLOCATED(AC1))THEN                           
       ALLOCATE(AC1(MDC,MSC,0:MT),STAT=ISTAT)                     
     ELSE IF(SIZE(AC1) == 0)THEN                                  
       DEALLOCATE(AC1)                                             
       ALLOCATE(AC1(MDC,MSC,0:MT),STAT=ISTAT)                     
     END IF                                                         
     IF(ISTAT /= 0)THEN                                         
       MSGSTR = 'Allocation problem: array AC1'                                       
       WRITE(*,*) MSGSTR
       CALL PSTOP
     END IF                                                         
     AC1 = 0.                                                       
   ELSE                                                              
     IF(.NOT.ALLOCATED(AC1)) ALLOCATE(AC1(0,0,0))                   
   ENDIF
   IF(LEVERR > MAXERR)THEN      
                                  
     WRITE (PRINTF, 6010) LEVERR
     IF(LEVERR < 4) WRITE (PRINTF, 6011)                           
6010 FORMAT(' ** No start of computation because of error level:',I3)
6011 FORMAT(' ** To ignore this error, change [maxerr] with the',     &
            ' SET command')                                          
   ELSE
!
     IF(ITEST >= 40)THEN                                           
       IF(NSTATC == 1)THEN                                         
         WRITE (PRINTF, '(" Type of computation: dynamic")')         
       ELSE                                                          
         IF(ONED)THEN                                              
           WRITE (PRINTF, '(" Type of computation: static 1-D")')    
         ELSE                                                        
           WRITE (PRINTF, '(" Type of computation: static 2-D")')    
         ENDIF                                                       
       ENDIF                                                         
     ENDIF
!
     IF(NSTATC == 1)THEN                                           
       IT0 = 0                                                       
       IF(ICOND == 1)THEN                                          
!
!        --- compute default initial conditions
!
         CALL SWINCO (SPCDIR, SPCSIG )    
!
!        --- reset ICOND to prevent second computation of
!            initial condition
         ICOND = 0                                                   
       ENDIF
     ELSE
       IT0 = 1
     ENDIF

     if(msr) print*,"IT0=",IT0,1
#      if defined (NETCDF_IO)     
     IF(RESTART == 'cold_start')THEN
     
     ELSE IF(RESTART == 'hot_start') THEN
#      if defined (NETCDF_IO)     
!JQI       NCFILE = "./netcdf_restart/restart.nc"
!JQI       CALL NCD_READ_GRID(NCFILE)
!JQI       ALLOCATE(AC2_TMP(MDC,MSC,MGL)); AC2_TMP = 0.0
!JQI       CALL NCD_READ_RST(NCFILE,AC2_TMP,IT0)

!JQI       IF(SERIAL) AC2 = AC2_TMP
!JQI#      if defined (MULTIPROCESSOR)
!JQI       IF(PAR)THEN              
!JQI         DO IDC = 1,MDC
!JQI           DO ISC = 1,MSC
!JQI             DO I=1,M
!JQI               AC2(IDC,ISC,I) = AC2_TMP(IDC,ISC,NGID(I))
!JQI	     END DO
!JQI	     DO I=1,NHN
!JQI	       AC2(IDC,ISC,I+M) = AC2_TMP(IDC,ISC,HN_LST(I))
!JQI	     END DO  
!JQI	   END DO
!JQI         END DO	     
!JQI       END IF
!JQI#      endif       	 
!JQI       DEALLOCATE(AC2_TMP)
#      endif       
       if(msr) print*,"IT0=",IT0,2
       IT0 = IT0+1
       TIMCO = TIMCO + DTW*IT0
     ELSE
       PRINT*,'RESTAR DEFINITION NOT CORRECT'
       PRINT*,'RESTAR==',RESTART
       CALL PSTOP
     END IF
#      endif       
                            
#    if defined(WAVE_SETUP)
     IF(LSETUP > 0)CALL ALLOC_VARS_WSU   
#    endif
!    --- loop over time steps                                        
     print*,'BEFORE ITW',IT0
     
     ITW = IT0
       
!JQI     DO 500 ITW = IT0, MTC 
!JQI       IF(MSR)PRINT*,'ITW = ',ITW                                           
!
!JQI       IF(LEVERR > MAXERR)THEN                                    
!JQI         WRITE (PRINTF, 6030) LEVERR                                 
!JQI         IF(LEVERR < 4) WRITE (PRINTF, 6011)                       
!JQI6030     FORMAT(' ** No continuation of computation because ',      &
!JQI!JQI                  'of error level:',I3)                              
!JQI         EXIT                                                        
!JQI       ENDIF                                                         
!
!      --- update boundary conditions and input fields
!
!JQI       CALL SNEXTI ( BSPECS ,BGRIDP ,AC1    ,SPCSIG ,SPCDIR ,       &
!JQI	             DEPTH  ,WLEVL  ,FRIC   ,UXB    ,UYB    ,       &
!JQI		     WXI    ,WYI    )                       

!JQI       PRINT*,'AFTER SNEXTI'
!
!JQI       print*,'COMPUT=',COMPUT,ITW
!JQI       IF (COMPUT /= 'NOCO' .AND. ITW > 0) THEN                      
!JQI         SAVITE = ITEST                                              
!JQI         IF (ICOTES > ITEST) ITEST = ICOTES
!
!        --- compute action density for current time step
!
!JQI         CALL SWCOMP( AC1, ITW, CROSS )                                      
!JQI         PRINT*,'AFTER SWCOMP'
!
!        --- set ICOND=4 for stationary computation, for next
!            (stationary) COMPUTE command                            
!JQI         ICOND = 4                                                   
!
!        --- check whether computed significant wave height at       
!            boundary differs from prescribed value given in         
!            boundary command values of incident Hs                  
!
!         if(mod(itw,CDF_INT) == 0)then
!JQI         if(mod(itw,10) == 0)then
!         if(mod(itw,30) == 0)then
!JQI           IF ( BNDCHK ) THEN                                          
!JQI             CALL HSOBND (COMPDA(1,JHSIBC))         
!JQI           ENDIF              
!JQI	 end if  

!
!--NETCDF OUTPUT---------------------------------------------------------------!
!
!JQI#        if defined (NETCDF_IO)
!JQI         IF(CDF_INT /= 0 .AND. CDF_OUT)THEN
!JQI           IF(MOD(ITW,CDF_INT)==0) CALL OUT_NETCDF(SPCDIR,SPCSIG,ITW)
!JQI         END IF
!JQI#        endif

!
!--RESTART OUTPUT--------------------------------------------------------------!
!
!JQI#        if defined (NETCDF_IO)
!JQI         IF(CDF_RST /= 0)THEN
!JQI           IF(MOD(ITW,CDF_RST) == 0)THEN
!JQI             IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    RESTART FILE'
!JQI             CALL OUT_NETCDF_RST(SPCDIR,SPCSIG,ITW)
!JQI           END IF
!JQI         END IF
!JQI#        endif
	                                          
!
!JQI         ITEST = SAVITE                                              

!JQI       ENDIF
!
!JQI       IF(ITW == IT0 .AND. .NOT.ALLOCATED(OURQT))THEN             
!JQI         ALLOCATE (OURQT(MAX_OUTP_REQ))                             
!JQI         OURQT = -9999.                                             
!JQI       ENDIF                                                         
!
!JQI       SAVITE = ITEST                                                
!JQI       IF(IOUTES > ITEST) ITEST = IOUTES

!
!JQI       IF(ERRPTS > 0) REWIND(ERRPTS)                               
!JQI       ITEST = SAVITE                                                

!      --- update time

!JQI       IF(NSTATM == 1)THEN                                         
!JQI         IF(NSTATC == 1 .AND. ITW < MTC) TIMCO = TIMCO + DTW         
!JQI!JQI  !JQI           CHTIME = DTTIWR(ITMOPT, TIMCO)                              
!JQI         IF(NSTATC == 1) WRITE (PRINTF, 222) CHTIME, TIMCO          
!JQI222      FORMAT(' Time of computation ->  ',A,' in sec:', F9.0)      
!JQI       ENDIF                                                         

!JQI500  CONTINUE

!JQI     IF(LEVERR > MAXERR) GOTO 900                                  

   END IF

!JQI900 CONTINUE

!JQI   DO IUNIT=1,HIOPEN                                                   
!JQI     INQUIRE(UNIT=IUNIT,OPENED=LOPEN)                                  
!JQI     IF(LOPEN .AND. IUNIT /= PRINTF) CLOSE(IUNIT)                       
!JQI   END DO                                                              
!
!JQI   INQUIRE(UNIT=PRINTF,OPENED=LOPEN)                                   
!JQI   IF(LOPEN) CLOSE(PRINTF)                                            
!
!  --- deallocate all allocated arrays                                 

!JQI   IF(ALLOCATED(AC1   )) DEALLOCATE(AC1   )                           
!JQI   IF(ALLOCATED(BGRIDP)) DEALLOCATE(BGRIDP)                           
!JQI   IF(ALLOCATED(BSPECS)) DEALLOCATE(BSPECS)                           
!JQI   IF(ALLOCATED(COMPDA)) DEALLOCATE(COMPDA)                           
!JQI   IF(ALLOCATED(CROSS )) DEALLOCATE(CROSS )                           
!JQI   IF(ALLOCATED(OURQT )) DEALLOCATE(OURQT )                           
!JQI   IF(ALLOCATED(BLKND )) DEALLOCATE(BLKND )                           

!JQI#  if !defined (EXPLICIT)
!JQI   CALL PETSc_CLEANUP
!JQI#  endif

   RETURN                                                              
   END SUBROUTINE SWMAIN_SETUP
 
!***********************************************************************
!                                                                      *
   SUBROUTINE SWMAIN_LOOP                                                   
!                                                                      *
!***********************************************************************
!
!     SWMAIN subroutine, calling SWINIT, SWREAD, SWCOMP 
!
!***********************************************************************
!
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif      
   USE TIMECOMM                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OUTP_DATA                                                       
   USE M_GENARR                                                        
#  if defined (EXPLICIT)
   USE MOD_ACTION_EX
#  else   
   USE MOD_PETSC
   USE MOD_ACTION_IM
#  endif   
!   USE MOD_USGRID 
   USE VARS_WAVE
#  if defined (SPHERICAL)   
   USE MOD_SPHERICAL
#  endif   
#  if defined (NETCDF_IO)
   USE MOD_NCDIO
#  endif
#  if defined (WAVE_SETUP)
   USE MOD_WAVESETUP                                                     
#  endif   
!
   IMPLICIT NONE

!   INTEGER :: IUNIT, IOSTAT, IT0, ITW, SAVITE, ILEN, INERR, IERR
   INTEGER :: IUNIT, IOSTAT , SAVITE, ILEN, INERR, IERR
   INTEGER :: ISTAT, IF1, IL1, IDC,ISC,IP   
!   CHARACTER :: PTYPE, PNAME *8, COMPUT *4, DTTIWR*18                     
   CHARACTER :: PTYPE, PNAME *8, DTTIWR*18                     
   CHARACTER*20 :: NUMSTR, CHARS(1)                                       
   CHARACTER*80 :: MSGSTR                                                 
   LOGICAL :: LOPEN     

   REAL, ALLOCATABLE :: FTEMP(:)
   INTEGER :: I
   CHARACTER(LEN=100) :: NCFILE

       IF(MSR)PRINT*,'ITW = ',ITW                                           
!
!print*,LEVERR,MAXERR
       IF(LEVERR > MAXERR)THEN                                    
         WRITE (PRINTF, 6030) LEVERR                                 
         IF(LEVERR < 4) WRITE (PRINTF, 6011)                       
6030     FORMAT(' ** No continuation of computation because ',      &
                  'of error level:',I3)                              
6011 FORMAT(' ** To ignore this error, change [maxerr] with the',     &
            ' SET command')                                          
         CALL PSTOP       !EXIT                                                        
       ENDIF                                                         
!
!      --- update boundary conditions and input fields
!
!JQI       CALL SNEXTI ( BSPECS ,BGRIDP ,AC1    ,SPCSIG ,SPCDIR ,       &
!JQIJQI       CALL SNEXTI (SPCSIG ,SPCDIR ,DEPTH ,WLEVL ,FRIC ,UXB ,UYB ,WXI ,WYI)                       
       CALL SNEXTI (SPCSIG ,SPCDIR ,DEPTH ,WLEVL ,FRIC ,UXB ,UYB)                       
!
       !print*,'COMPUT=',COMPUT,ITW
       IF (COMPUT /= 'NOCO' .AND. ITW > 0) THEN                      
         SAVITE = ITEST                                              
         IF (ICOTES > ITEST) ITEST = ICOTES
!
!        --- compute action density for current time step
!
!         CALL SWCOMP( AC1, ITW, CROSS )                                      
         CALL SWCOMP(ITW)                                      
         !PRINT*,'AFTER SWCOMP'
!
!        --- set ICOND=4 for stationary computation, for next
!            (stationary) COMPUTE command                            
         ICOND = 4                                                   
!
!        --- check whether computed significant wave height at       
!            boundary differs from prescribed value given in         
!            boundary command values of incident Hs                  
!
!         if(mod(itw,CDF_INT) == 0)then
!         if(mod(itw,10) == 0)then
!         if(mod(itw,30) == 0)then
!           IF ( BNDCHK ) THEN                                          
!             CALL HSOBND (COMPDA(1,JHSIBC),DEPTH)         
!           ENDIF              
!	 end if  
         !PRINT*,'AFTER HSOBND'
!
!--NETCDF OUTPUT---------------------------------------------------------------!
!
! print*,'CDF_INT,ITW=',ITW
!#        if defined (NETCDF_IO)
         !IF(CDF_INT /= 0 .AND. CDF_OUT)THEN
         CALL SWANOUT(COMPDA(1,JDP2))
         !END IF
!#        endif

!
!--RESTART OUTPUT--------------------------------------------------------------!
!
#        if defined (NETCDF_IO)
!JQI         IF(CDF_RST /= 0)THEN
!JQI           IF(MOD(ITW,CDF_RST) == 0)THEN
!JQI             IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    RESTART FILE'
!JQI             CALL OUT_NETCDF_RST(SPCDIR,SPCSIG,ITW)
!JQI           END IF
!JQI         END IF
#        endif
	                                          
!
         ITEST = SAVITE                                              

       ENDIF
!
       IF(ITW == IT0 .AND. .NOT.ALLOCATED(OURQT))THEN             
         ALLOCATE (OURQT(MAX_OUTP_REQ))                             
         OURQT = -9999.                                             
       ENDIF                                                         
!
       SAVITE = ITEST                                                
       IF(IOUTES > ITEST) ITEST = IOUTES

!
       IF(ERRPTS > 0) REWIND(ERRPTS)                               
       ITEST = SAVITE                                                

!      --- update time

       IF(NSTATM == 1)THEN                                         
         IF(NSTATC == 1 .AND. ITW < MTC) TIMCO = TIMCO + DTW         
!JQI           CHTIME = DTTIWR(ITMOPT, TIMCO)                              
         IF(NSTATC == 1) WRITE (PRINTF, 222) CHTIME, TIMCO          
222      FORMAT(' Time of computation ->  ',A,' in sec:', F9.0)      
       ENDIF                                                         
       !PRINT*,'FINISHED SWANMAINLOOP'
   RETURN                                                              
   END SUBROUTINE SWMAIN_LOOP
 
!***********************************************************************
!                                                                      *
   SUBROUTINE SWMAIN_CLOSE                                                   
!                                                                      *
   USE OCPCOMM4                                                        
   USE VARS_WAVE
   USE MOD_PETSC

   IMPLICIT NONE

   INTEGER :: IUNIT
   LOGICAL :: LOPEN     

   DO IUNIT=1,HIOPEN                                                   
     INQUIRE(UNIT=IUNIT,OPENED=LOPEN)                                  
     IF(LOPEN .AND. IUNIT /= PRINTF) CLOSE(IUNIT)                       
   END DO                                                              
!
   INQUIRE(UNIT=PRINTF,OPENED=LOPEN)                                   
   IF(LOPEN) CLOSE(PRINTF)                                            
!
!  --- deallocate all allocated arrays                                 

   IF(ALLOCATED(AC1   )) DEALLOCATE(AC1   )                           
   IF(ALLOCATED(BGRIDP)) DEALLOCATE(BGRIDP)                           
   IF(ALLOCATED(BSPECS)) DEALLOCATE(BSPECS)                           
   IF(ALLOCATED(COMPDA)) DEALLOCATE(COMPDA)                           
   IF(ALLOCATED(CROSS )) DEALLOCATE(CROSS )                           
   IF(ALLOCATED(OURQT )) DEALLOCATE(OURQT )                           
   IF(ALLOCATED(BLKND )) DEALLOCATE(BLKND )                           

#  if !defined (EXPLICIT)
   CALL PETSc_CLEANUP
#  endif
   
   RETURN                                                              
   END SUBROUTINE SWMAIN_CLOSE
 
!***********************************************************************
!                                                                      *
   SUBROUTINE SWINIT (INERR)                                           
!                                                                      *
!***********************************************************************
!
!     Initialize several variables and arrays
!
!***********************************************************************

   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE TIMECOMM                                                        
   USE M_SNL4                                                          
   USE M_BNDSPEC                                                       
   USE ALL_VARS, ONLY: MGL  

   IMPLICIT NONE                                       

   INTEGER :: INERR
   LOGICAL :: STPNOW                                                      

   INTEGER :: IGRID,IVT,IVTYPE,MXOUTAR
   
   VERTXT = BLANK                                                      
   VERNUM = 40.51
   WRITE (VERTXT, '(F5.2)') VERNUM                                     
!
   CALL OCPINI ('swaninit', .TRUE.,INERR) 
   IF(INERR > 0) RETURN                                              
   IF(STPNOW()) RETURN                                                
!
   WRITE (PRINTF, 6010) VERTXT                                         
6010 FORMAT (/,20X,'---------------------------------------',          &
             /,20X,'                 SWAN',                            &
	     /,20X,'SIMULATION OF WAVES IN NEAR SHORE AREAS',          &
	     /,20X,'         VERSION NUMBER ', A,                      &
	     /,20X,'---------------------------------------',//)
!
!  ***** initial values for common variables *****
!  ***** names *****
   PROJID = 'SWAN'
   PROJNR = BLANK
   PROJT1 = BLANK
   PROJT2 = BLANK
   PROJT3 = BLANK
   FNEST  = BLANK
   FBCR   = BLANK
   FBCL   = BLANK
   UH     = 'm'
   UV     = 'm/s'
   UT     = 'sec'
   UL     = 'm'
   UET    = 'm3/s'
   UDI    = 'degr'
   UST    = 'm2/s2'
   UF     = 'N/m2'
   UP     = 'W/m'
   UAP    = 'W/m2'
   UDL    = 'm2/s'
!  ***** physical parameters *****
   GRAV_W   = 9.81
   WLEV   = 0.
   CASTD  = 0.               ! const. air-sea temp diff                
   PI_W     = 4.*ATAN(1.)                                                
   PI2_W    = 2.*PI_W
   UNDFLW = 1.E-15
   DNORTH = 90.                                                        
   DEGRAD = PI_W/180.
   RHO_W    = 1025.
!  power of tail in spectrum, 1: E with f, 2: E with k,
!                             3: A with f, 4: A with k
   PWTAIL(1) = 4.
   PWTAIL(2) = 2.5
   PWTAIL(3) = PWTAIL(1)+1.                                            
   PWTAIL(4) = 3.
!  ***** number of computational grid points ****                      
   MCGRD   = MGL                                                  
!   MCGRDGL = 1                                                         
   NGRBND  = 0                                                         
!   NGRBGL  = 0                                                         
!  time of computation                                                 
   TIMCO = -1.E10                                                      
   CHTIME = '    '                                                     
!  boundary conditions                                                 
   NBFILS = 0                                                          
   NBSPEC = 0                                                          
   NBGRPT = 0                                                          
!   NBGGL  = 0                                                          
   FSHAPE = 2                                                          
   DSHAPE = 2                                                          
   PSHAPE(1) = 3.3                                                     
   PSHAPE(2) = 0.1                                                     
!  ***** input grids *****
   DO IGRID = 1, NUMGRD
     XPG(IGRID)    = 0.
     YPG(IGRID)    = 0.
     ALPG(IGRID)   = 0.
     COSPG(IGRID)  = 1.                                                
     SINPG(IGRID)  = 0.
!     DXG(IGRID)    = 0.
!     DYG(IGRID)    = 0.
!     MXG(IGRID)    = 0
!     MYG(IGRID)    = 0
     LEDS(IGRID)   = 0
     STAGX(IGRID)  = 0.                                                
     STAGY(IGRID)  = 0.                                                
     EXCFLD(IGRID) = -1.E20                                            
     IFLDYN(IGRID) = 0                                                 
     IFLTIM(IGRID) = -1.E20                                            
   END DO
!  ***** computational grid *****
!   MXC    = 0
!   MYC    = 0
!   MXCGL  = 0                                                          
!   MYCGL  = 0                                                          
!   MXF    = 1                                                          
!   MXL    = 1                                                          
!   MYF    = 1                                                          
!   MYL    = 1                                                          
   MSC    = 0
   MDC    = 0
   MTC    = 1
   ICOMP  = 1
   ALPC   = 0.
   FULCIR = .TRUE.
   SPDIR1 = 0.
!  number of points needed in computational stencil:
   ICMAX  = 3
!  ***** numerical scheme *****
   NCOR   = 1
   NSTATM = -1                                                         
   NSTATC = -1                                                         
   NCOMPT = 0                                                          
!
!  initialise number of iterations stationary and nonstationary        
   MXITST = 15                                                         
   MXITNS = 1                                                          
   ITERMX = MXITST                                                     
   ICUR   = 0
   IDIF   = 0
   IINC   = 0
!
!  --- meaning IREFR:                                                  
!      IREFR = -1: limiter on Ctheta activated                         
!      IREFR =  1: No limiter on Ctheta                                
!      IREFR =  0: No refraction                                       
!
   IREFR  = 1 !1                                                          
   ITFRE  = 1
   IWIND  = 0
   IGEN   = 3                                                          
   IQUAD  = 2                                                          
   IWCAP  = 1                                                          
   ISURF  = 1
   IBOT   = 0
   ITRIAD = 0
   VARWI  = .FALSE.
   VARFR  = .FALSE.
   VARWLV = .FALSE.                                                    
   VARAST = .FALSE.       ! True means spatially variable air-sea t.d. 
#  if defined (VEGETATION)
   VARNPL = .FALSE.
#  endif
   U10    = 0.
   WDIP   = 0.
   INRHOG = 0                                                          
   DEPMIN = 0.08    !0.05
   SY0    = 3.3
   SIGMAG = 0.1
   XOFFS  = 0.
   YOFFS  = 0.
   LXOFFS = .FALSE.
   DYNDEP = .FALSE.                                                    
   NWAMN = 0
   MXOUTAR = 0
!
   FBS%NBS = -999                                                      
!
!  Set the defaults for the MDIA:                                      

   MDIA  = 6                                                           
   ALLOCATE(LAMBDA(MDIA),CNL4_1(MDIA),CNL4_2(MDIA))                   
   LAMBDA = (/0.08,0.09,0.11,0.15,0.16,0.29/)                         
   CNL4_1 = (/8.77,-13.82,10.02,-15.92,14.41,0.65/)                    
   CNL4_2 = CNL4_1                                                     
   CNL4_1 = CNL4_1 * ((2.*PI_W)**9)                                      
   CNL4_2 = CNL4_2 * ((2.*PI_W)**9)                                      
!
!  *** Initial conditions ***
   ICOND = 0                                                          
!
   BNAUT  = .FALSE.                                                    
   BNDCHK = .TRUE.                                                     
   BRESCL = .TRUE.                                                     
   ONED   = .FALSE.                                                    
   ACUPDA = .TRUE.                                                     
   HSRERR = 0.1                                                        
!
!  higher order propagation and spherical coordinates                  
!
   PROJ_METHOD = 0                                                     
   PROPSS = 2                                                          
   PROPSN = 3                                                          
   PROPSC = 1                                                          
   PROPSL = 1                                                          
   PROPFL = 0                                                          
   WAVAGE = 0.                                                         
   KSPHER = 0                                                          
   KREPTX = 0                                                          
   REARTH2 = 2.E7/PI_W    !Jianzhong Ge                                                  
   LENDEG = 2.E7/180.                                                  
!
!  *** setup flag ***                                                  
   LSETUP = 0                                                          
!
!  *** flag for setup convergence                                      
!
   CSETUP = .TRUE.                                                     
!
!  PSETUP(1) is currently unused, but can be used as setup nesting flag
!  PSETUP(2) is the user defined correction for the level of the setup
!
   PSETUP(1) = 0.0                                                     
   PSETUP(2) = 0.0                                                     
!
!  *** ACCURACY criterion ***
!
!  *** relative error in significant wave height and mean period ***
   PNUMS(1)  = 0.02                                                    
!  *** absolute error in significant wave heigth (m) ***
   PNUMS(2)  = 0.03
!  *** absolute error in mean wave period (s) ***
   PNUMS(3)  = 0.3
!  *** total number of wet gridpoints were accuracy has ***
!  *** been reached                                     ***
   PNUMS(4)  = 98.00
!
!  *** DIFFUSION schemes ***
!
!  *** Numerical diffusion over theta ***
   PNUMS(6)  = 0.5                                                     
!  *** Numerical diffusion over sigma ***
   PNUMS(7)  = 0.5                                                     
!  *** Explicit or implicit scheme in frequency space ***
!  *** default = implicit : PNUMS(8) = 1              ***
   PNUMS(8) = 1.
!  *** diffusion coefficient for explicit scheme ***
   PNUMS(9) = 0.01
!
!  *** parameters for the SIP solver                        ***
!
!  *** Required accuracy to terminate the solver            ***
!  ***                                                      ***
!  ***  || Ax-b ||  <  eps2 * || b ||                       ***
!  ***                                                      ***
!  ***  eps2 = PNUMS(12)                                    ***
!
!  *** PNUMS(13) output for the solver. Possible values:    ***
!  ***     <0  : no output                                  ***
!  ***      0  : only fatal errors will be printed          ***
!  ***      1  : additional information about the iteration ***
!  ***           is printed                                 ***
!  ***      2  : gives a maximal amount of output           ***
!  ***           concerning the iteration process           ***
!
!  *** PNUMS(14) : maximum number of iterations             ***
!
   PNUMS(12) = 1.E-4
   PNUMS(13) = 0.
   PNUMS(14) = 20.
!
!  For the setup calculation, next parameters for the solver are used:
!
!  PNUMS(23) : required accuracy to terminate the solver
!  PNUMS(24) : output for the solver (see PNUMS(13) for meanings)
!  PNUMS(25) : maximum number of iterations
!
   PNUMS(23) = 1.E-6                                                   
   PNUMS(24) = 0.                                                      
   PNUMS(25) = 1000.                                                   
!
!  Maximum growth in spectral bin
!  The value is the default in the command GEN3 KOM
!
   PNUMS(20) = 0.1                                                     
!
!  Added coefficient for use with limiter on action (Qb switch)        
!
   PNUMS(28) = 1.                                                      
!
!  *** set the values of PNUMS that are not used equal 0. ***
!
   PNUMS(5)  = 0.
!
!  The allowed global errors in the iteration procedure:               
!  PNUMS(15) for Hs and PNUMS(16) for Tm01                             
!
   PNUMS(15) = 0.02                                                    
   PNUMS(16) = 0.02                                                    
!
!  coefficient for limitation of Ctheta                                
!  default no limitation on refraction                                 
!
   PNUMS(17) = -1.                                                     
!
!  Limitation on Froude number; current velocity is reduced if greater
!  than Pnums(18)*Sqrt(grav*depth)
!
   PNUMS(18) = 0.8                                                     
!
!  *** CFL criterion for explicit scheme in frequency space ***
!
   PNUMS(19) = 0.5 * sqrt (2.)
!
!  --- coefficient for type stopping criterion                         
!
   PNUMS(21) = 0.                                                      
!
!  --- under-relaxation factor
!
   PNUMS(30) = 0.00                                                    
!
!  *** (1) and (2): Komen et al. (1984) formulation ***
!
   PWCAP(1)  = 2.36E-5
   PWCAP(2)  = 3.02E-3
   PWCAP(9)  = 2.                                                      
   PWCAP(10) = 0.                                                      
   PWCAP(11) = 1.                                                      
!
!  *** (3): Coefficient for Janssen(1989,1991) formulation ***
!  ** according to Komen et al. (1994) ***
!
   PWCAP(3)  = 4.5
   PWCAP(4)  = 0.5
!
!  *** (5): Coefficient for Longuet-Higgins ***
!
   PWCAP(5) = 1.
!
!  *** (6): ALPHA in Battjes/Janssen ***
!
   PWCAP(6) = 0.88
   PWCAP(7) = 1.
   PWCAP(8) = 0.75
!
!  *** (12): Proportionality coefficient in cumulative steepness       
!            method                                                    
!      (13): Power of cosine term in cumulative steepness method       
!
   PWCAP(12) = 0.5                                                     
   PWCAP(13) = 2.0                                                     
!
!  PBOT(1)   = 0.005          modified
   PBOT(1)   = 0.0                                                     
   PBOT(2)   = 0.015                                                   
   PBOT(3)   = 0.067
   PBOT(4)   = -0.08
   PBOT(5)   = 0.05
!
   PSURF(1)  = 1.0
   PSURF(2)  = 0.73                                                    
!
!  triad interactions
   PTRIAD(1)  = 0.1                                                    
   PTRIAD(2)  = 5.0                                                    
   PTRIAD(3)  = 10.                                                    
   PTRIAD(4)  = 0.2                                                    
   PTRIAD(5)  = 0.01                                                   
!
!  quadruplet interactions
   PQUAD(1) = 0.25                                                     
   PQUAD(2) = 3.E7                                                     
   PQUAD(3) = 5.5                                                      
   PQUAD(4) = 0.833                                                    
   PQUAD(5) = -1.25                                                    
!
   PWIND(1)  = 188.0
   PWIND(2)  = 0.59
   PWIND(3)  = 0.12
   PWIND(4)  = 250.0
   PWIND(5)  = 0.0023
   PWIND(6)  = -0.223
   PWIND(7)  = 0.
   PWIND(8)  = -0.56
   PWIND(10) = 0.0036
   PWIND(11) = 0.00123
   PWIND(12) = 1.0
   PWIND(13) = 0.13
!  *** Janssen (1991) wave growth model ***
!  *** alpha ***
   PWIND(14) = 0.01
!   PWIND(14) = 0.0144
!  *** Charnock: Von Karman constant ***
   PWIND(15) = 0.41
!  *** rho air (density) ****
   PWIND(16) = 1.28
!  *** rho water (density) ***
   PWIND(17) = RHO_W
   PWIND(9)  = PWIND(16) / RHO_W
!  Coefficient in front of A term in 3d gen. growth term
!  default is 0; can be made non-zero in command GEN3 or GROWTH
   PWIND(31) = 0.                                                      
!
!  --- coefficients for diffraction approximation                      
!
   IDIFFR    = 0                                                       
   PDIFFR(:) = 0.                                                      
!
!  pointers in array COMPDA
   JDISS = 2
   JUBOT = 3
   JQB   = 4
   JSTP  = 5
   JDHS  = 6
   JDP1  = 7
   JDP2  = 8
   JVX1  = 9
   JVY1  = 10
   JVX2  = 11
   JVY2  = 12
   JVX3  = 13
   JVY3  = 14
   JDP3  = 15
   JWX2  = 16                                                          
   JWY2  = 17                                                          
   JWX3  = 18                                                          
   JWY3  = 19                                                          
   JDTM  = 20                                                          
   JLEAK = 21                                                          
   JWLV1 = 22                                                          
   JWLV3 = 23                                                          
   JWLV2 = 24                                                          
   JHSIBC = 25                                                         
   JHS    = 26                                                         
   JURSEL = 27                                                         
   JPBOT  = 28                           
#  if defined (VEGETATION)                              
   JNPLA2 = 29
   JNPLA3 = 30
   MCMVAR = 30      !28 before vegetation
#  else
   MCMVAR = 28
#  endif                                                         
!  subarray sequence number 1 is used only for unused subarrays        
   JFRC2 = 1                                                           
   JFRC3 = 1                                                           
   JUSTAR= 1                                                           
   JZEL  = 1                                                           
   JTAUW = 1                                                           
   JCDRAG= 1                                                           
!  added for air-sea temp. diff.:
   JASTD2= 1                                                           
   JASTD3= 1                                                           
!
!  *** added for wave setup ***                                        
!
   JSETUP = 1                                                          
   JDPSAV = 1                                                          
!
!  Next pointers are for the plot of the source terms (SWTSDA)         
!
   JPWNDA = 1
   JPWNDB = 2
   JPWCAP = 3
   JPBTFR = 4
   JPWBRK = 5
   JP4S   = 6
   JP4D   = 7
   JPTRI  = 8
#  if defined (VEGETATION)
   JPVEGT = 9
   MTSVAR = 9
#  else
   MTSVAR = 8
#  endif
!
!  ***** test output control *****
   ITEST  = 1                                                          
   INTES  = 0                                                          
   ICOTES = 0                                                          
   IOUTES = 0                                                          
   LTRACE = .FALSE.
   TESTFL = .FALSE.
   NPTST  = 0
   NPTSTA = 1
   LXDMP  = -1
   LYDMP  = 0
   NEGMES = 0
   MAXMES = 200
   IFPAR = 0                                                           
   IFS1D = 0                                                           
   IFS2D = 0                                                           
!  number of obstacles initialised at 0                                
   NUMOBS = 0                                                          
!  ***** output *****
   IUBOTR = 0
!  ***** plot output *****
!
   DO IVT = 1, NMOVAR
     OVKEYW(IVT) = 'XXXX'                                              
   ENDDO
!
!  properties of output variables
!
   IVTYPE = 1
!  keyword used in SWAN command
   OVKEYW(IVTYPE) = 'XP'                                               
!  short name
   OVSNAM(IVTYPE) = 'Xp'
!  long name
   OVLNAM(IVTYPE) = 'X user coordinate'
!  unit name
   OVUNIT(IVTYPE) = UL
!  type (scalar/vector etc.)
   OVSVTY(IVTYPE) = 1
!  lower and upper limit
   OVLLIM(IVTYPE) = -1.E10
   OVULIM(IVTYPE) = 1.E10
!  lowest and highest expected value
   OVLEXP(IVTYPE) = -1.E10
   OVHEXP(IVTYPE) = 1.E10
!  exception value
   OVEXCV(IVTYPE) = -1.E10
!
   IVTYPE = 2
   OVKEYW(IVTYPE) = 'YP'                                               
   OVSNAM(IVTYPE) = 'Yp'
   OVLNAM(IVTYPE) = 'Y user coordinate'
   OVUNIT(IVTYPE) = UL
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = -1.E10
   OVULIM(IVTYPE) = 1.E10
   OVLEXP(IVTYPE) = -1.E10
   OVHEXP(IVTYPE) = 1.E10
   OVEXCV(IVTYPE) = -1.E10
!
   IVTYPE = 3
   OVKEYW(IVTYPE) = 'DIST'                                             
   OVSNAM(IVTYPE) = 'Dist'
   OVLNAM(IVTYPE) = 'distance along output curve'
   OVUNIT(IVTYPE) = UL
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.E10
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.E10
   OVEXCV(IVTYPE) = -99.
!
   IVTYPE = 4
   OVKEYW(IVTYPE) = 'DEP'                                              
   OVSNAM(IVTYPE) = 'Depth'
   OVLNAM(IVTYPE) = 'Depth'
   OVUNIT(IVTYPE) = UH
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = -1.E4
   OVULIM(IVTYPE) = 1.E4
   OVLEXP(IVTYPE) = -100.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -99.
!
   IVTYPE = 5
   OVKEYW(IVTYPE) = 'VEL'                                              
   OVSNAM(IVTYPE) = 'Vel'
   OVLNAM(IVTYPE) = 'Current velocity'
   OVUNIT(IVTYPE) = UV
   OVSVTY(IVTYPE) = 3
   OVLLIM(IVTYPE) = -100.
   OVULIM(IVTYPE) = 100.
   OVLEXP(IVTYPE) = -2.
   OVHEXP(IVTYPE) = 2.
   OVEXCV(IVTYPE) = 0.
!
   IVTYPE = 6
   OVKEYW(IVTYPE) = 'UBOT'                                             
   OVSNAM(IVTYPE) = 'Ubot'
   OVLNAM(IVTYPE) = 'Orbital velocity at the bottom'
   OVUNIT(IVTYPE) = UV
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 10.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -10.
!
   IVTYPE = 7
   OVKEYW(IVTYPE) = 'DISS'                                             
   OVSNAM(IVTYPE) = 'Dissip'
   OVLNAM(IVTYPE) = 'Energy dissipation'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 8
   OVKEYW(IVTYPE) = 'QB'                                               
   OVSNAM(IVTYPE) = 'Qb'
   OVLNAM(IVTYPE) = 'Fraction breaking waves'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -1.
!
   IVTYPE = 9
   OVKEYW(IVTYPE) = 'LEA'                                              
   OVSNAM(IVTYPE) = 'Leak'
   OVLNAM(IVTYPE) = 'Energy leak over spectral boundaries'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 10
   OVKEYW(IVTYPE) = 'HS'                                               
   OVSNAM(IVTYPE) = 'Hsig'
   OVLNAM(IVTYPE) = 'Significant wave height'
   OVUNIT(IVTYPE) = UH
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 100.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 10.
   OVEXCV(IVTYPE) = -9.
!                                                                      
   IVTYPE = 11
   OVKEYW(IVTYPE) = 'TM01'                                             
   OVSNAM(IVTYPE) = 'Tm01'                                             
   OVLNAM(IVTYPE) = 'Average absolute wave period'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 12
   OVKEYW(IVTYPE) = 'RTP'                                              
   OVSNAM(IVTYPE) = 'RTpeak'                                           
   OVLNAM(IVTYPE) = 'Relative peak period'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 13
   OVKEYW(IVTYPE) = 'DIR'                                              
   OVSNAM(IVTYPE) = 'Dir'
   OVLNAM(IVTYPE) = 'Average wave direction'
   OVUNIT(IVTYPE) = UDI                                                
   OVSVTY(IVTYPE) = 2
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 360.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 360.
   OVEXCV(IVTYPE) = -999.
!
   IVTYPE = 14
   OVKEYW(IVTYPE) = 'PDI'                                              
   OVSNAM(IVTYPE) = 'PkDir'
   OVLNAM(IVTYPE) = 'direction of the peak of the spectrum'
   OVUNIT(IVTYPE) = UDI                                                
   OVSVTY(IVTYPE) = 2
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 360.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 360.
   OVEXCV(IVTYPE) = -999.
!
   IVTYPE = 15
   OVKEYW(IVTYPE) = 'TDI'                                              
   OVSNAM(IVTYPE) = 'TDir'
   OVLNAM(IVTYPE) = 'direction of the energy transport'
   OVUNIT(IVTYPE) = UDI                                                
   OVSVTY(IVTYPE) = 2
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 360.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 360.
   OVEXCV(IVTYPE) = -999.
!
   IVTYPE = 16
   OVKEYW(IVTYPE) = 'DSPR'                                             
   OVSNAM(IVTYPE) = 'Dspr'
   OVLNAM(IVTYPE) = 'directional spreading'
   OVUNIT(IVTYPE) = UDI                                                
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 360.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 60.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 17
   OVKEYW(IVTYPE) = 'WLEN'                                             
   OVSNAM(IVTYPE) = 'Wlen'
   OVLNAM(IVTYPE) = 'Average wave length'
   OVUNIT(IVTYPE) = UL
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 200.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 18
   OVKEYW(IVTYPE) = 'STEE'                                             
   OVSNAM(IVTYPE) = 'Steepn'
   OVLNAM(IVTYPE) = 'Wave steepness'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 19
   OVKEYW(IVTYPE) = 'TRA'                                              
   OVSNAM(IVTYPE) = 'Transp'
   OVLNAM(IVTYPE) = 'Wave energy transport'
   OVUNIT(IVTYPE) = 'm3/s'                                             
   OVSVTY(IVTYPE) = 3
   OVLLIM(IVTYPE) = -100.
   OVULIM(IVTYPE) = 100.
   OVLEXP(IVTYPE) = -10.
   OVHEXP(IVTYPE) = 10.
   OVEXCV(IVTYPE) = 0.
!
   IVTYPE = 20
   OVKEYW(IVTYPE) = 'FOR'                                              
   OVSNAM(IVTYPE) = 'WForce'
   OVLNAM(IVTYPE) = 'Wave driven force per unit surface'
   OVUNIT(IVTYPE) = UF                                                 
   OVSVTY(IVTYPE) = 3
   OVLLIM(IVTYPE) = -1.E5
   OVULIM(IVTYPE) =  1.E5
   OVLEXP(IVTYPE) = -10.
   OVHEXP(IVTYPE) =  10.
   OVEXCV(IVTYPE) = 0.
!
   IVTYPE = 21
   OVKEYW(IVTYPE) = 'AAAA'                                             
   OVSNAM(IVTYPE) = 'AcDens'
   OVLNAM(IVTYPE) = 'spectral action density'
   OVUNIT(IVTYPE) = 'm2s'
   OVSVTY(IVTYPE) = 5
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -99.
!
   IVTYPE = 22
   OVKEYW(IVTYPE) = 'EEEE'                                             
   OVSNAM(IVTYPE) = 'EnDens'
   OVLNAM(IVTYPE) = 'spectral energy density'
   OVUNIT(IVTYPE) = 'm2'
   OVSVTY(IVTYPE) = 5
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -99.
!
   IVTYPE = 23
   OVKEYW(IVTYPE) = 'AAAA'                                             
   OVSNAM(IVTYPE) = 'Aux'
   OVLNAM(IVTYPE) = 'auxiliary variable'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = -1.E10
   OVULIM(IVTYPE) = 1.E10
   OVLEXP(IVTYPE) = -1.E10
   OVHEXP(IVTYPE) = 1.E10
   OVEXCV(IVTYPE) = -1.E10
!
   IVTYPE = 24
   OVKEYW(IVTYPE) = 'XC'                                               
   OVSNAM(IVTYPE) = 'Xc'
   OVLNAM(IVTYPE) = 'X computational grid coordinate'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 25
   OVKEYW(IVTYPE) = 'YC'                                               
   OVSNAM(IVTYPE) = 'Yc'
   OVLNAM(IVTYPE) = 'Y computational grid coordinate'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 26
   OVKEYW(IVTYPE) = 'WIND'                                             
   OVSNAM(IVTYPE) = 'Windv'
   OVLNAM(IVTYPE) = 'Wind velocity at 10 m above sea level'
   OVUNIT(IVTYPE) = UV                                                 
   OVSVTY(IVTYPE) = 3
   OVLLIM(IVTYPE) = -100.
   OVULIM(IVTYPE) = 100.
   OVLEXP(IVTYPE) = -50.
   OVHEXP(IVTYPE) = 50.
   OVEXCV(IVTYPE) = 0.
!
   IVTYPE = 27
   OVKEYW(IVTYPE) = 'FRC'                                              
   OVSNAM(IVTYPE) = 'FrCoef'
   OVLNAM(IVTYPE) = 'Bottom friction coefficient'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!                                                                     
   IVTYPE = 28
   OVKEYW(IVTYPE) = 'RTM01'                                           
   OVSNAM(IVTYPE) = 'RTm01'                                           
   OVLNAM(IVTYPE) = 'Average relative wave period'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 29                                                        
   OVKEYW(IVTYPE) = 'EEEE'                                             
   OVSNAM(IVTYPE) = 'EnDens'
   OVLNAM(IVTYPE) = 'energy density integrated over direction'         
   OVUNIT(IVTYPE) = 'm2'
   OVSVTY(IVTYPE) = 5
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -99.
!
   IVTYPE = 30                                                         
   OVKEYW(IVTYPE) = 'DHS'                                              
   OVSNAM(IVTYPE) = 'dHs'
   OVLNAM(IVTYPE) = 'difference in Hs between iterations'
   OVUNIT(IVTYPE) = UH
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 100.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 31                                                         
   OVKEYW(IVTYPE) = 'DRTM01'                                           
   OVSNAM(IVTYPE) = 'dTm'
   OVLNAM(IVTYPE) = 'difference in Tm between iterations'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 100.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 2.
   OVEXCV(IVTYPE) = -9.
!                                                                      
   IVTYPE = 32
   OVKEYW(IVTYPE) = 'TM02'                                             
   OVSNAM(IVTYPE) = 'Tm02'
   OVLNAM(IVTYPE) = 'Zero-crossing period'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!                                                                      
   IVTYPE = 33
   OVKEYW(IVTYPE) = 'FSPR'                                             
   OVSNAM(IVTYPE) = 'FSpr'                                             
   OVLNAM(IVTYPE) = 'Frequency spectral width (Kappa)'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 34                                                         
   OVKEYW(IVTYPE) = 'URMS'                                             
   OVSNAM(IVTYPE) = 'Urms'
   OVLNAM(IVTYPE) = 'RMS of orbital velocity at the bottom'
   OVUNIT(IVTYPE) = UV
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 10.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 35                                                         
   OVKEYW(IVTYPE) = 'UFRI'                                             
   OVSNAM(IVTYPE) = 'Ufric'
   OVLNAM(IVTYPE) = 'Friction velocity'
   OVUNIT(IVTYPE) = UV
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 10.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 36                                                         
   OVKEYW(IVTYPE) = 'ZLEN'                                             
   OVSNAM(IVTYPE) = 'Zlen'
   OVLNAM(IVTYPE) = 'Zero velocity thickness of boundary layer'
   OVUNIT(IVTYPE) = UL
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 37                                                         
   OVKEYW(IVTYPE) = 'TAUW'                                             
   OVSNAM(IVTYPE) = 'TauW'
   OVLNAM(IVTYPE) = '    '
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 10.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 38                                                         
   OVKEYW(IVTYPE) = 'CDRAG'                                            
   OVSNAM(IVTYPE) = 'Cdrag'
   OVLNAM(IVTYPE) = 'Drag coefficient'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!
!  *** wave-induced setup ***                                          
!
   IVTYPE = 39                                                         
   OVKEYW(IVTYPE) = 'SETUP'                                            
   OVSNAM(IVTYPE) = 'Setup'                                            
   OVLNAM(IVTYPE) = 'Setup due to waves'                               
   OVUNIT(IVTYPE) = 'm'                                                
   OVSVTY(IVTYPE) = 1                                                  
   OVLLIM(IVTYPE) = -1.                                                
   OVULIM(IVTYPE) = 1.                                                 
   OVLEXP(IVTYPE) = -1.                                                
   OVHEXP(IVTYPE) = 1.                                                 
   OVEXCV(IVTYPE) = -9.                                                
!
   IVTYPE = 40                                                         
   OVKEYW(IVTYPE) = 'TIME'                                             
   OVSNAM(IVTYPE) = 'Time'
   OVLNAM(IVTYPE) = 'Date-time'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -99999.
!
   IVTYPE = 41                                                         
   OVKEYW(IVTYPE) = 'TSEC'                                             
   OVSNAM(IVTYPE) = 'Tsec'
   OVLNAM(IVTYPE) = 'Time in seconds from reference time'
   OVUNIT(IVTYPE) = 's'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 100000.
   OVLEXP(IVTYPE) = -100000.
   OVHEXP(IVTYPE) = 1000000.
   OVEXCV(IVTYPE) = -99999.
!                                                        
   IVTYPE = 42
   OVKEYW(IVTYPE) = 'PER'                                              
   OVSNAM(IVTYPE) = 'Period'                                           
   OVLNAM(IVTYPE) = 'Average absolute wave period'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!                                                        
   IVTYPE = 43
   OVKEYW(IVTYPE) = 'RPER'                                             
   OVSNAM(IVTYPE) = 'RPeriod'                                          
   OVLNAM(IVTYPE) = 'Average relative wave period'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 44                                                         
   OVKEYW(IVTYPE) = 'HSWE'                                             
   OVSNAM(IVTYPE) = 'Hswell'
   OVLNAM(IVTYPE) = 'Wave height of swell part'
   OVUNIT(IVTYPE) = UH
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 100.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 10.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 45
   OVKEYW(IVTYPE) = 'URSELL'                                           
   OVSNAM(IVTYPE) = 'Ursell'
   OVLNAM(IVTYPE) = 'Ursell number'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 46
   OVKEYW(IVTYPE) = 'ASTD'                                             
   OVSNAM(IVTYPE) = 'ASTD'
   OVLNAM(IVTYPE) = 'Air-Sea temperature difference'
   OVUNIT(IVTYPE) = 'K'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = -50.
   OVULIM(IVTYPE) =  50.
   OVLEXP(IVTYPE) = -10.
   OVHEXP(IVTYPE) =  10.
   OVEXCV(IVTYPE) = -99.
!
   IVTYPE = 47                                                         
   OVKEYW(IVTYPE) = 'TMM10'
   OVSNAM(IVTYPE) = 'Tm_10'
   OVLNAM(IVTYPE) = 'Average absolute wave period'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 48                                                         
   OVKEYW(IVTYPE) = 'RTMM10'
   OVSNAM(IVTYPE) = 'RTm_10'
   OVLNAM(IVTYPE) = 'Average relative wave period'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 49
   OVKEYW(IVTYPE) = 'DIFPAR'                                           
   OVSNAM(IVTYPE) = 'DifPar'
   OVLNAM(IVTYPE) = 'Diffraction parameter'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = -50.
   OVULIM(IVTYPE) =  50.
   OVLEXP(IVTYPE) = -10.
   OVHEXP(IVTYPE) =  10.
   OVEXCV(IVTYPE) = -99.
!
   IVTYPE = 50                                                         
   OVKEYW(IVTYPE) = 'TMBOT'
   OVSNAM(IVTYPE) = 'TmBot'
   OVLNAM(IVTYPE) = 'Bottom wave period'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 51                                                         
   OVKEYW(IVTYPE) = 'WATL'
   OVSNAM(IVTYPE) = 'Watlev'
   OVLNAM(IVTYPE) = 'Water level'
   OVUNIT(IVTYPE) = UH
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = -1.E4
   OVULIM(IVTYPE) = 1.E4
   OVLEXP(IVTYPE) = -100.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -99.
!
   IVTYPE = 52                                                         
   OVKEYW(IVTYPE) = 'BOTL'
   OVSNAM(IVTYPE) = 'Botlev'
   OVLNAM(IVTYPE) = 'Bottom level'
   OVUNIT(IVTYPE) = UH
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = -1.E4
   OVULIM(IVTYPE) = 1.E4
   OVLEXP(IVTYPE) = -100.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -99.
!
   IVTYPE = 53                                                         
   OVKEYW(IVTYPE) = 'TPS'
   OVSNAM(IVTYPE) = 'TPsmoo'
   OVLNAM(IVTYPE) = 'Relative peak period (smooth)'
   OVUNIT(IVTYPE) = UT
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 54
   OVKEYW(IVTYPE) = 'DISB'                                             
   OVSNAM(IVTYPE) = 'Sfric'                                            
   OVLNAM(IVTYPE) = 'Bottom friction dissipation'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 55
   OVKEYW(IVTYPE) = 'DISSU'                                            
   OVSNAM(IVTYPE) = 'Ssurf'                                            
   OVLNAM(IVTYPE) = 'Surf breaking dissipation'                        
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 56
   OVKEYW(IVTYPE) = 'DISW'                                             
   OVSNAM(IVTYPE) = 'Swcap'                                            
   OVLNAM(IVTYPE) = 'Whitecapping dissipation'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 57
   OVKEYW(IVTYPE) = 'DISV'                                             
   OVSNAM(IVTYPE) = 'Sveg'                                             
   OVLNAM(IVTYPE) = 'Vegetation dissipation'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 58                                                         
   OVKEYW(IVTYPE) = 'QP'
   OVSNAM(IVTYPE) = 'Qp'
   OVLNAM(IVTYPE) = 'Peakedness'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 59                                                         
   OVKEYW(IVTYPE) = 'BFI'
   OVSNAM(IVTYPE) = 'BFI'
   OVLNAM(IVTYPE) = 'Benjamin-Feir index'
   OVUNIT(IVTYPE) = ' '
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 1000.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 60
   OVKEYW(IVTYPE) = 'GENE'                                             
   OVSNAM(IVTYPE) = 'Genera'
   OVLNAM(IVTYPE) = 'Energy generation'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 61
   OVKEYW(IVTYPE) = 'GENW'                                             
   OVSNAM(IVTYPE) = 'Swind'
   OVLNAM(IVTYPE) = 'Wind source term'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 62
   OVKEYW(IVTYPE) = 'REDI'                                             
   OVSNAM(IVTYPE) = 'Redist'
   OVLNAM(IVTYPE) = 'Energy redistribution'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 63
   OVKEYW(IVTYPE) = 'REDQ'                                             
   OVSNAM(IVTYPE) = 'Snl4'
   OVLNAM(IVTYPE) = 'Total absolute 4-wave interaction'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 64
   OVKEYW(IVTYPE) = 'REDT'                                             
   OVSNAM(IVTYPE) = 'Snl3'
   OVLNAM(IVTYPE) = 'Total absolute 3-wave interaction'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 65
   OVKEYW(IVTYPE) = 'PROPA'                                            
   OVSNAM(IVTYPE) = 'Propag'
   OVLNAM(IVTYPE) = 'Energy propagation'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 66
   OVKEYW(IVTYPE) = 'PROPX'                                            
   OVSNAM(IVTYPE) = 'Propxy'
   OVLNAM(IVTYPE) = 'xy-propagation'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 67
   OVKEYW(IVTYPE) = 'PROPT'                                            
   OVSNAM(IVTYPE) = 'Propth'
   OVLNAM(IVTYPE) = 'theta-propagation'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 68
   OVKEYW(IVTYPE) = 'PROPS'                                            
   OVSNAM(IVTYPE) = 'Propsi'
   OVLNAM(IVTYPE) = 'sigma-propagation'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 69
   OVKEYW(IVTYPE) = 'RADS'                                             
   OVSNAM(IVTYPE) = 'Radstr'
   OVLNAM(IVTYPE) = 'Radiation stress'
   OVUNIT(IVTYPE) = 'm2/s'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 0.1
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 70
   OVKEYW(IVTYPE) = 'NPL'                                            
   OVSNAM(IVTYPE) = 'Nplant'
   OVLNAM(IVTYPE) = 'Plants per m2'
   OVUNIT(IVTYPE) = '1/m2'
   OVSVTY(IVTYPE) = 1
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 1000.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 100.
   OVEXCV(IVTYPE) = -9.
!
   IVTYPE = 71                                                   
   OVKEYW(IVTYPE) = 'DIRBOT'                                         
   OVSNAM(IVTYPE) = 'DirBot'
   OVLNAM(IVTYPE) = 'Average bottom wave direction'
   OVUNIT(IVTYPE) = UDI                                              
   OVSVTY(IVTYPE) = 2
   OVLLIM(IVTYPE) = 0.
   OVULIM(IVTYPE) = 360.
   OVLEXP(IVTYPE) = 0.
   OVHEXP(IVTYPE) = 360.
   OVEXCV(IVTYPE) = -999.
!
!  various parameters for computation of output quantities             
!
!  reference time for TSEC
   OUTPAR(1) = 0.
!  power in expression for PER and RPER
!  previous name: SPCPOW
   OUTPAR(2) = 1.
!  power in expression for WLEN
!  previous name: AKPOWR
   OUTPAR(3) = 1.
!  indicator for direction
!  =0: direction always w.r.t. user coordinates; =1: dir w.r.t. frame
   OUTPAR(4) = 0.
!  frequency limit for swell
   OUTPAR(5) = 0.1
!
   RETURN
   END SUBROUTINE SWINIT
!
!***********************************************************************
!                                                                      *
  SUBROUTINE SWPREP( BSPECS, BGRIDP, CROSS, SPCDIR, SPCSIG )                
!                                                                      *
!***********************************************************************
!
!     do some preparations before computation is started
!
!       ----------------------------------------------------------------
!       Compute origin of computational grid in problem grid XPC, YPC
!       Compute origin of problem grid in computationl grid coordinates
!       Check input fields
!       Close files containing stationary input fields
!       Compute crossing of comp.grid lines with obstacles
!       ----------------------------------------------------------------
!***********************************************************************
!
  USE OCPCOMM4                                                        
  USE SWCOMM1                                                         
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE M_BNDSPEC                                                       
  USE M_DIFFR
  USE ALL_VARS
      
  IMPLICIT NONE                                                         

  REAL    :: SPCDIR(MDC,6)  ,    SPCSIG(MSC)                             
  INTEGER :: CROSS(2,M)
  INTEGER :: BGRIDP(6*NBGRPT)                                          
  REAL    :: BSPECS(MDC,MSC,NBSPEC,2)                                  
  INTEGER, ALLOCATABLE :: CROSSGL(:,:)  
  INTEGER :: IFLD,IP,INDX,ITMP1,ITMP2,ITMP3
  INTEGER :: II,IX,IY,IIXX,IIYY,I1,I2,NBS,INDXGR,IG
  REAL    :: ALTMP,ALBC,PPTAIL
!
  TYPE(BSDAT) , POINTER :: CURRBS                                     
  TYPE(BGPDAT), POINTER :: CURBGP                                     
!
!     coefficients for transformation from user coordinates to comp. coord.
!
  XCP  =  XPC
  YCP  =  XPC
  ALCP = -ALPC
!
! wind direction w.r.t. computational grid
!
  ALTMP = (WDIP) / PI2_W                                                
  WDIC  = PI2_W * (ALTMP - NINT(ALTMP))
!
  IF(LEDS(2) == 2)THEN
    IF(LEDS(3) /= 2) CALL MSGERR (3, 'VY not read, while VX is read')
!       ALBC  = ALPC - ALPG(2)
    ALBC  = - ALPG(2)
    COSVC = COS(ALBC)
    SINVC = SIN(ALBC)
  ENDIF
!
  IF(LEDS(4) == 2) VARFR = .TRUE.
!
  IF(LEDS(5) == 2)THEN
    IF(LEDS(6) /= 2)                                          &
      CALL MSGERR (3, 'WY not read, while WX is read')
    VARWI = .TRUE.
!       ALBC  = ALPC - ALPG(5)
    ALBC  = - ALPG(5)
    COSWC = COS(ALBC)
    SINWC = SIN(ALBC)
  ENDIF
!
  IF(LEDS(7) == 2) VARWLV = .TRUE.                                   
!
  IF(IFLDYN(1) == 1 .OR. IFLDYN(7) == 1) DYNDEP = .TRUE.             
!
!     close input files containing stationary input fields                
!
  DO IFLD = 1, NUMGRD
    IF(IFLDYN(IFLD) == 0 .AND. IFLNDS(IFLD) /= 0)THEN               
      CLOSE (IFLNDS(IFLD))
      IFLNDS(IFLD) = 0
    ENDIF
  ENDDO
!
!     computation of tail factors for moments of action spectrum
!     IP=0: action int. IP=1: energy int. IP=2: first moment of energy etc.
!
  DO IP = 0, 3
    PPTAIL = PWTAIL(1) - REAL(IP)
    PWTAIL(5+IP) = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))        
  ENDDO
!
!     *** for the four sweeps find                           ***
!     *** obstacles crossing the points in the stencil       ***
!
  IF (NUMOBS > 0)THEN
!        DO INDX = 1, M                                               
!	  CROSS(1,INDX) = 0                                              
!          CROSS(2,INDX) = 0                                              
!        ENDDO                                                            
!
!        ITMP1  = MXC                                                      
!        ITMP2  = MYC                                                      
!        ITMP3  = MCGRD                                                    
!        MXC    = MXCGL                                                    
!        MYC    = MYCGL                                                    
!        MCGRD  = MCGRDGL                                                  
!        ALLOCATE(CROSSGL(2,MCGRDGL))                                      
!        CROSSGL = 0                                                       
!JQI        CALL OBSTMOVE (XGRDGL, YGRDGL, KGRPGL)                            
!JQI        CALL SWOBST (XGRDGL, YGRDGL, KGRPGL, CROSSGL)                     
!        MXC    = ITMP1                                                    
!        MYC    = ITMP2                                                    
!        MCGRD  = ITMP3                                                    
!
!        II = 1                                                            
!!        DO IX = MXF, MXL                                                  
!           DO IG = 1, MCGRD                                               
!JQI              INDX = KGRPGL(IG)                                        
!JQI              IF(INDX /= 1)THEN                                       
!JQI                 II = II + 1                                              
!JQI                 CROSS(1:2,II) = CROSSGL(1:2,INDX)                        
!JQI              END IF                                                      
!           END DO                                                         
!!        END DO                                                            
!        DEALLOCATE(CROSSGL)                                               
!
  ENDIF
!
  CURRBS => FBS                                                       
  DO                                                                  
    NBS = CURRBS%NBS  
    IF(NBS == -999) EXIT                                            
    SPPARM(1:4) = CURRBS%SPPARM(1:4)                                 
    FSHAPE      = CURRBS%FSHAPE                                      
    DSHAPE      = CURRBS%DSHAPE                                      
    CALL SSHAPE (BSPECS(1,1,NBS,1), SPCSIG, SPCDIR,         &                 
                 FSHAPE, DSHAPE)    
    IF (.NOT.ASSOCIATED(CURRBS%NEXTBS)) EXIT                         
    CURRBS => CURRBS%NEXTBS                                          
  END DO                                                              

!JQI   BGRIDP = 0 
!JQI   IF(NBGGL == 0) NBGGL  = NBGRPT                                     
!JQI   NBGRPT = 0                                                          
!JQI   DO IG = 1,MCGRD 
!JQI     INDX = KGRPGL(IG)                                          
!JQI     CURBGP => FBGP                                                
!JQI     DO II = 1, NBGGL                                              
!JQI       INDXGR = CURBGP%BGP(1)
!JQI       IF(INDXGR == INDX)THEN                                 
!JQI         NBGRPT = NBGRPT + 1                                     
!JQI         BGRIDP(6*NBGRPT-5) = KGRPNT(IG)    !(IX-MXF+1,IY-MYF+1)          
!JQI         BGRIDP(6*NBGRPT-5) = IG           
!JQI         BGRIDP(6*NBGRPT-4) = CURBGP%BGP(2)                      
!JQI         BGRIDP(6*NBGRPT-3) = CURBGP%BGP(3)                      
!JQI         BGRIDP(6*NBGRPT-2) = CURBGP%BGP(4)                      
!JQI         BGRIDP(6*NBGRPT-1) = CURBGP%BGP(5)                      
!JQI         BGRIDP(6*NBGRPT  ) = CURBGP%BGP(6)                      
!JQI       END IF                                                     
!JQI       IF(.NOT.ASSOCIATED(CURBGP%NEXTBGP)) EXIT                  
!JQI       CURBGP => CURBGP%NEXTBGP                                   
!JQI     END DO                                                        
!JQI   END DO   
!
!     --- allocate arrays for diffraction and set prop scheme to BSBT     
!                                                                         
  IF(IDIFFR == 1)THEN                                             
    IF(.NOT.ALLOCATED(DIFPARAM)) ALLOCATE(DIFPARAM(1:M))        
    IF(.NOT.ALLOCATED(DIFPARDX)) ALLOCATE(DIFPARDX(1:M))        
    IF(.NOT.ALLOCATED(DIFPARDY)) ALLOCATE(DIFPARDY(1:M))        
    PROPSN = 1                                                       
    PROPSS = 1                                                       
  ELSE                                                                
    IF(.NOT.ALLOCATED(DIFPARAM)) ALLOCATE(DIFPARAM(0))              
    IF(.NOT.ALLOCATED(DIFPARDX)) ALLOCATE(DIFPARDX(0))              
    IF(.NOT.ALLOCATED(DIFPARDY)) ALLOCATE(DIFPARDY(0))              
  END IF                                                              
!
  RETURN
  END SUBROUTINE SWPREP
 
!**********************************************************************
!
   SUBROUTINE FLFILE (IGR1, IGR2,                                     &
                      ARR, ARR2, JX1, JX2, JX3, JY1, JY2, JY3,        &
                      COSFC, SINFC, COMPDA,                           &
                      XCGRID, YCGRID,                                 &
                      KGRPNT, IERR)
!  (This subroutine has not been used and tested yet)
!
!**********************************************************************

   USE TIMECOMM                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE M_PARALL 
!   USE MOD_UNSTRUCT_GRID                                                       

   IMPLICIT NONE
!
!  2. PURPOSE
!
!     Update boundary conditions, update nonstationary input fields
!
!  4. Argument list
!
!     ARR      real  i  array holding values read from file (x-comp)
!     ARR2     real  i  array holding values read from file (y-comp)
!     TIMR2    real i/o time of last reading of input field
!     INTRV    real  i  time interval between input fields
!     TMENDR   real  i  end time of input field
!     IGR1     int   i  location in array COMPDA for interpolated input field data (x-comp)
!     IGR2     int   i  location in array COMPDA for interpolated input field data (y-comp)
!                       for a scalar field IGR2=0
!     JX1      int   i  location in array COMPDA for interpolated input field data (x-comp)
!     JX2      int   i  location in array COMPDA for interpolated input field data (x-comp)
!     JX3      int   i  location in array COMPDA for interpolated input field data (x-comp)
!     JY1      int   i  location in array COMPDA for interpolated input field data (y-comp)
!     JY2      int   i  location in array COMPDA for interpolated input field data (y-comp)
!     JY3      int   i  location in array COMPDA for interpolated input field data (y-comp)
!     COSFC    real  i  cos of angle between input grid and computational grid
!     SINFC    real  i  sin of angle between input grid and computational grid
!     COMPDA   real i/o array holding values for computational grid points
!     XCGRID   real  i  x-coordinate of computational grid points
!     YCGRID   real  i  y-coordinate of computational grid points
!     KGRPNT   int   i  indirect addresses of computational grid points
!     NHDF     int   i  number of heading lines for a data file
!     NHDT     int   i  number of heading lines per time step
!     NHDC     int   i  number of heading lines before second component of vector field
!     IDLA     int   i  lay-out identifier for a data file
!     IDFM     int   i  format identifier for a data file
!     DFORM    char  i  format to read a data file
!     VFAC     real  i  multiplication factor applied to values from data file
!     IERR     int   o  error status: 0=no error, 9=end-of-file
!
!
   LOGICAL STPNOW                                                      
!
!  9. Structure
!
!     --------------------------------------------------------------
!     for all comp. grid points do
!         copy new values to old
!     --------------------------------------------------------------
!     repeat
!         if present time > time of last reading
!         then read new values from file
!              update time of last reading
!              interpolate values to computational grid
!         else exit from repeat
!     --------------------------------------------------------------
!     for all comp. grid points do
!         interpolate new values
!     --------------------------------------------------------------
!
! 10. SOURCE
!
!****************************************************************
!
   INTEGER    KGRPNT(MXC,MYC),                                     &
              IGR1, IGR2, JX1, JX2, JX3, JY1, JY2, JY3, IERR
!
   REAL       COMPDA(MCGRD,MCMVAR),                                &
              XCGRID(MCGRD), YCGRID(MCGRD),                            &
              COSFC, SINFC
   REAL       ARR(*), ARR2(*)
!
!     local variables
!
   INTEGER    IENT, INDX, IX, IY
!     INDX       counter of comp. grid points
!     IX         index in x-dir of comput grid point
!     IY         index in y-dir of comput grid point
!
   REAL       SVALQI
!     SVALQI     real function giving interpolated value of an input array
!
!JQI   REAL       TIMR1, XP, YP, UU, VV, VTOT, W1, W3,                &
   REAL       XP, YP, UU, VV, VTOT, W1, W3,                &
              SIZE1, SIZE2, SIZE3
   REAL(DP)       TIMR1
!     TIMR1      time of one but last input field
!     XP         x-coord of one comput grid point
!     YP         y-coord of one comput grid point
!     UU         x-component of vector, or scalar value
!     VV         y-component of vector
!     VTOT       length of vector
!     W1         weighting coeff for interpolation in time
!     W3         weighting coeff for interpolation in time
!     DIRE       direction of interpolated vector
!     SIZE1      length of vector at time TIMR1
!     SIZE2      length of vector at time TIMCO
!     SIZE3      length of vector at time TIMR2
!
   INTEGER IG
   SAVE IENT
   DATA IENT /0/
   CALL STRACE (IENT, 'FLFILE')
!
   IERR = 0
!
   IF(JX1 > 1)THEN
     DO INDX = 2, MCGRD
       COMPDA(INDX,JX1)=COMPDA(INDX,JX2)
     ENDDO
   ENDIF
   IF(IGR2 > 0 .AND. JY1 > 1)THEN
     DO INDX = 2, MCGRD
       COMPDA(INDX,JY1)=COMPDA(INDX,JY2)
     ENDDO
   ENDIF
   TIMR1 = TIMCO - DTW
!
200 IF(TIMCO <= IFLTIM(IGR1)) GOTO 400
   TIMR1 = IFLTIM(IGR1)
   IFLTIM(IGR1) = IFLTIM(IGR1) + IFLINT(IGR1)
   IF(IFLTIM(IGR1) > IFLEND(IGR1))THEN
     IFLTIM(IGR1) = 1.E10
     IF(IGR2 > 0) IFLTIM(IGR2) = IFLTIM(IGR1)
     GOTO 400
   ENDIF
   IF(IFLNDS(IGR1) > 0)THEN                                         
     IF(INODE == MASTER)THEN
       CALL INAR2D( ARR, MXG(IGR1), MYG(IGR1),                       &
                    IFLNDF(IGR1),                                    &
		    IFLNDS(IGR1), IFLIFM(IGR1), IFLFRM(IGR1),        &
		    IFLIDL(IGR1), IFLFAC(IGR1),                      &
		    IFLNHD(IGR1), IFLNHF(IGR1))
       IF(STPNOW()) RETURN                                           
     END IF
!JQI     CALL SWBROADC(IFLIDL(IGR1),1,SWINT)                              
     IF(IFLIDL(IGR1) < 0)THEN
!      end of file was encountered
       IFLTIM(IGR1) = 1.E10
       IF(IGR2 > 0) IFLTIM(IGR2) = IFLTIM(IGR1)
       GOTO 400
     ELSE
!JQI       CALL SWBROADC(ARR,MXG(IGR1)*MYG(IGR1),SWREAL)                   
     ENDIF
   ELSE                                                                
     IF(ITEST >= 20)THEN                                            
       CALL MSGERR (1,                                           &
                   'no read of input field because unit nr=0')                   
       WRITE (PRINTF, 208) IGR1
208    FORMAT (' field nr.', I2)
     ENDIF
   ENDIF                                                               
   IF(IGR2 > 0)THEN
     IFLTIM(IGR2) = IFLTIM(IGR1)
     IF(IFLNDS(IGR2) > 0)THEN                                       
       IF(INODE == MASTER)THEN                                      
         CALL INAR2D( ARR2, MXG(IGR2), MYG(IGR2),                   &
	              IFLNDF(IGR2),                                 &
		      IFLNDS(IGR2), IFLIFM(IGR2), IFLFRM(IGR2),     &
		      IFLIDL(IGR2), IFLFAC(IGR2), IFLNHD(IGR2), 0)
         IF (STPNOW()) RETURN                                         
       END IF
!JQI       CALL SWBROADC(ARR2,MXG(IGR2)*MYG(IGR2),SWREAL)                  
     ENDIF                                                             
   ENDIF
!  Interpolation over the computational grid
   DO IG = 1, MCGRD
     IF(IGR2 == 0)THEN
!JQI       COMPDA(INDX,JX3) = ARR
     ELSE
!JQI       COMPDA(INDX,JX3) = ARR
!JQI       COMPDA(INDX,JY3) = ARR2
     ENDIF
   END DO
   
   GOTO 200
!
!   Interpolation in time
!
400  W3 = (TIMCO-TIMR1) / (IFLTIM(IGR1)-TIMR1)
     W1 = 1.-W3
     IF(ITEST >= 60) WRITE(PRTEST,402) IGR1,                          &
             TIMCO,IFLTIM(IGR1),W1,W3,JX1,JY1,JX2,JY2,JX3,JY3
402  FORMAT (' input field', I2, ' interp at ', 2F9.0, 2F8.3, 6I3)

   DO INDX = 2, MCGRD
     UU = W1 * COMPDA(INDX,JX2) + W3 * COMPDA(INDX,JX3)
     IF(IGR2 <= 0)THEN
       COMPDA(INDX,JX2) = UU
     ELSE
       VV = W1 * COMPDA(INDX,JY2) + W3 * COMPDA(INDX,JY3)
       VTOT = SQRT (UU*UU + VV*VV)
!
!      procedure to prevent loss of magnitude due to interpolation
!
       IF(VTOT > 0.)THEN
         SIZE1 = SQRT(COMPDA(INDX,JX2)**2 + COMPDA(INDX,JY2)**2)
         SIZE3 = SQRT(COMPDA(INDX,JX3)**2 + COMPDA(INDX,JY3)**2)
         SIZE2 = W1*SIZE1 + W3*SIZE3
!        SIZE2 is to be length of vector
         COMPDA(INDX,JX2) = SIZE2*UU/VTOT
         COMPDA(INDX,JY2) = SIZE2*VV/VTOT
       ELSE
         COMPDA(INDX,JX2) = UU
         COMPDA(INDX,JY2) = VV
       ENDIF
     ENDIF
   END DO	

   RETURN
!
   END SUBROUTINE FLFILE

!********************************************************************
!
!**********************************************************************
!
!JQIJQI  SUBROUTINE FLFILE1 (IGR1, IGR2, ARR, ARR2, JX1, JX2, JX3, JY1,     &
  SUBROUTINE FLFILE1 (IGR1, IGR2, JX1, JX2, JX3, JY1,     &
                      JY2, JY3, COSFC, SINFC, IERR)
!
!**********************************************************************
!
!     Update boundary conditions, update nonstationary input fields
!
!**********************************************************************

  USE TIMECOMM                                                        
  USE OCPCOMM4                                                        
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
# if defined (MULTIPROCESSOR)
  USE MOD_PAR
# endif  
  USE VARS_WAVE, ONLY : M,MT,INP_WI_NTIME,COMPDA,PAR,MYID,NPROCS,    &
                        WaveTime,IMDTW,NTVE,NBVE,UWWIND,VWWIND
  USE MOD_TIME
  USE MOD_FORCE, ONLY : UPDATE_WIND2WAVE
  USE CONTROL, ONLY : WIND_ON
  USE ALL_VARS, ONLY : UUWIND,VVWIND

  IMPLICIT NONE

  LOGICAL STPNOW                                                      
!
  INTEGER :: IGR1, IGR2, JX1, JX2, JX3, JY1, JY2, JY3, IERR
  REAL    :: COSFC, SINFC
!JQIJQI  REAL    :: ARR(MT,INP_WI_NTIME), ARR2(MT,INP_WI_NTIME)
  INTEGER :: DT_WI1,DT_WI2
!
!     local variables
!
  INTEGER :: IENT, INDX, IX, IY
  REAL    :: SVALQI
!     SVALQI     real function giving interpolated value of an input array
!
!JQI  REAL    :: TIMR1, XP, YP, UU, VV, VTOT, W1, W3, SIZE1, SIZE2, SIZE3
  REAL    :: XP, YP, UU, VV, VTOT, W1, W3, SIZE1, SIZE2, SIZE3
  REAL(DP) :: TIMR1

  INTEGER :: IG,JJ

  REAL, ALLOCATABLE :: COMPDA_TMP1(:),COMPDA_TMP2(:)
!
  IERR = 0
  DT_WI1 = 0
  DT_WI2 = 0
!
  IF(JX1 > 1)THEN
!    DO INDX = 2, MCGRD
    DO INDX = 1, MT
      COMPDA(INDX,JX1)=COMPDA(INDX,JX2)
    ENDDO
  ENDIF
  IF(IGR2 > 0 .AND. JY1 > 1)THEN
!    DO INDX = 2, MCGRD
    DO INDX = 1, MT
      COMPDA(INDX,JY1)=COMPDA(INDX,JY2)
    ENDDO
  ENDIF

! TIMR1 = TIMCO - DTW
  IF(IFLTIM(IGR1) <= IFLEND(IGR1))THEN
    TIMR1 = IFLTIM(IGR1)-IFLINT(IGR1)
  ELSE
    TIMR1 = IFLEND(IGR1)
  END IF    

!   print*,"TIMCO=",TIMCO,"IFLTIM=",IFLTIM(IGR1),"IFLEND=",IFLEND(IGR1)
!   print*,"IFLINT=",IFLINT(IGR1), "DTW=",DTW

200 IF(TIMCO <= IFLTIM(IGR1)) GOTO 400
   TIMR1 = IFLTIM(IGR1)
   IFLTIM(IGR1) = IFLTIM(IGR1) + IFLINT(IGR1)
   IF(IFLTIM(IGR1) > IFLEND(IGR1))THEN
     TIMR1 = IFLEND(IGR1)
     IFLTIM(IGR1) = 1.E10
     IF(IGR2 > 0) IFLTIM(IGR2) = IFLTIM(IGR1)
     GOTO 400
   ENDIF
!JQI   DT_WI1 = INT(IFLTIM(IGR1)/IFLINT(IGR1))

   IF(IGR2 > 0)THEN
     IFLTIM(IGR2) = IFLTIM(IGR1)
!JQI     DT_WI2 = INT(IFLTIM(IGR2)/IFLINT(IGR2))
   ENDIF
!  Interpolation over the computational grid
   DO IG = 1, MT
     IF(IGR2 == 0)THEN
!JQIJQI       COMPDA(IG,JX3) = ARR(IG,DT_WI1)
     ELSE
!JQIJQI       COMPDA(IG,JX3) = ARR(IG,DT_WI1)
!JQIJQI       COMPDA(IG,JY3) = ARR2(IG,DT_WI2)
     ENDIF
   END DO
   
   GOTO 200
   
!
!   Interpolation in time
!
400  W3 = (TIMCO-TIMR1) / (IFLTIM(IGR1)-TIMR1)
     W1 = 1.-W3
     
!     WaveTime = WaveTime + IMDTW
     IF(VARWI) CALL UPDATE_WIND2WAVE(WaveTime,UWWIND,VWWIND)
     WaveTime = WaveTime + IMDTW
     
!JQI  print*,"TIMCO=",TIMCO,"IFLTIM=",IFLTIM(IGR1),"TIMR1=",TIMR1
!JQI  print*,"W3=",W3,"W1=",W1
     
  DO INDX = 1, M
    UU = 0.0_SP
    DO JJ = 1,NTVE(INDX)
      UU = UU + UWWIND(NBVE(INDX,JJ))
!JQIJQI      UU = UU + UUWIND(NBVE(INDX,JJ))
    END DO
    UU = UU/NTVE(INDX)      
!    UU = W1 * COMPDA(INDX,JX2) + W3 * COMPDA(INDX,JX3)
    IF(IGR2 <= 0)THEN
      COMPDA(INDX,JX2) = UU
    ELSE
      VV = 0.0_SP
      DO JJ = 1,NTVE(INDX)
        VV = VV + VWWIND(NBVE(INDX,JJ))
!JQIJQI        VV = VV + VVWIND(NBVE(INDX,JJ))
      END DO
      VV = VV/NTVE(INDX)	
!      VV = W1 * COMPDA(INDX,JY2) + W3 * COMPDA(INDX,JY3)
!JQI      VTOT = SQRT (UU*UU + VV*VV)
!
!     procedure to prevent loss of magnitude due to interpolation
!
!JQI      IF(VTOT > 0.)THEN
!JQI        SIZE1 = SQRT(COMPDA(INDX,JX2)**2 + COMPDA(INDX,JY2)**2)
!JQI        SIZE3 = SQRT(COMPDA(INDX,JX3)**2 + COMPDA(INDX,JY3)**2)
!JQI        SIZE2 = W1*SIZE1 + W3*SIZE3
!JQI!       SIZE2 is to be length of vector
!JQI        COMPDA(INDX,JX2) = SIZE2*UU/VTOT
!JQI        COMPDA(INDX,JY2) = SIZE2*VV/VTOT
!JQI      ELSE
        COMPDA(INDX,JX2) = UU
        COMPDA(INDX,JY2) = VV
!JQI      ENDIF
    ENDIF
  END DO	

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    ALLOCATE(COMPDA_TMP1(0:MT));   COMPDA_TMP1 = 0.0
    ALLOCATE(COMPDA_TMP2(0:MT));   COMPDA_TMP2 = 0.0
    DO IG = 1,M
      COMPDA_TMP1(IG)=COMPDA(IG,JX2)
      COMPDA_TMP2(IG)=COMPDA(IG,JY2)
    END DO  
!--------------------Jianzhong----------------------
!    CALL EXCHANGE(NC,MT,1,MYID,NPROCS,COMPDA_TMP1)
!    CALL EXCHANGE(NC,MT,1,MYID,NPROCS,COMPDA_TMP2)

    CALL AEXCHANGE(NC,MYID,NPROCS,COMPDA_TMP1)
    CALL AEXCHANGE(NC,MYID,NPROCS,COMPDA_TMP2)
!---------------------------------------------------
    DO IG = 1,MT
      COMPDA(IG,JX2) = COMPDA_TMP1(IG)
      COMPDA(IG,JY2) = COMPDA_TMP2(IG)
    END DO
    DEALLOCATE(COMPDA_TMP1,COMPDA_TMP2)
  END IF    	  
# endif

  RETURN
!
  END SUBROUTINE FLFILE1

!********************************************************************
!********************************************************************
!
# if defined (WAVE_ONLY)
!**********************************************************************
!
  SUBROUTINE FLFILE1_ICE (IGR1, IGR2, JX1, JX2, JX3, JY1,     &
                      JY2, JY3, COSFC, SINFC, IERR)
!
!**********************************************************************
!
!     Update boundary conditions, update nonstationary input fields
!
!**********************************************************************

  USE TIMECOMM                                                        
  USE OCPCOMM4                                                        
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
# if defined (MULTIPROCESSOR)
  USE MOD_PAR
# endif  
  USE VARS_WAVE, ONLY : M,MT,INP_WI_NTIME,COMPDA,PAR,MYID,NPROCS,    &
                        WaveTime,IMDTW,NTVE,NBVE,UWWIND,VWWIND
  USE MOD_TIME
  USE MOD_FORCE, ONLY : UPDATE_WIND2WAVE
  USE CONTROL, ONLY : WIND_ON
  USE ALL_VARS, ONLY : UUWIND,VVWIND,CICE_TEST
  USE MOD_NORTHPOLE
  IMPLICIT NONE

  LOGICAL STPNOW                                                      
!
  INTEGER :: IGR1, IGR2, JX1, JX2, JX3, JY1, JY2, JY3, IERR
  REAL    :: COSFC, SINFC
!JQIJQI  REAL    :: ARR(MT,INP_WI_NTIME), ARR2(MT,INP_WI_NTIME)
  INTEGER :: DT_WI1,DT_WI2
!
!     local variables
!
  INTEGER :: IENT, INDX, IX, IY
  REAL    :: SVALQI
!     SVALQI     real function giving interpolated value of an input array
!
!JQI  REAL    :: TIMR1, XP, YP, UU, VV, VTOT, W1, W3, SIZE1, SIZE2, SIZE3
  REAL    :: XP, YP, UU, VV, VTOT, W1, W3, SIZE1, SIZE2, SIZE3
!========================zhangyang=========================================
  REAL    :: DIR1, DIR2, UUVV,DIFF,CICE_TMP
!========================zhangyang=========================================
  REAL(DP) :: TIMR1

  INTEGER :: IG,JJ

  REAL, ALLOCATABLE :: COMPDA_TMP1(:),COMPDA_TMP2(:)
!
  IERR = 0
  DT_WI1 = 0
  DT_WI2 = 0
!
  IF(JX1 > 1)THEN
!    DO INDX = 2, MCGRD
    DO INDX = 1, MT
      COMPDA(INDX,JX1)=COMPDA(INDX,JX2)
    ENDDO
  ENDIF
  IF(IGR2 > 0 .AND. JY1 > 1)THEN
!    DO INDX = 2, MCGRD
    DO INDX = 1, MT
      COMPDA(INDX,JY1)=COMPDA(INDX,JY2)
    ENDDO
  ENDIF

! TIMR1 = TIMCO - DTW
  IF(IFLTIM(IGR1) <= IFLEND(IGR1))THEN
    TIMR1 = IFLTIM(IGR1)-IFLINT(IGR1)
  ELSE
    TIMR1 = IFLEND(IGR1)
  END IF    

!   print*,"TIMCO=",TIMCO,"IFLTIM=",IFLTIM(IGR1),"IFLEND=",IFLEND(IGR1)
!   print*,"IFLINT=",IFLINT(IGR1), "DTW=",DTW

200 IF(TIMCO <= IFLTIM(IGR1)) GOTO 400
   TIMR1 = IFLTIM(IGR1)
   IFLTIM(IGR1) = IFLTIM(IGR1) + IFLINT(IGR1)
   IF(IFLTIM(IGR1) > IFLEND(IGR1))THEN
     TIMR1 = IFLEND(IGR1)
     IFLTIM(IGR1) = 1.E10
     IF(IGR2 > 0) IFLTIM(IGR2) = IFLTIM(IGR1)
     GOTO 400
   ENDIF
!JQI   DT_WI1 = INT(IFLTIM(IGR1)/IFLINT(IGR1))

   IF(IGR2 > 0)THEN
     IFLTIM(IGR2) = IFLTIM(IGR1)
!JQI     DT_WI2 = INT(IFLTIM(IGR2)/IFLINT(IGR2))
   ENDIF
!  Interpolation over the computational grid
   DO IG = 1, MT
     IF(IGR2 == 0)THEN
!JQIJQI       COMPDA(IG,JX3) = ARR(IG,DT_WI1)
     ELSE
!JQIJQI       COMPDA(IG,JX3) = ARR(IG,DT_WI1)
!JQIJQI       COMPDA(IG,JY3) = ARR2(IG,DT_WI2)
     ENDIF
   END DO
   
   GOTO 200
   
!
!   Interpolation in time
!
400  W3 = (TIMCO-TIMR1) / (IFLTIM(IGR1)-TIMR1)
     W1 = 1.-W3
     
!     WaveTime = WaveTime + IMDTW
     IF(VARWI) CALL UPDATE_WIND2WAVE(WaveTime,UWWIND,VWWIND)
     WaveTime = WaveTime + IMDTW
     
!JQI  print*,"TIMCO=",TIMCO,"IFLTIM=",IFLTIM(IGR1),"TIMR1=",TIMR1
!JQI  print*,"W3=",W3,"W1=",W1
     
  DO INDX = 1, M
    UU = 0.0_SP
    VV = 0.0_SP
    DO JJ = 1,NTVE(INDX)
!    DO JJ = 1,NTVE(INDX)
!      UU = UU + UWWIND(NBVE(INDX,JJ))
!    END DO
!    UU = UU/NTVE(INDX)      
!    IF(IGR2 <= 0)THEN
!      COMPDA(INDX,JX2) = UU
!    ELSE
!      VV = 0.0_SP
!      DO JJ = 1,NTVE(INDX)
!        VV = VV + VWWIND(NBVE(INDX,JJ))
!      END DO
!      VV = VV/NTVE(INDX)	
!=========================zhangyang=====start==========================
        DIR1=ATAN2(VWWIND(NBVE(INDX,JJ)),UWWIND(NBVE(INDX,JJ)))
        DIFF=VX(INDX)-XC(NBVE(INDX,JJ))
      IF(DIFF >  180.0)THEN
      DIFF = -360.0+DIFF
      ELSE IF(DIFF < -180.0)THEN
      DIFF =  360.0+DIFF
      END IF
      DIR2=DIR1-DIFF*DEG2RAD
    IF(DIR2>(PI*2)) DIR2=DIR2-(PI*2)
    IF(DIR2<0) DIR2=DIR2+(PI*2)
    UUVV=SQRT(VWWIND(NBVE(INDX,JJ))**2+UWWIND(NBVE(INDX,JJ))**2)
    UU = UU+UUVV*COS(DIR2)
    VV = VV+UUVV*SIN(DIR2)



        !UU = UU + UWWIND(NBVE(INDX,JJ))
        !VV = VV + VWWIND(NBVE(INDX,JJ))
!JQIJQI        UU = UU + UUWIND(NBVE(INDX,JJ))
!JQIJQI        VV = VV + VVWIND(NBVE(INDX,JJ))
!IF(INDX==2435.OR.INDX==2465.OR.INDX==2439.OR.INDX==2437)PRINT *,
!"INDX",INDX,"UU",UU,"VV",VV
      END DO
IF(ICEIN_ON)THEN
      CICE_TMP=CICE_TEST(INDX)
ELSE
      CICE_TMP=0.0
END IF
      UU = UU/NTVE(INDX)*(1.-CICE_TMP)
     VV = VV/NTVE(INDX)*(1.-CICE_TMP)
!=========================zhangyang==end=============================
!=========================zhangyang===start============================
    IF(INDX==NODE_NORTHPOLE)THEN
    UU = 0.0_SP
    VV = 0.0_SP
!PRINT *, INDX    
    DO JJ = 1,NTVE(INDX)
    DIR1=ATAN2(VWWIND(NBVE(INDX,JJ)),UWWIND(NBVE(INDX,JJ)))
    DIR2=DIR1+XC(NBVE(INDX,JJ))*DEG2RAD+(PI/2)
!    DIR2=DIR1+XC(NBVE(INDX,JJ))*DEG2RAD
    IF(DIR2>(PI*2)) DIR2=DIR2-(PI*2)
    IF(DIR2<0) DIR2=DIR2+(PI*2)
    UUVV=SQRT(VWWIND(NBVE(INDX,JJ))**2+UWWIND(NBVE(INDX,JJ))**2)
    UU = UU+UUVV*COS(DIR2)
    VV = VV+UUVV*SIN(DIR2)
!PRINT *,
!"DIR1",DIR1,"DIR2",DIR2,"XC",XC(NBVE(INDX,JJ)),"XCRAD",XC(NBVE(INDX,JJ))*DEG2RAD
!PRINT *, "UWINDSPEED",UWWIND(NBVE(INDX,JJ)),"VWINDSPEED",VWWIND(NBVE(INDX,JJ))
!PRINT *, "UUVV",UUVV,"UU",UU,"VV",VV
    END DO
     UU = UU/NTVE(INDX)*(1.-CICE_TMP)
     VV = VV/NTVE(INDX)*(1.-CICE_TMP)
    END IF
!=========================zhangyang========end=======================
      



        COMPDA(INDX,JX2) = UU
        COMPDA(INDX,JY2) = VV
!JQI      ENDIF
!    ENDIF
  END DO	

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    ALLOCATE(COMPDA_TMP1(0:MT));   COMPDA_TMP1 = 0.0
    ALLOCATE(COMPDA_TMP2(0:MT));   COMPDA_TMP2 = 0.0
    DO IG = 1,M
      COMPDA_TMP1(IG)=COMPDA(IG,JX2)
      COMPDA_TMP2(IG)=COMPDA(IG,JY2)
    END DO  
!--------------------Jianzhong----------------------
!    CALL EXCHANGE(NC,MT,1,MYID,NPROCS,COMPDA_TMP1)
!    CALL EXCHANGE(NC,MT,1,MYID,NPROCS,COMPDA_TMP2)

    CALL AEXCHANGE(NC,MYID,NPROCS,COMPDA_TMP1)
    CALL AEXCHANGE(NC,MYID,NPROCS,COMPDA_TMP2)
!---------------------------------------------------
    DO IG = 1,MT
      COMPDA(IG,JX2) = COMPDA_TMP1(IG)
      COMPDA(IG,JY2) = COMPDA_TMP2(IG)
    END DO
    DEALLOCATE(COMPDA_TMP1,COMPDA_TMP2)
  END IF    	  
# endif

  RETURN
!
  END SUBROUTINE FLFILE1_ICE

!********************************************************************
# endif 
!
      SUBROUTINE ERRCHK                                                   
!
!****************************************************************
!
      USE OCPCOMM4                                                        
      USE SWCOMM1                                                         
      USE SWCOMM2                                                         
      USE SWCOMM3                                                         
      USE SWCOMM4                                                         
      USE M_GENARR  
      
      IMPLICIT NONE                                                      
!
!  0. Authors
!
!  1. Updates
!
!  2. Purpose
!
!     Check all possible combinations of physical processes if
!     they are being activated and change value of settings if
!     necessary
!
!  3. Method
!
!     0      MESSAGE
!     1      WARNING
!     2      ERROR REPAIRABLE
!     3      SEVERE ERROR (calculation continues, however problems
!                          may arise)
!     4      TERMINATION ERROR (calculation is terminated )
!
!  4. Argument variables (updated 30.72)
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     MSGSTR:     string to pass message to call MSGERR
!
      CHARACTER*20 NUMSTR, CHARS(3)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     MSGERR : Handles error messeges according to severity
!     NUMSTR : Converts integer/real to string
!     TXPBLA : Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWANCOM
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ------------------------------------------------------------
!     End of the subroutine ERRCHK
!     ------------------------------------------------------------
!
! 13. Source text
!
      LOGICAL  STPNOW                                                     
      LOGICAL  EQREAL  
      INTEGER :: GAMMA,II                                                   
      INTEGER, SAVE :: IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'ERRCHK')

!     --- choose scheme for stationary or nonstationary computation

!      PRINT*,'NSTATC = ',NSTATC
      IF (NSTATC.EQ.1) THEN                                               
        PROPSC = PROPSN                                                   
      ELSE                                                                
        PROPSC = PROPSS                                                   
      ENDIF                                                               
!
!     -----------------------------------------------------------------
!
!     *** WARNINGS AND ERROR MESSAGES ***
!
!     -----------------------------------------------------------------
!
!     *** WAM cycle 3 physics ***
!
      IF((IWIND == 3 .OR. IWIND == 5) .AND. IWCAP /= 1             &
         .AND. IWCAP /= 7)THEN
        CALL MSGERR(1,'Activate whitecapping mechanism according to')
        CALL MSGERR(1,'Komen et al. (1984) for wind option G3/YAN')       
      ENDIF

      IF(IWCAP == 7 .AND. IWIND /= 5)THEN                            
        CALL MSGERR(1,'Activate wind option Yan (1987) in case of')       
        CALL MSGERR(1,'Alves & Banner (2003) white-capping method')       
      END IF                                                              
!
!     *** WAM cycle 4 physics ***
!
      IF(IWIND == 4 .AND. IWCAP /= 2)THEN
        CALL MSGERR(1,'Activate whitecapping mechanism according to')
        CALL MSGERR(1,'Janssen (1991) for wind option JANS        ')
      ENDIF
!
!     *** check numerical scheme in presence of a current ***
!
      IF(ICUR == 1)THEN
        IF(PNUMS(6) == 0.)THEN
          CALL MSGERR(1,'In presence of a current it is recommended to')
          CALL MSGERR(1,'use an implicit upwind scheme in theta space ')
          CALL MSGERR(1,'-> set CDD = 1.')
          WRITE(PRINTF,*)
        ENDIF
        IF(PNUMS(7) == 0.)THEN
          CALL MSGERR(1,'In presence of a current it is recommended to')
          CALL MSGERR(1,'use an implicit upwind scheme in sigma space ')
          CALL MSGERR(1,'-> set CSS = 1.')
          WRITE(PRINTF,*)
        ENDIF
      END IF
!
!     check combination of REPeating option and grid type and dimension   
!
!JQI      IF(KREPTX > 0)THEN                                               
!       --- "GE" changed to "EQ" since OPTG.EQ.3 FOR CURVILINEAR          
!JQI        IF(OPTG == 3)                                             &
!JQI	  CALL MSGERR (3, 'Curvilinear grid cannot be REPeating')
!JQI        IF(PROPSC == 1 .AND. MXC < 1)                             &
!JQI	  CALL MSGERR (3, 'MXC must be >=1 for REPeating option')
!JQI        IF(PROPSC == 2 .AND. MXC < 2)                             &
!JQI	  CALL MSGERR (3, 'MXC must be >=2 for REPeating option')           
!JQI        IF(PROPSC == 3 .AND. MXC < 3)                             &
!JQI	  CALL MSGERR (3, 'MXC must be >=3 for REPeating option')           
!JQI      ENDIF
!
#     if defined (SPHERICAL)
      IF(OPTG /= 1 .AND. KSPHER /= 0)THEN                               
         CALL MSGERR (3, 'Spherical coordinates must be given in')        
         CALL MSGERR (3, 'uniform, recti-linear computational grid')      
      END IF                                                              
#     endif
!
      IF(PROPSC == 2 .AND. NSTATC > 0)THEN                             
        CALL MSGERR (3, 'SORDUP scheme only in stationary run')           
      ENDIF
      IF(PROPSC == 3 .AND. NSTATC == 0)THEN                             
        CALL MSGERR (3, 'S&L scheme not in stationary run')               
      ENDIF
!
!     --- A warning about curvilinear and S&L scheme                      
      IF((PROPSC == 3).AND.(OPTG == 3))THEN                             
         CALL MSGERR(1,'the S&L scheme (higher order nonstationary')      
         CALL MSGERR(1,'IS NOT fully implemented for curvilinear')        
         CALL MSGERR(1,'coordinates. This may or may not be noticeable')  
         CALL MSGERR(1,'in simulations. DX and DY are approximated')      
         CALL MSGERR(1,'with DX and DY of two nearest cells.')            
         CALL MSGERR(1,'Note that SORDUP (higher order stationary)')      
         CALL MSGERR(1,'IS fully implemented for curvilinear coord.')     
         CALL MSGERR(1,'and differences are usually negligible.')         
      ENDIF                                                               
!
!     Here the various problems with quadruplets are checked
!     The combination of quadruplets and sectors is an error
!     in the calculation of quadruplets when the SECTOR option is
!     used in the CGRID command. This error should be corrected
!     in the future
!
      IF(IWIND == 3 .OR. IWIND == 4)THEN                                
        IF(IQUAD == 0)THEN
          CALL MSGERR(2,'Quadruplets should be activated when SWAN  ')    
          CALL MSGERR(2,'is running in a third generation mode and  ')    
          CALL MSGERR(2,'wind is present                            ')    
        ENDIF
      ENDIF
!
      IF(IQUAD >= 1)THEN                                              
!
       IF(.NOT. FULCIR)THEN                                             
        IF((SPDIR2-SPDIR1) < (PI_W/12.))THEN                           
          CALL MSGERR(2,'A combination of using quadruplets with a'    )  
          CALL MSGERR(2,'sector of less than 30 degrees should be'     )  
          CALL MSGERR(2,'avoided at all times, it is likely to produce')  
          CALL MSGERR(2,'unreliable results and unexpected errors.'    )  
          CALL MSGERR(2,'Refer to the manual (CGRID) for details'      )  
        ELSE                                                              
          CALL MSGERR(1,'It is not recommended to use quadruplets'     )  
          CALL MSGERR(1,'in combination with calculations on a sector.')  
          CALL MSGERR(1,'Refer to the manual (CGRID) for details'      )  
        END IF                                                            
       END IF                                                             
!
       IF(IWIND == 0)THEN                                               
         CALL MSGERR(2,'It is not recommended to use quadruplets'     )   
         CALL MSGERR(2,'in combination with zero wind conditions.'    )   
       END IF                                                             
!
       IF(MSC == 3)THEN
         CALL MSGERR(4,'Do not activate quadruplets for boundary ')
         CALL MSGERR(4,'option BIN -> use other option           ')
         WRITE(PRINTF,*)
       END IF
!
      END IF                                                              
!
!     Check whether limiter should be de-activated                        
!
      IF(PNUMS(20) < 100.)THEN                                       
         IF(IGEN == 3)THEN
!        --- third generation mode
            IF(IQUAD == 0 .AND. ITRIAD == 0)THEN
               CALL MSGERR(1,'Limiter is de-activated in case of')
               CALL MSGERR(1,'no wave-wave interactions')
               PNUMS(20) = 1.E+20
            END IF
         ELSE IF(IGEN < 3)THEN
!        --- first or second generation mode
            IF(ITRIAD == 0)THEN
               CALL MSGERR(1,'Limiter is de-activated in case of')
               CALL MSGERR(1,'first or second generation mode')
               PNUMS(20) = 1.E+20
            END IF
         END IF
      END IF
!
!     Check resolution in frequency-space when DIA is used                
!
      IF(IQUAD > 0 .AND. IQUAD <= 3 .OR. IQUAD == 8)THEN               
         GAMMA = EXP(ALOG(SHIG/SLOW)/REAL(MSC-1))                         
         IF(ABS(GAMMA-1.1) > 0.055)THEN                                
            CALL MSGERR(1,                                        &
	        'relative frequency resolution (df/f) deviates more')    
            CALL MSGERR(1,                                        &
	        'than 5% from 10%-resolution. This may be problematic')  
            CALL MSGERR(1,                                        &
	        'when quadruplets are approximated by means of DIA.')    
         END IF                                                           
      END IF                                                              
!
!     When using triads MSC must be less than 200!                        
!
      IF((ITRIAD > 0).AND.(MSC > 200))THEN                            
         CALL MSGERR(4,'When triads are active the number of     ')       
         CALL MSGERR(4,'directions chosen in the CGRID command   ')       
         CALL MSGERR(4,'must be less than 200                    ')       
      END IF                                                              
!
!     When CSM is applied only DIA should be used                         
!
      IF(IWCAP == 6 .AND. IQUAD > 3 .AND. IQUAD /= 8)THEN              
         CALL MSGERR(3,                                            &
	     'In case of cumulative steepness method, only DIA '// &
	     'technique is supported')                                    
      END IF                                                              
!
!     When Alves and Banner and XNL are applied change the parameters     
!
      IF(IWCAP == 7 .AND. (IQUAD == 51 .OR. IQUAD == 52 .OR.        &
         IQUAD == 53))THEN                          
         IF(EQREAL(PWCAP( 1), 5.0E-5)) PWCAP( 1) = 5.0E-5                
         IF(EQREAL(PWCAP(12),1.75E-3)) PWCAP(12) = 1.95E-3               
         IF(EQREAL(PWCAP(10),    4.0)) PWCAP(10) = 4.                    
         IF(EQREAL(PWCAP( 9),    0.0)) PWCAP( 9) = 0.                    
         IF(EQREAL(PWCAP(11),    0.0)) PWCAP(11) = 0.                    
      END IF                                                              
!OMP!
!OMP      IF ( IQUAD.EQ.51 .OR. IQUAD.EQ.52 .OR. IQUAD.EQ.53 ) THEN           
!OMP         CALL MSGERR(4,'XNL not supported within OpenMP environment')     
!OMP      END IF                                                              
!
      IF(ITEST >= 120)THEN
        WRITE(PRINTF,3000) IWIND ,IQUAD, ICUR, IWCAP, MSC
3000    FORMAT(' ERRCHK : IWIND QUAD CUR WCAP MSC   : ',5I4)
        IF(IWIND > 0)THEN                                            
          DO II = 1, MWIND
            WRITE(PRINTF,30) II,PWIND(II)
 30         FORMAT(' PWIND(',I2,') = ',E11.4)
          ENDDO
        ENDIF
      END IF
!
      RETURN
      END SUBROUTINE ERRCHK

!
!***********************************************************************
!                                                                      *
  SUBROUTINE SWINCO (SPCDIR ,SPCSIG )                            
!                                                                      *
!***********************************************************************
!
!     Imposing of wave initial conditions at a computational grid
!
!  Method
!
!     The initial conditions are given using the following equation       30.70
!     for dimensionless Hs as function of dimensionless fetch:            30.70
!
!     Hs = 0.00288 f**(0.45)                                              40.00
!     Tp = 0.46    f**(0.27)                                              40.00
!     average direction = wind direction
!     directional distribution: Cos**2
!
!     after computation of the integral parameters the subroutine SSHAPE  30.70
!     is used to compute the spectrum
!
!***********************************************************************
!
  USE OCPCOMM4                                                        
  USE SWCOMM1                                                         
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
  USE SWCOMM4 
  USE ALL_VARS   
  USE VARS_WAVE, ONLY : SERIAL,PAR,ART,   &
                       IOBCN_W,I_OBC_N_W,AC2,COMPDA
                                                     
# if defined (EXPLICIT)
  USE MOD_ACTION_EX
# else            
  USE MOD_ACTION_IM
# endif      
# if defined (MULTIPROCESSOR)
  USE MOD_PAR
# endif
 
  IMPLICIT NONE

  REAL    :: SPCDIR(MDC,6)                                               
  REAL    :: SPCSIG(MSC)                                                 
  REAL    :: FDLSS,  TPDLSS, HSDLSS                                      
  INTEGER :: IN,IENT,IG,INX,JJ,IS,ID,IERR,IOB,IOB_NODE
  REAL    :: TLEN,TDXY,COSYG,FETCH,WX,WY,WSLOC,WDLOC
  REAL(SP), ALLOCATABLE :: ART_TMP(:)
!
  REAL    :: SPPARM_TMP(4)
!
  LOGICAL  INTERN

# if defined (MULTIPROCESSOR)
  REAL(SP), ALLOCATABLE :: AC2_TMP(:)
  INTEGER :: ISC,IDC,IP
# endif
!
! *** Fetch Computation 'the mean delta' ***
  TLEN = 0.0
	
  ALLOCATE(ART_TMP(NGL)); ART_TMP = 0.0
  IF(SERIAL) ART_TMP = ART
# if defined (MULTIPROCESSOR)
!-------------------Jianzhong---------------------------
!  IF(PAR)CALL GATHER(LBOUND(ART,1), UBOUND(ART,1), N,NGL,1,MYID,NPROCS, &
!                     EMAP,ART, ART_TMP)	
  IF(PAR)CALL ACOLLECT(MYID,MSRID,NPROCS,EMAP,ART,ART_TMP)
!-------------------------------------------------------
# endif		    
  DO IN = 1, NGL
    TLEN = TLEN + ART_TMP(IN)
  END DO
  TLEN = SQRT(TLEN)
  DEALLOCATE(ART_TMP)
  TDXY  = SQRT(REAL(NGL))                                     

  FETCH = TLEN/TDXY
# if defined (MULTIPROCESSOR)
  CALL MPI_BCAST(FETCH,1,MPI_F,0,MPI_FVCOM_GROUP,IERR)      
  CALL MPI_BCAST(TLEN,1,MPI_F,0,MPI_FVCOM_GROUP,IERR)      
# endif
  SY0   = 3.3                                                         

  SPPARM_TMP(1:4) = SPPARM(1:4)
                                                                    
  DO IG = 1, M         
!   check if the point is a true internal point                     
    INTERN = .TRUE.
	  
    DO IOB = 1,IOBCN_W
      IOB_NODE = I_OBC_N_W(IOB)
      IF(IG == IOB_NODE) INTERN = .FALSE.
    END DO  
	  
    IF(INTERN)THEN                                                
      TESTFL = .FALSE.                                              
      IF (VARWI) THEN
        WX  = COMPDA(IG,JWX2)
        WY  = COMPDA(IG,JWY2)
!
!       *** Local wind speed and direction ***
        WSLOC = SQRT(WX*WX + WY*WY)
        IF(WX /= 0. .OR. WY /= 0.)THEN
          WDLOC = ATAN2(WY,WX)
        ELSE
          WDLOC = 0.
        END IF
      ELSE
!     uniform wind field
        WSLOC = U10                                                
        WDLOC = WDIP                                               
      END IF
!
      IF(WSLOC > 1.E-10)THEN                                  
!
! Dimensionless Hs and Tp calculated according to K.K. Kahma & C.J. Calkoen,
! (JPO, 1992) and Pierson-Moskowitz for limit values.
!
!       calculate dimensionless fetch:
        FDLSS = GRAV_W * FETCH / (WSLOC*WSLOC)                       
!
!       calculate dimensionless significant wave height:
        HSDLSS = MIN (0.21, 0.00288*FDLSS**0.45)                   
        SPPARM(1) = HSDLSS * WSLOC**2 / GRAV_W
!       calculate dimensionless peak period:
        TPDLSS = MIN (1./0.13, 0.46*FDLSS**0.27)                   
        SPPARM(2) = WSLOC * TPDLSS / GRAV_W                          
        IF(SPPARM(1) < 0.05) SPPARM(2) = 2.                      
        SPPARM(3) = 180. * WDLOC / PI_W
        SPPARM(4) = 2.
      ELSE
        SPPARM(1) = 0.02
        SPPARM(2) = 2.                                             
        SPPARM(3) = 0.
        SPPARM(4) = 0.
      END IF
      CALL SSHAPE (AC2(1,1,IG), SPCSIG, SPCDIR, 2, 2) 
    END IF
  END DO
  
# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    ALLOCATE(AC2_TMP(0:MT));   AC2_TMP = 0.0
    DO ISC = 1,MSC
      DO IDC = 1,MDC
        DO IP = 1,MT
          AC2_TMP(IP)=AC2(IDC,ISC,IP)
        END DO  
        CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,AC2_TMP)
        CALL AEXCHANGE(NC,MYID,NPROCS,AC2_TMP)
        DO IP = 1,MT
          AC2(IDC,ISC,IP) = AC2_TMP(IP)
        END DO
      END DO
    END DO
    DEALLOCATE(AC2_TMP)
  END IF    	  

!JQI  IF(PAR)CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
# endif
      
  SPPARM(1:4) = SPPARM_TMP(1:4)

  RETURN
  END SUBROUTINE SWINCO
 

!*********************************************************************
!                                                                    *
!  SUBROUTINE SNEXTI (BSPECS ,BGRIDP ,AC1    ,SPCSIG ,SPCDIR ,        &
!JQIJQI  SUBROUTINE SNEXTI (SPCSIG ,SPCDIR ,DEPTH ,WLEVL ,FRIC ,UXB ,UYB ,WXI ,WYI)                     
  SUBROUTINE SNEXTI (SPCSIG ,SPCDIR ,DEPTH ,WLEVL ,FRIC ,UXB ,UYB)                     
! (This subroutine has not been fully tested. The codes not tested 
!  are commented using !JQI)
!                                                                    *
!*********************************************************************
!
!     Updates boundary conditions and input fields
!
!*********************************************************************
!
  USE TIMECOMM                                                    
  USE OCPCOMM4                                                   
  USE SWCOMM1                                                    
  USE SWCOMM2                                                    
  USE SWCOMM3                                                    
  USE SWCOMM4                                                    
  USE M_BNDSPEC                                                  
  USE VARS_WAVE
# if defined (WAVE_CURRENT_INTERACTION)
  USE MOD_WAVE_CURRENT_INTERACTION
# endif  

# if defined (MULTIPROCESSOR)
  USE MOD_PAR
  USE MOD_NESTING, ONLY : NESTING_ON_WAVE, NESTING_TYPE_WAVE, SET_VAR_WAVE,SET_VAR_WAVE_AC2
# endif
      
  IMPLICIT NONE                                         

!  INTEGER :: BGRIDP(*)      
!  REAL    :: AC1(MDC,MSC,0:MT)                                         
!JQI  REAL    :: BSPECS(MDC,MSC,NBSPEC,2)                                   
  REAL    :: SPCDIR(MDC,6)                                              
  REAL    :: SPCSIG(MSC)                                                
!JQIJQI  REAL    :: DEPTH(*), WLEVL(*), FRIC(*), UXB(*), UYB(*), WXI(*), WYI(*)                                             
  REAL    :: DEPTH(*), WLEVL(*), FRIC(*), UXB(*), UYB(*)                                             

  LOGICAL STPNOW                                                      
!
!
  CHARACTER   PTYPE*1                                                 
  INTEGER     IERR
  REAL        TSTVAL(10)                                              
  TYPE(BSPCDAT), POINTER :: CURBFL
      
  INTEGER :: IOB,IOB_NODE,IBGRID,INDXGR,IS,ID,IG
  INTEGER :: K1,K2
  REAL    :: W1,W2,ETOT,SIG2,HS,XP,YP,DEP,WLVL,DEPW,UU,VV
  REAL    :: VTOT,CGMAX,CGFACT
  
  REAL(SP), ALLOCATABLE :: AC2_TTMP(:,:,:)
      
!
! **   All action densities are shifted from array T+DTW
! **   to the array at time T
!
  IF(NSTATC == 1)THEN                                               
    IF(ITERMX > 1 .OR. PROPSC == 3)THEN             
      DO IG = 1, M                                               
        DO IS = 1, MSC
          DO ID = 1, MDC
            AC1(ID,IS,IG) = AC2(ID,IS,IG)   
	  END DO
	END DO
      END DO	    
    END IF                                                             
  END IF

! --- update boundary conditions
  DO IOB = 1, IOBCN_W
    IOB_NODE = I_OBC_N_W(IOB)
    CALL SSHAPE(AC2(1,1,IOB_NODE),SPCSIG,SPCDIR,FSHAPE,DSHAPE)
  END DO    
	   
# if defined (MULTIPROCESSOR)
  IF(NESTING_ON_WAVE)THEN
   IF(NESTING_TYPE_WAVE == 'wave parameters')THEN
  
    CALL SET_VAR_WAVE(WaveTime,HSC1=HSC1)
    CALL SET_VAR_WAVE(WaveTime,TPEAK=TPEAK)
    CALL SET_VAR_WAVE(WaveTime,DIRDEG1=DIRDEG1)

    DO IOB = 1, IOBCN_W
      IOB_NODE = I_OBC_N_W(IOB)
      SPPARM(1) = HSC1(IOB_NODE)
      SPPARM(2) = TPEAK(IOB_NODE)
      SPPARM(3) = DIRDEG1(IOB_NODE)
      SPPARM(4) = 20.0_SP 
      CALL SSHAPE(AC2(1,1,IOB_NODE),SPCSIG,SPCDIR,FSHAPE,DSHAPE)
    END DO    
	   
   ELSE IF(NESTING_TYPE_WAVE == 'spectral density')THEN
  
    ALLOCATE(AC2_TTMP(MDC,MSC,0:MT));  AC2_TTMP = 0.0_SP
    CALL SET_VAR_WAVE_AC2(WaveTime,AC2=AC2_TTMP)
	 
    DO IOB = 1, IOBCN_W
      IOB_NODE = I_OBC_N_W(IOB)
      AC2(:,:,IOB_NODE) = AC2_TTMP(:,:,IOB_NODE)
    END DO
    DEALLOCATE(AC2_TTMP)  
  
   ELSE
    CALL FATAL_ERROR("THE PARAMETER NESTING_TYPE_WAVE SHOULD BE 'wave parameters' or 'spectral density'")
   END IF
  END IF
# endif

!     --- determine spectra on boundary points of grid

!JQI    IF ( NBGRPT.GT.0 ) THEN
!JQI      IF (ITEST.GE.80) WRITE(PRTEST,*) ' number of boundary points ', &
!JQI	                                  NBGRPT
!JQI      DO IBGRID = 1, NBGRPT
!JQI        INDXGR = BGRIDP(6*IBGRID-5)
!JQI        IF (BGRIDP(6*IBGRID-4).EQ.1) THEN
!             obtain spectrum in boundary point from interpolation in space
!JQI          W1 = 0.001 * REAL(BGRIDP(6*IBGRID-3))
!JQI          K1 = BGRIDP(6*IBGRID-2)
!JQI          W2 = 1.-W1
!JQI          K2 = BGRIDP(6*IBGRID)
!JQI          CALL SINTRP (W1, W2, BSPECS(1,1,K1,1), BSPECS(1,1,K2,1),    &
!JQI	                   AC2(1,1,INDXGR), SPCDIR, SPCSIG)
!           --- store Hs from boundary condition in array HSOBND
!JQI          ETOT = 0.
!JQI          DO IS = 1, MSC
!JQI            SIG2 = SPCSIG(IS) ** 2
!JQI            DO ID = 1, MDC
!JQI              ETOT = ETOT + SIG2 * AC2(ID,IS,INDXGR)
!JQI            ENDDO
!JQI          ENDDO
!JQI          IF (ETOT.GT.0.) THEN
!JQI            HS = 4. * SQRT(FRINTF*DDIR*ETOT)
!JQI          ELSE
!JQI            HS = 0.
!JQI          ENDIF
!JQI          COMPDA(INDXGR,JHSIBC) = HS
!JQI        END IF
!JQI      END DO
!JQI    END IF
!
!     --- update input fields (wind, water level etc.)
!
!     fields 5 and 6: wind
      IF(IFLDYN(5) == 1)THEN
!JQIJQI        CALL FLFILE1(5, 6, WXI, WYI, 0, JWX2, JWX3, 0,           &
# if !defined (WAVE_ONLY)
        CALL FLFILE1(5, 6, 0, JWX2, JWX3, 0,           &
	             JWY2, JWY3, COSWC, SINWC, IERR)
# else
        IF(ICEIN_ON)THEN
	  CALL FLFILE1_ICE(5, 6, 0, JWX2, JWX3, 0,           &
	               JWY2, JWY3, COSWC, SINWC, IERR)
        ELSE
          CALL FLFILE1(5, 6, 0, JWX2, JWX3, 0,           &
	               JWY2, JWY3, COSWC, SINWC, IERR)
	END IF		       
# endif
      ELSE 
        WaveTime = WaveTime + IMDTW
      END IF

!     field 4: friction coeff.
      IF(IFLDYN(4) == 1)THEN
!JQI     CALL FLFILE ( 4, 0, FRIC, 0.,                            &
!JQI	               0, JFRC2, JFRC3, 0, 0, 0,                  &
!JQI		       1., 0.,                                    &
!JQI		       COMPDA, XCGRID, YCGRID,                    &
!JQI		       KGRPNT, IERR)
      ENDIF

!     field 7: water level
      IF(IFLDYN(7) == 1)THEN
!JQI     CALL FLFILE ( 7, 0, WLEVL, 0.,                            &
!JQI	               JWLV1, JWLV2, JWLV3, 0, 0, 0,               &
!JQI		       1., 0.,                                     &
!JQI		       COMPDA, XCGRID, YCGRID,                     &
!JQI		       KGRPNT, IERR)
!        Add bottom level to obtain depth
!JQI     DO IG = 1, MCGRD
!JQI       IF(IG > 1)THEN
!JQI         XP = XCGRID(IG)
!JQI         YP = YCGRID(IG)
!JQI         DEP = SVALQI (XP, YP, 1, DEPTH, 1 ,IX ,IY) 
!JQI         DEP = DEPTH(IG)                 
!JQI         COMPDA(IG,JDP1) = COMPDA(IG,JDP2)
!JQI         WLVL = COMPDA(IG,JWLV2)
!JQI         DEPW = DEP + WLVL + WLEV
!JQI         COMPDA(IG,JDP2) = DEPW
!JQI         IF (LSETUP > 0) THEN                                       
!JQI           COMPDA(IG,JDPSAV) = COMPDA(IG,JDP2)                   
!JQI         ENDIF                                                       
!JQI       ENDIF
!JQI     END DO
        !LWU COMPDA(IG,JDP2) = D(IG)
      ENDIF
#  if  defined (WAVE_CURRENT_INTERACTION)
       !if(msr) print*,'add fvcom elevaton replace in COMPDA'
!JQI       if(ITW>=360) then
!!JQIJQI       IF(IINT >= IRAMP) THEN
        DO IG = 1, MT
        COMPDA(IG,JDP1) = COMPDA(IG,JDP2)
        COMPDA(IG,JDP2) = D(IG)
        DEPTH(IG) = H(IG)
        ENDDO
!!JQIJQI       endif
#  endif
!     field 2 and 3: current velocity
      IF(IFLDYN(2) == 1)THEN
!JQI     CALL FLFILE ( 2, 3, UXB, UYB,                                    &
!JQI	               JVX1, JVX2, JVX3, JVY1, JVY2, JVY3,                &
!JQI		       COSVC, SINVC,                                      &
!JQI		       COMPDA, XCGRID, YCGRID, KGRPNT, IERR)
!       reduce current velocity if Froude number is larger than PNUMS(18)
!JQI     DO IG = 1, MCGRD                                                  
!JQI       IF (IG > 1) THEN                                           
!JQI         DEPW = COMPDA(IG,JDP2)
!JQI         IF (DEPW > 0.) THEN
!JQI           UU = COMPDA(IG,JVX2)
!JQI           VV = COMPDA(IG,JVY2)
!JQI           VTOT = SQRT (UU*UU + VV*VV)
!JQI           CGMAX = PNUMS(18)*SQRT(GRAV_W*DEPW)
!JQI           IF (VTOT > CGMAX) THEN
!JQI             CGFACT = CGMAX / VTOT
!JQI             COMPDA(IG,JVX2) = UU * CGFACT
!JQI             COMPDA(IG,JVY2) = VV * CGFACT
!                write IX,IY to error points file
!JQI             IF (ERRPTS > 0 .AND. INODE == MASTER) THEN               
!                  WRITE (ERRPTS, 211) IX+MXF-1, IY+MYF-1, 1
! 211              FORMAT (I4, 1X, I4, 1X, I2)
!JQI             ENDIF
!JQI           ENDIF
!JQI         ENDIF
!JQI       ENDIF                                                         
!JQI     ENDDO                                                           
      ENDIF
# if defined (WAVE_CURRENT_INTERACTION)
      !if(msr) print*, 'add fvcom current in COMPDA'
      COMPDA(:,JVX2) = 0.0_SP
      COMPDA(:,JVY2) = 0.0_SP
      CALL CURRENT2WAVE
      !PRINT*,'DOPPLER==',UDOP(1),VDOP(1)
!JQI      IF(ITW>=360) then
      IF(IINT >= IRAMP) THEN
       DO IG =1,MT
          COMPDA(IG,JVX2) = UDOP(IG)
          COMPDA(IG,JVY2) = VDOP(IG)
       ENDDO
      ENDIF
# endif
      RETURN
      END SUBROUTINE SNEXTI

!***********************************************************************
!                                                                      *
  SUBROUTINE SWRBC                                          

! (This subroutine has not been fully tested. The codes not tested
!  are commented using !JQI)
!                                                                      *
!***********************************************************************
!
  USE OCPCOMM4                                                        
  USE SWCOMM1                                                         
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE M_GENARR
  USE MOD_PREC
  USE VARS_WAVE, ONLY : M,MT,COMPDA,WaveTime,UWWIND,VWWIND,NTVE,NBVE
  USE MOD_FORCE, ONLY : UPDATE_WIND2WAVE
  USE CONTROL, ONLY : WIND_ON
  USE ALL_VARS, ONLY : UUWIND,VVWIND
# if defined (VEGETATION)
  USE MOD_VEGETATION, ONLY : VEG,nveg,pdens
# endif
      
  IMPLICIT NONE

  INTEGER :: IENT,IG,INDX,JJ
  REAL    :: DEP,DEPW,UU,VV,VTOT,CGMAX,CGFACT
      
!
!     *** The arrays start to be filled in the second value ***
!     *** because in COMPDA(1,"variable"), is the default   ***
!     *** value for land points    version 30.21            ***
!
!     ***  Default values for land point  ***
  COMPDA(1,JDP2) = -1.                                                
  IF(JDP1 > 1) COMPDA(1,JDP1) = -1.                                 
  IF(JDP3 > 1) COMPDA(1,JDP3) = -1.                                 
  IF(VARWLV)THEN
    COMPDA(1,JWLV2) = 0.                                              
!   next two lines only for nonstat water level
    IF(JWLV1 > 1) COMPDA(1,JWLV1) = 0.                              
    IF(JWLV3 > 1) COMPDA(1,JWLV3) = 0.                              
  END IF
  IF(ICUR > 0)THEN
    COMPDA(1,JVX2) = 0.                                               
    COMPDA(1,JVY2) = 0.                                               
    IF(JVX1 > 1) COMPDA(1,JVX1) = 0.                                
    IF(JVY1 > 1) COMPDA(1,JVY1) = 0.                                
    IF(JVX3 > 1) COMPDA(1,JVX3) = 0.                                
    IF(JVY3 > 1) COMPDA(1,JVY3) = 0.                                
  END IF
  IF(VARFR)THEN
    COMPDA(1,JFRC2) = 0.                                              
    COMPDA(1,JFRC3) = 0.                                              
  END IF
  IF(VARWI)THEN
    COMPDA(1,JWX2) = 0.                                               
    COMPDA(1,JWY2) = 0.                                               
    IF(JWX3 > 1) COMPDA(1,JWX3) = 0.                                
    IF(JWY3 > 1) COMPDA(1,JWY3) = 0.                                
  END IF
  IF(VARAST)THEN
    COMPDA(1,JASTD2) = 0.                                             
    COMPDA(1,JASTD3) = 0.                                             
  END IF
# if defined (VEGETATION)
  IF (VARNPL) THEN
    COMPDA(1,JNPLA2) = 0.
    COMPDA(1,JNPLA3) = 0.
  ENDIF
# endif

  IF(VARWI) CALL UPDATE_WIND2WAVE(WaveTime,UWWIND,VWWIND)

  DO INDX = 1, M      !MT      
!
!   ***** compute depth and water level *****
!
    DEP = DEPTH(INDX)

    IF(VARWLV)THEN                                                
!JQI      WLVL = SVALQI (XP, YP, 7, WLEVL, 1 ,IX ,IY)                   
!JQI      COMPDA(INDX,JWLV2) = WLVL
!JQI      IF (JWLV1.GT.1) COMPDA(INDX,JWLV1) = WLVL                     
!JQI      IF (JWLV3.GT.1) COMPDA(INDX,JWLV3) = WLVL                     
!JQI      DEP = DEP + WLVL
    END IF
!   add constant water level
    DEPW = DEP + WLEV                                               
    COMPDA(INDX,JDP2) = DEPW
!   ***In this step the water level at T+DTW is copied to the  ***
!   ***water level at T (Only for first time computation)     ***
    IF(JDP1 > 1) COMPDA(INDX,JDP1) = DEPW                         
    IF(JDP3 > 1) COMPDA(INDX,JDP3) = DEPW                         
!
!   ***** compute current velocity *****
!
    IF(ICUR == 1)THEN                      
      IF(DEPW > 0.)THEN
        UU  = UXB(INDX,1)
        VV  = UYB(INDX,1)
        VTOT = SQRT (UU*UU + VV*VV)
        CGMAX = PNUMS(18)*SQRT(GRAV_W*DEPW)
        IF(VTOT > CGMAX)THEN
          CGFACT = CGMAX / VTOT
          UU = UU * CGFACT
          VV = VV * CGFACT
        END IF
        COMPDA(INDX,JVX2) =  UU
        COMPDA(INDX,JVY2) =  VV
      ELSE
        COMPDA(INDX,JVX2) =  0.
        COMPDA(INDX,JVY2) =  0.
      END IF
      IF(JVX1 > 1) COMPDA(INDX,JVX1) = COMPDA(INDX,JVX2)          
      IF(JVY1 > 1) COMPDA(INDX,JVY1) = COMPDA(INDX,JVY2)          
      IF(JVX3 > 1) COMPDA(INDX,JVX3) = COMPDA(INDX,JVX2)          
      IF(JVY3 > 1) COMPDA(INDX,JVY3) = COMPDA(INDX,JVY2)          
    END IF
!
!     ***** compute variable friction coefficient *****
!
    IF(VARFR)THEN
!JQI     FRI = SVALQI (XP, YP, 4, FRIC, 1 ,IX ,IY)                    
!JQI     COMPDA(INDX,JFRC2) = FRI                                     
!JQI     COMPDA(INDX,JFRC3) = FRI                                     
    END IF
!
!     ***** compute variable wind velocity *****
!
    IF(VARWI)THEN
      UU = 0.0_SP
      VV = 0.0_SP
      DO JJ = 1, NTVE(INDX)
        UU = UU + UWWIND(NBVE(INDX,JJ))
        VV = VV + VWWIND(NBVE(INDX,JJ))
!JQIJQI        UU = UU + UUWIND(NBVE(INDX,JJ))
!JQIJQI        VV = VV + VVWIND(NBVE(INDX,JJ))
      END DO
      UU = UU/NTVE(INDX)	
      VV = VV/NTVE(INDX)	
!      UU  = WXI(INDX,1)                     
!      VV  = WYI(INDX,1)                     
      COMPDA(INDX,JWX2) =  UU
      COMPDA(INDX,JWY2) =  VV
      IF(JWX3 > 1) COMPDA(INDX,JWX3) = COMPDA(INDX,JWX2)          
      IF(JWY3 > 1) COMPDA(INDX,JWY3) = COMPDA(INDX,JWY2)          
    END IF
!
!   ***** compute variable air-sea temperature difference *****         
!
    IF(VARAST)THEN
!JQI             ASTD = SVALQI (XP, YP, 10, ASTDF, 1 ,IX ,IY)                 
!JQI             COMPDA(INDX,JASTD2) = ASTD                                   
!JQI             COMPDA(INDX,JASTD3) = ASTD                                   
    END IF
  END DO
!
!   *** initialise setup and saved depth ***                            
!
  IF(LSETUP > 0)THEN                                               
    DO INDX = 1, M                                                  
      COMPDA(INDX,JSETUP) =  0.                                       
      COMPDA(INDX,JDPSAV) = COMPDA(INDX,JDP2)                         
    END DO                                                             
  END IF                                                               
!
!   --- initialize HSIBC                                                
  COMPDA(:,JHSIBC) = 0.                                               
!
!   --- initialize ZELEN and USTAR                                      
  COMPDA(:,JZEL  ) = 2.E-33                                           
  COMPDA(:,JUSTAR) = 1.E-15                                           
!
  RETURN
  END SUBROUTINE SWRBC

!***********************************************************************
# if defined (WAVE_ONLY)
!                                                                      *
  SUBROUTINE SWRBC_ICE                                          

! (This subroutine has not been fully tested. The codes not tested
!  are commented using !JQI)
!                                                                      *
!***********************************************************************
!
  USE OCPCOMM4                                                        
  USE SWCOMM1                                                         
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE M_GENARR
  USE MOD_PREC
  USE VARS_WAVE, ONLY : M,MT,COMPDA,WaveTime,UWWIND,VWWIND,NTVE,NBVE
  USE MOD_FORCE, ONLY : UPDATE_WIND2WAVE
  USE CONTROL, ONLY : WIND_ON
  USE ALL_VARS, ONLY : UUWIND,VVWIND,CICE_TEST
  USE MOD_NORTHPOLE    
  IMPLICIT NONE

  INTEGER :: IENT,IG,INDX,JJ
  REAL    :: DEP,DEPW,UU,VV,VTOT,CGMAX,CGFACT
!========================zhangyang=========================================
  REAL    :: DIR1, DIR2, UUVV,DIFF,CICE_TMP
!========================zhangyang=========================================
      
!
!     *** The arrays start to be filled in the second value ***
!     *** because in COMPDA(1,"variable"), is the default   ***
!     *** value for land points    version 30.21            ***
!
!     ***  Default values for land point  ***
  COMPDA(1,JDP2) = -1.                                                
  IF(JDP1 > 1) COMPDA(1,JDP1) = -1.                                 
  IF(JDP3 > 1) COMPDA(1,JDP3) = -1.                                 
  IF(VARWLV)THEN
    COMPDA(1,JWLV2) = 0.                                              
!   next two lines only for nonstat water level
    IF(JWLV1 > 1) COMPDA(1,JWLV1) = 0.                              
    IF(JWLV3 > 1) COMPDA(1,JWLV3) = 0.                              
  END IF
  IF(ICUR > 0)THEN
    COMPDA(1,JVX2) = 0.                                               
    COMPDA(1,JVY2) = 0.                                               
    IF(JVX1 > 1) COMPDA(1,JVX1) = 0.                                
    IF(JVY1 > 1) COMPDA(1,JVY1) = 0.                                
    IF(JVX3 > 1) COMPDA(1,JVX3) = 0.                                
    IF(JVY3 > 1) COMPDA(1,JVY3) = 0.                                
  END IF
  IF(VARFR)THEN
    COMPDA(1,JFRC2) = 0.                                              
    COMPDA(1,JFRC3) = 0.                                              
  END IF
  IF(VARWI)THEN
    COMPDA(1,JWX2) = 0.                                               
    COMPDA(1,JWY2) = 0.                                               
    IF(JWX3 > 1) COMPDA(1,JWX3) = 0.                                
    IF(JWY3 > 1) COMPDA(1,JWY3) = 0.                                
  END IF
  IF(VARAST)THEN
    COMPDA(1,JASTD2) = 0.                                             
    COMPDA(1,JASTD3) = 0.                                             
  END IF

  IF(VARWI) CALL UPDATE_WIND2WAVE(WaveTime,UWWIND,VWWIND)

  DO INDX = 1, M      !MT      
!
!   ***** compute depth and water level *****
!
    DEP = DEPTH(INDX)

    IF(VARWLV)THEN                                                
!JQI      WLVL = SVALQI (XP, YP, 7, WLEVL, 1 ,IX ,IY)                   
!JQI      COMPDA(INDX,JWLV2) = WLVL
!JQI      IF (JWLV1.GT.1) COMPDA(INDX,JWLV1) = WLVL                     
!JQI      IF (JWLV3.GT.1) COMPDA(INDX,JWLV3) = WLVL                     
!JQI      DEP = DEP + WLVL
    END IF
!   add constant water level
    DEPW = DEP + WLEV                                               
    COMPDA(INDX,JDP2) = DEPW
!   ***In this step the water level at T+DTW is copied to the  ***
!   ***water level at T (Only for first time computation)     ***
    IF(JDP1 > 1) COMPDA(INDX,JDP1) = DEPW                         
    IF(JDP3 > 1) COMPDA(INDX,JDP3) = DEPW                         
!
!   ***** compute current velocity *****
!
    IF(ICUR == 1)THEN                      
      IF(DEPW > 0.)THEN
        UU  = UXB(INDX,1)
        VV  = UYB(INDX,1)
        VTOT = SQRT (UU*UU + VV*VV)
        CGMAX = PNUMS(18)*SQRT(GRAV_W*DEPW)
        IF(VTOT > CGMAX)THEN
          CGFACT = CGMAX / VTOT
          UU = UU * CGFACT
          VV = VV * CGFACT
        END IF
        COMPDA(INDX,JVX2) =  UU
        COMPDA(INDX,JVY2) =  VV
      ELSE
        COMPDA(INDX,JVX2) =  0.
        COMPDA(INDX,JVY2) =  0.
      END IF
      IF(JVX1 > 1) COMPDA(INDX,JVX1) = COMPDA(INDX,JVX2)          
      IF(JVY1 > 1) COMPDA(INDX,JVY1) = COMPDA(INDX,JVY2)          
      IF(JVX3 > 1) COMPDA(INDX,JVX3) = COMPDA(INDX,JVX2)          
      IF(JVY3 > 1) COMPDA(INDX,JVY3) = COMPDA(INDX,JVY2)          
    END IF
!
!     ***** compute variable friction coefficient *****
!
    IF(VARFR)THEN
!JQI     FRI = SVALQI (XP, YP, 4, FRIC, 1 ,IX ,IY)                    
!JQI     COMPDA(INDX,JFRC2) = FRI                                     
!JQI     COMPDA(INDX,JFRC3) = FRI                                     
    END IF
!
!     ***** compute variable wind velocity *****
!
    IF(VARWI)THEN
      UU = 0.0_SP
      VV = 0.0_SP
      DO JJ = 1, NTVE(INDX)
!=========================zhangyang=====start==========================
        DIR1=ATAN2(VWWIND(NBVE(INDX,JJ)),UWWIND(NBVE(INDX,JJ)))
        DIFF=VX(INDX)-XC(NBVE(INDX,JJ))
      IF(DIFF >  180.0)THEN
      DIFF = -360.0+DIFF
      ELSE IF(DIFF < -180.0)THEN
      DIFF =  360.0+DIFF
      END IF
      DIR2=DIR1-DIFF*DEG2RAD
    IF(DIR2>(PI*2)) DIR2=DIR2-(PI*2)
    IF(DIR2<0) DIR2=DIR2+(PI*2)
    UUVV=SQRT(VWWIND(NBVE(INDX,JJ))**2+UWWIND(NBVE(INDX,JJ))**2)
    UU = UU+UUVV*COS(DIR2)
    VV = VV+UUVV*SIN(DIR2)



        !UU = UU + UWWIND(NBVE(INDX,JJ))
        !VV = VV + VWWIND(NBVE(INDX,JJ))
!JQIJQI        UU = UU + UUWIND(NBVE(INDX,JJ))
!JQIJQI        VV = VV + VVWIND(NBVE(INDX,JJ))
!IF(INDX==2435.OR.INDX==2465.OR.INDX==2439.OR.INDX==2437)PRINT *,
!"INDX",INDX,"UU",UU,"VV",VV
      END DO
!=========================zhangyang==end=============================
IF(ICEIN_ON)THEN
      CICE_TMP=CICE_TEST(INDX)
ELSE
      CICE_TMP=0.0
END IF
      UU = UU/NTVE(INDX)*(1.-CICE_TMP)
      VV = VV/NTVE(INDX)*(1.-CICE_TMP)
!=========================zhangyang=====start==========================
    IF(INDX==NODE_NORTHPOLE)THEN
    UU = 0.0_SP
    VV = 0.0_SP
!PRINT *, INDX    
    DO JJ = 1,NTVE(INDX)
    DIR1=ATAN2(VWWIND(NBVE(INDX,JJ)),UWWIND(NBVE(INDX,JJ)))
    DIR2=DIR1+XC(NBVE(INDX,JJ))*DEG2RAD+(PI/2)
!    DIR2=DIR1+XC(NBVE(INDX,JJ))*DEG2RAD
    IF(DIR2>(PI*2)) DIR2=DIR2-(PI*2)
    IF(DIR2<0) DIR2=DIR2+(PI*2)
    UUVV=SQRT(VWWIND(NBVE(INDX,JJ))**2+UWWIND(NBVE(INDX,JJ))**2)
    UU = UU+UUVV*COS(DIR2)
    VV = VV+UUVV*SIN(DIR2)
!PRINT *,
!"DIR1",DIR1,"DIR2",DIR2,"XC",XC(NBVE(INDX,JJ)),"XCRAD",XC(NBVE(INDX,JJ))*DEG2RAD
!PRINT *, "UWINDSPEED",UWWIND(NBVE(INDX,JJ)),"VWINDSPEED",VWWIND(NBVE(INDX,JJ))
!PRINT *, "UUVV",UUVV,"UU",UU,"VV",VV
    END DO
     UU = UU/NTVE(INDX)*(1.-CICE_TMP)
     VV = VV/NTVE(INDX)*(1.-CICE_TMP)
    END IF
!=========================zhangyang==========end=====================

!      UU  = WXI(INDX,1)                     
!      VV  = WYI(INDX,1)                     
      COMPDA(INDX,JWX2) =  UU
      COMPDA(INDX,JWY2) =  VV
      IF(JWX3 > 1) COMPDA(INDX,JWX3) = COMPDA(INDX,JWX2)          
      IF(JWY3 > 1) COMPDA(INDX,JWY3) = COMPDA(INDX,JWY2)          
    END IF
!
!   ***** compute variable air-sea temperature difference *****         
!
    IF(VARAST)THEN
!JQI             ASTD = SVALQI (XP, YP, 10, ASTDF, 1 ,IX ,IY)                 
!JQI             COMPDA(INDX,JASTD2) = ASTD                                   
!JQI             COMPDA(INDX,JASTD3) = ASTD                                   
    END IF
  END DO

# if defined (VEGETATION)
!
!     ***** compute number of plants per square meter *****               40.55
!
          IF (VARNPL) THEN
!JZGE        XNPL = SVALQI (XP, YP, 11, NPLAF, 1 ,IX ,IY)               
             COMPDA(INDX,JNPLA2) = sum(VEG % plant(INDX,1:nveg,pdens))
             COMPDA(INDX,JNPLA3) = sum(VEG % plant(INDX,1:nveg,pdens))
          ENDIF
# endif
!
!   *** initialise setup and saved depth ***                            
!
  IF(LSETUP > 0)THEN                                               
    DO INDX = 1, M                                                  
      COMPDA(INDX,JSETUP) =  0.                                       
      COMPDA(INDX,JDPSAV) = COMPDA(INDX,JDP2)                         
    END DO                                                             
  END IF                                                               
!
!   --- initialize HSIBC                                                
  COMPDA(:,JHSIBC) = 0.                                               
!
!   --- initialize ZELEN and USTAR                                      
  COMPDA(:,JZEL  ) = 2.E-33                                           
  COMPDA(:,JUSTAR) = 1.E-15                                           
!
  RETURN
  END SUBROUTINE SWRBC_ICE
# endif
# endif
