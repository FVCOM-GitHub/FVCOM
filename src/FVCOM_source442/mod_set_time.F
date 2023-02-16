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

MODULE MOD_SET_TIME
  USE CONTROL
  USE MOD_NCTOOLS
  USE MOD_INPUT
  USE MOD_TIME
# if defined (WAVE_CURRENT_INTERACTION)
  USE VARS_WAVE, ONLY : WaveTime
# endif  
  IMPLICIT NONE
  
  INTERFACE SET_STARTUP_FILE_STACK
     MODULE PROCEDURE SET_STARTUP_FILE_STACK_BY_TIME
     MODULE PROCEDURE SET_STARTUP_FILE_STACK_BY_CYCLE
  END INTERFACE
  
CONTAINS

  SUBROUTINE SETUP_TIME
    !===============================================================================!
    ! SETUP_TIME
    !===============================================================================!
    IMPLICIT NONE
    INTEGER STATUS
    TYPE(TIME)    :: RTIME
    TYPE(NCFILE), POINTER   :: NCF
    TYPE(NCVAR),  POINTER   :: VAR
    LOGICAL FOUND

    TYPE(TIME) :: INTERVAL_TIME, NEXT_TIME,PREV_TIME
    INTEGER :: STKLEN, STKMAX,PREVSTK,NEXTSTK

    INTEGER(itime) :: Dummy
    INTEGER :: I
    CHARACTER(LEN=4) :: BFLAG,EFLAG,FLAG


    ! SET THE TIME STEP
    CALL SET_MODEL_TIMESTEP !from runfile (imdte, dte, imdti, dti)


    SELECT CASE(STARTUP_TYPE)
       !=================================================
       ! HOTSTART 
    CASE("hotstart")
       !=================================================
       if(DBG_SET(dbg_log)) then 
          WRITE(IPT,*)'!               SET TIME FOR HOTSTART                            !'
          WRITE(IPT,*)'!                                                                !'
       end if
       IF(.NOT. ASSOCIATED(NC_START)) CALL FATAL_ERROR&
            & ('STARUP FILE IS NOT ASSOCIATED IN SETUP_TIME!')
       
 
       SELECT CASE(USE_REAL_WORLD_TIME)
       CASE(.TRUE.)

          ! GET THE START TIME
          StartTime = READ_DATETIME(START_DATE,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string START_DATE: "//trim(START_DATE))
          RUNFILE_STARTTIME = STARTTIME

          ! GET THE END TIME
          EndTime = READ_DATETIME(END_DATE,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string END_DATE: "&
               &//trim(END_DATE))
          ! SANITY CHECK
          if(StartTime .GT. EndTime) &
               & Call Fatal_Error("Runfile Start_Date exceeds or equal to End_Date")
          
          ! FIND THE START TIME IN THE RESTART FILE AND SET ISTART
          CALL SET_STARTUP_FILE_STACK(StartTime,IINT)

          ! ADVANCE THE TO THE FIRST TIME STEP OF THE MODEL FROM THE
          ! INITIAL CONDITION
          ISTART = IINT +1

          !CALCULATE THE NUMBER OF STEPS AND IEND
          NSTEPS = CALCULATE_NUMBER_OF_TIMESTEPS(StartTime,EndTime)
          IEND = ISTART + NSTEPS

          ! GET THE REFERENCE DATE TIME
          IF(DATE_REFERENCE /= 'default')THEN
	    ReferenceDate = READ_DATETIME(DATE_REFERENCE,DATE_FORMAT,TIMEZONE,status)
            if (status == 0 ) &
               & Call Fatal_Error("Could not read the date string DATE_REFERENCE: "//trim(DATE_REFERENCE))
	  ELSE
	    ReferenceDate%MJD   = 0
	    ReferenceDate%MuSOD = 0
	  END IF  

       CASE (.FALSE.)

          ! GET THE START AND END INFORMATION
          CALL IDEAL_TIME_STRING2TIME(START_DATE,BFLAG,StartTIME,IINT)
          CALL IDEAL_TIME_STRING2TIME(END_DATE,EFLAG,EndTIME,IEND)
          
          ! SANITY CHECK
          IF (BFLAG /= EFLAG) CALL FATAL_ERROR&
               ('IDEALIZED MODEL TIME SPECIFICATION IS INCORRENT',&
               &'BEGIN AND END CAN BE IN EITHER CYCLES OR TIME BUT NOT MIXED',&
               & trim(start_date),trim(end_date) )
          
          IF (BFLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

             !CALCULATE THE NUMBER OF STEPS
             NSTEPS = CALCULATE_NUMBER_OF_TIMESTEPS(StartTime,EndTime)

             ! FIND THE INITIAL CONDITION IN THE RESTART FILE
             CALL SET_STARTUP_FILE_STACK(StartTime,IINT)

             ! ADVANCE THE TO THE FIRST TIME STEP OF THE MODEL FROM THE
             ! INITIAL CONDITION
             ISTART = IINT +1

             ! GET IEND
             IEND = ISTART + NSTEPS
             
             
          ELSE IF(BFLAG == 'step') THEN ! IF START AND END IINT WERE SPECIFIED

             ! CALCULATE NSTEPS
             NSTEPS = IEND - IINT  
             
             ! SANITY CHECK
             IF(NSTEPS .LT. 0) CALL FATAL_ERROR&
                  &('Number of steps can not be less than zero')


             ! FIND THE START IINT CYCLE IN THE RESTART FILE AND GET
             ! THE START TIME
             CALL SET_STARTUP_FILE_STACK(IINT,StartTime)

             ! CALCULATE THE END TIME
             EndTime = StartTime + IMDTI * nSteps

             ! ADVANCE THE TO THE FIRST TIME STEP OF THE MODEL FROM THE
             ! INITIAL CONDITION
             ISTART = IINT + 1


          ELSE
             CALL FATAL_ERROR('IDEAL_TIME_STRING2TIME returned invalid flag')

          END IF

       END SELECT

       ! SET THE INITIAL MODEL TIME
       IntTime = StartTime
       ExtTime = StartTime
       RUNFILE_STARTTIME = STARTTIME

#      if defined (WAVE_CURRENT_INTERACTION)
       WaveTime = StartTime
#      endif  
          

       ! INITIALIZE IEXT
       IEXT = 1

       CALL REPORT_TIME_SETUP

! SET THE OUT PUT TIME FOR THE DIFFERENT FILES
       IF (NCAV_ON) THEN
       
          CALL GET_FIRST_OUTPUT_TIME(TRIM(NCAV_FIRST_OUT),NEXT_TIME)
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(NCAV_OUT_INTERVAL),INTERVAL_TIME)

          NEXTSTK = 0
          PREVSTK = 0 
          PREV_TIME = ZEROTIME
          STKLEN = 0
          STKMAX = NCAV_OUTPUT_STACK


          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TIME AVERAGE OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_AVG,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)


          if (STARTTIME > NEXT_TIME)&
               & Call Fatal_Error("The time of the first averaged output &
               &must be greater than or equal to the start time")
          
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first averaged output &
               &must be less than or equal to the end time")

       END IF

       IF (NC_ON) THEN
       
          CALL GET_FIRST_OUTPUT_TIME(TRIM(NC_FIRST_OUT),NEXT_TIME)
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(NC_OUT_INTERVAL),INTERVAL_TIME)

          NEXTSTK = 0
          PREVSTK = 0 
          PREV_TIME = ZEROTIME
          STKLEN = 0
          STKMAX = NC_OUTPUT_STACK


          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NETCDF OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_DAT,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)



          if (STARTTIME > NEXT_TIME)&
               & Call Fatal_Error("The time of the first file output &
               &must be greater than or equal to the start time")
          
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first file output &
               &must be less than or equal to the end time")

       END IF

       IF (NCSF_ON) THEN
       
          CALL GET_FIRST_OUTPUT_TIME(TRIM(NCSF_FIRST_OUT),NEXT_TIME)
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(NCSF_OUT_INTERVAL),INTERVAL_TIME)

          NEXTSTK = 0
          PREVSTK = 0 
          PREV_TIME = ZEROTIME
          STKLEN = 0
          STKMAX = NCSF_OUTPUT_STACK


          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! SURFACE NETCDF OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_SF,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)



          if (STARTTIME > NEXT_TIME)&
               & Call Fatal_Error("The time of the first surface file output &
               &must be greater than or equal to the start time")
          
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first surface file output &
               &must be less than or equal to the end time")

       END IF

       IF (RST_ON) THEN
       
 
       
          CALL GET_FIRST_OUTPUT_TIME(TRIM(RST_FIRST_OUT),NEXT_TIME)
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(RST_OUT_INTERVAL),INTERVAL_TIME)

          NEXTSTK = 0
          PREVSTK = 0 
          PREV_TIME = ZEROTIME
          STKLEN = 0
          STKMAX = RST_OUTPUT_STACK


          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! RESTART OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_RST,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)



          if (STARTTIME > NEXT_TIME)&
               & Call Fatal_Error("The time of the first restart output &
               &must be greater than or equal to the start time")
          
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first restart output &
               &must be less than or equal to the end time")

       END IF


!=================================================
! CRASHSTART 
    CASE("crashrestart") ! OR FORECAST MODE
!=================================================
       if(DBG_SET(dbg_log)) then 
          WRITE(IPT,*)'!               SET TIME FOR CRASHRESTART                        !'
          WRITE(IPT,*)'!                                                                !'
       end if

       IF(.NOT. ASSOCIATED(NC_START)) CALL FATAL_ERROR&
            & ('STARUP FILE IS NOT ASSOCIATED IN SETUP_TIME!')
       
       IF(FORECAST_MODE) THEN
          ! GET THE START TIME FROM THE NAMELIST
          StartTime = READ_DATETIME(START_DATE,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string START_DATE: "//trim(START_DATE))

          CALL SET_STARTUP_FILE_STACK(StartTime,IINT)

       ELSE
          ! OTHERWISE GET THE LAST TIME IN THE RESTART FILE
          CALL  SET_LAST_FILE_TIME(NC_START,StartTime, IINT)
       END IF

       SELECT CASE(USE_REAL_WORLD_TIME)
       CASE(.TRUE.)

          RUNFILE_StartTime = READ_DATETIME(START_DATE,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string START_DATE: "//trim(START_DATE))


          ! GET THE END TIME
          EndTime = READ_DATETIME(END_DATE,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string END_DATE: "&
               &//trim(END_DATE))

          ! SANITY CHECK
          if(StartTime .GT. EndTime) &
               & Call Fatal_Error("Runfile Start_Date exceeds or equal to End_Date")
          
          ! ADVANCE TO THE FIRST TIME STEP OF THE MODEL FROM THE
          ! INITIAL CONDITION
          ISTART = IINT +1

          !CALCULATE THE NUMBER OF STEPS AND IEND
          NSTEPS = CALCULATE_NUMBER_OF_TIMESTEPS(StartTime,EndTime)
          IEND = ISTART + NSTEPS

          ! GET THE REFERENCE DATE TIME
          IF(DATE_REFERENCE /= 'default')THEN
            ReferenceDate = READ_DATETIME(DATE_REFERENCE,DATE_FORMAT,TIMEZONE,status)
            if (status == 0 ) &
               & Call Fatal_Error("Could not read the date string DATE_REFERENCE: "//trim(DATE_REFERENCE))
	  ELSE
	    ReferenceDate%MJD   = 0
	    ReferenceDate%MuSOD = 0
	  END IF  

       CASE (.FALSE.)

          ! GET THE START AND END INFORMATION
          CALL IDEAL_TIME_STRING2TIME(START_DATE,BFLAG,RUNFILE_STARTTIME,ISTART)


          CALL IDEAL_TIME_STRING2TIME(END_DATE,EFLAG,EndTIME,IEND)

          ! SANITY CHECK
          IF (BFLAG /= EFLAG) CALL FATAL_ERROR&
               ('IDEALIZED MODEL TIME SPECIFICATION IS INCORRENT',&
               &'BEGIN AND END CAN BE IN EITHER CYCLES OR TIME BUT NOT MIXED',&
               & trim(start_date),trim(end_date) )
          
          IF (EFLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

             !CALCULATE THE NUMBER OF STEPS
             NSTEPS = CALCULATE_NUMBER_OF_TIMESTEPS(StartTime,EndTime)

             ! ADVANCE THE TO THE FIRST TIME STEP OF THE MODEL FROM THE
             ! INITIAL CONDITION
             ISTART = IINT +1

             ! GET IEND
             IEND = ISTART + NSTEPS
             
             
          ELSE IF(EFLAG == 'step') THEN ! IF START AND END IINT WERE SPECIFIED

             ! CALCULATE THE RUNFILE START TIME
             RUNFILE_STARTTIME = IMDTI * ISTART

             ! CALCULATE NSTEPS
             NSTEPS = IEND - IINT 
             
             ! SANITY CHECK
             IF(NSTEPS .LT. 0) CALL FATAL_ERROR&
                  &('Number of steps can not be less than zero')

             IF(ISTART .LT. 0) CALL FATAL_ERROR&
                  &('Starting time step can not be less than zero')

             ! CALCULATE THE END TIME
             EndTime = StartTime + IMDTI * nSteps

             ! ADVANCE THE TO THE FIRST TIME STEP OF THE MODEL FROM THE
             ! INITIAL CONDITION
             ISTART = IINT + 1

          ELSE
             CALL FATAL_ERROR('IDEAL_TIME_STRING2TIME returned invalid flag')

          END IF

       END SELECT

       ! SET THE INITIAL MODEL TIME
       IntTime = StartTime
       ExtTime = StartTime

#      if defined (WAVE_CURRENT_INTERACTION)
       WaveTime = StartTime
#      endif  

       ! INITIALIZE IEXT
       IEXT = 1

       CALL REPORT_TIME_SETUP

! SET THE OUT PUT TIME FOR THE DIFFERENT FILES
       IF (NCAV_ON) THEN
                 
          ! GET THE OUT PUT INTERVAL FROM RUN FILE
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(NCAV_OUT_INTERVAL),INTERVAL_TIME)

          ! FIRST OUPUT TIME: Must open existing file and find out!
          NCF => NEW_FILE()
          NCF%FNAME=NC_AVG%FNAME
          
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)

          ! GET THE CURRENT STACK LENGTH
          stklen = NCF%FTIME%stk_len
          
          !!! DAS{ 5.28.14

          ! GET THE TIME OF THE FIRST OUTPUT
          IF (STKLEN .GT. 1) THEN

             CALL UPDATE_FILE_BRACKET(NCF,StartTime,STATUS)
             IF(status == 0) THEN

                !PREV_TIME = NCF%FTIME%PREV_IO
                !PREVSTK= NCF%FTIME%PREV_STKCNT
                !
                !! SET THE NEXT OUTPUT
                !Nextstk = NCF%FTIME%NEXT_STKCNT 
                !NEXT_TIME = NCF%FTIME%NEXT_IO
                
                ! The last output will be the previous
                !PREV_TIME = NCF%FTIME%NEXT_IO 
                !PREVSTK   = NCF%FTIME%NEXT_STKCNT 

                !Nextstk   = NCF%FTIME%NEXT_STKCNT + 1
                !NEXT_TIME = NCF%FTIME%NEXT_IO + INTERVAL_TIME
                !!! END DAS

                IF(NCF%FTIME%NEXT_IO == StartTime) then
                   NEXT_TIME = StartTime
                   Nextstk = NCF%FTIME%NEXT_STKCNT + 1
                   
                   PREV_TIME = NCF%FTIME%PREV_IO
                   PREVSTK = NCF%FTIME%PREV_STKCNT + 1

                else if (NCF%FTIME%PREV_IO == StartTime) then
                   NEXT_TIME = StartTime
                   Nextstk = NCF%FTIME%NEXT_STKCNT
                   
                   PREV_TIME = StartTime - INTERVAL_TIME
                   PREVSTK   = NCF%FTIME%PREV_STKCNT

                else
                   call fatal_error('Something is very wrong with update_file_bracket!')

                end if

             ELSE

                call fatal_error('Start time is before or after average output ???')

                ! THERE IS NO SAVED OUTPUT PAST THE RESART TIME
                !NEXT_TIME = GET_FILE_TIME_NCF(NCF,stklen)
                !NEXTSTK= STKLEN

                !prevstk = stklen - 1
                !PREV_TIME = NEXT_TIME - INTERVAL_TIME
                !
                !IF(StartTime > PREV_TIME) call Fatal_error&
                !     &("COULD NOT FIND VALID BRACKET FOR START TIME IN THE DATA FILE",&
                !     & "Try comparing the time in the restart and the data file?")

             END IF

          ELSE IF(STKLEN .eq. 1)THEN
              
             ! CHECK TO SEE WHETHER THIS OUTPUT IS AHEAD OF THE START TIME
             !PREV_TIME = GET_FILE_TIME_NCF(NCF,stklen)
             !IF (PREV_TIME <= StartTime) THEN
             !   PREVSTK= STKLEN
             !   
             !   ! SET THE NEXT OUTPUR
             !   Nextstk = stklen + 1
             !   NEXT_TIME = PREV_TIME + INTERVAL_TIME
             !ELSE
             !   NEXTSTK=STKLEN
             !   NEXT_TIME= PREV_TIME
             !   PREV_TIME= ZEROTIME
             !END IF

             NEXT_TIME = GET_FILE_TIME_NCF(NCF,stklen)
             if (NEXT_TIME == StartTime) then
                nextstk = stklen + 1

                prev_time = zerotime
                prevstk = 0


             else

                call fatal_error('Cant crash restart when the restart time is not an average output time!')

             end if



          ELSE
             ! IF THERE ARE NOT YET ANY TIMESTEPS IN THE FILE
             PREV_TIME=ZEROTIME
             PREVSTK=0
             NEXTSTK=1
             CALL GET_FIRST_OUTPUT_TIME(TRIM(NCAV_FIRST_OUT),NEXT_TIME)
          END IF

          !DAS} 5.28.14


          IF(NEXT_TIME < StartTime) THEN 
             CALL PRINT_TIME(StartTime,IPT,"Start Time")
             CALL PRINT_TIME(NEXT_TIME,IPT,"First AVG OUTPUT")

             CALL FATAL_ERROR&
               & ("CAN NOT CONNECT TO EXISTING FILE - THERE IS A TIME GAP",&
               "BETWEEN THE EXISTING OUTPUT AND THE START TIME?")

          END IF

          
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first averaged output &
               &must be less than or equal to the end time")


          ! FINISHED WITH FILE
          CALL KILL_FILE(NCF)
          

          ! SET THE STACK MAX
          STKMAX = NCAV_OUTPUT_STACK


          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TIME AVERAGE OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_AVG,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)


           Call Print_File(NC_AVG)
           
       END IF

       IF (NC_ON) THEN
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(NC_OUT_INTERVAL),INTERVAL_TIME)

          ! FIRST OUPUT TIME: Must open existing file and find out!
          NCF => NEW_FILE()
          NCF%FNAME=NC_DAT%FNAME
          
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)

          ! GET THE CURRENT STACK LENGTH
          stklen = NCF%FTIME%stk_len

          ! GET THE TIME OF THE FIRST OUTPUT
          IF (STKLEN .GT. 1) THEN

             CALL UPDATE_FILE_BRACKET(NCF,StartTime,STATUS)
             IF(status == 0) THEN
                PREV_TIME = NCF%FTIME%PREV_IO
                PREVSTK= NCF%FTIME%PREV_STKCNT
                
                ! SET THE NEXT OUTPUT
                Nextstk = NCF%FTIME%NEXT_STKCNT 
                NEXT_TIME = NCF%FTIME%NEXT_IO
                
             ELSE

                ! THERE IS NO SAVED OUTPUT PAST THE RESART TIME
                PREV_TIME = GET_FILE_TIME_NCF(NCF,stklen)
                PREVSTK= STKLEN

                Nextstk = stklen + 1
                NEXT_TIME = PREV_TIME + INTERVAL_TIME

                IF(StartTime > NEXT_TIME) call Fatal_error&
                     &("COULD NOT FIND VALID BRACKET FOR START TIME IN THE DATA FILE",&
                     & "Try comparing the time in the restart and the data file?")

             END IF

          ELSE IF(STKLEN .eq. 1) THEN
             
             PREV_TIME = GET_FILE_TIME_NCF(NCF,stklen)
             IF (PREV_TIME <= StartTime) THEN
                PREVSTK= STKLEN
                
                ! SET THE NEXT OUTPUR
                Nextstk = stklen + 1
                NEXT_TIME = PREV_TIME + INTERVAL_TIME
             ELSE
                NEXTSTK=STKLEN
                NEXT_TIME= PREV_TIME
                PREV_TIME= ZEROTIME
             END IF

          ELSE
             ! IF THERE ARE NOT YET ANY TIMESTEPS IN THE FILE
             PREV_TIME=ZEROTIME
             PREVSTK=0
             NEXTSTK=1
             CALL GET_FIRST_OUTPUT_TIME(TRIM(NC_FIRST_OUT),NEXT_TIME)
          END IF
         
          IF(NEXT_TIME < StartTime) THEN 
             CALL PRINT_TIME(StartTime,IPT,"Start Time")
             CALL PRINT_TIME(NEXT_TIME,IPT,"First DATA OUTPUT")

             CALL FATAL_ERROR&
               & ("CAN NOT CONNECT TO EXISTING FILE - THERE IS A TIME GAP",&
               "BETWEEN THE EXISTING OUTPUT AND THE START TIME?")

          END IF

          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first file output &
               &must be less than or equal to the end time")

          
          ! SET THE STACK MAX
          STKMAX = NC_OUTPUT_STACK

          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NETCDF OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_DAT,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)
          
          ! FINISHED WITH FILE
          CALL KILL_FILE(NCF)

       END IF

       IF (NCSF_ON) THEN
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(NCSF_OUT_INTERVAL),INTERVAL_TIME)

          ! FIRST OUPUT TIME: Must open existing file and find out!
          NCF => NEW_FILE()
          NCF%FNAME=NC_SF%FNAME
          
          Call NC_OPEN(NCF)
          CALL NC_LOAD(NCF)

          ! GET THE CURRENT STACK LENGTH
          stklen = NCF%FTIME%stk_len

          ! GET THE TIME OF THE FIRST OUTPUT
          IF (STKLEN .GT. 1) THEN

             CALL UPDATE_FILE_BRACKET(NCF,StartTime,STATUS)
             IF(status == 0) THEN
                PREV_TIME = NCF%FTIME%PREV_IO
                PREVSTK= NCF%FTIME%PREV_STKCNT
                
                ! SET THE NEXT OUTPUT
                Nextstk = NCF%FTIME%NEXT_STKCNT 
                NEXT_TIME = NCF%FTIME%NEXT_IO
                
             ELSE

                ! THERE IS NO SAVED OUTPUT PAST THE RESART TIME
                PREV_TIME = GET_FILE_TIME_NCF(NCF,stklen)
                PREVSTK= STKLEN

                Nextstk = stklen + 1
                NEXT_TIME = PREV_TIME + INTERVAL_TIME

                IF(StartTime > NEXT_TIME) call Fatal_error&
                     &("COULD NOT FIND VALID BRACKET FOR START TIME IN THE DATA FILE",&
                     & "Try comparing the time in the restart and the data file?")

             END IF

          ELSE IF(STKLEN .eq. 1) THEN
             
             PREV_TIME = GET_FILE_TIME_NCF(NCF,stklen)
             IF (PREV_TIME <= StartTime) THEN
                PREVSTK= STKLEN
                
                ! SET THE NEXT OUTPUR
                Nextstk = stklen + 1
                NEXT_TIME = PREV_TIME + INTERVAL_TIME
             ELSE
                NEXTSTK=STKLEN
                NEXT_TIME= PREV_TIME
                PREV_TIME= ZEROTIME
             END IF

          ELSE
             ! IF THERE ARE NOT YET ANY TIMESTEPS IN THE FILE
             PREV_TIME=ZEROTIME
             PREVSTK=0
             NEXTSTK=1
             CALL GET_FIRST_OUTPUT_TIME(TRIM(NCSF_FIRST_OUT),NEXT_TIME)
          END IF
         
          IF(NEXT_TIME < StartTime) THEN 
             CALL PRINT_TIME(StartTime,IPT,"Start Time")
             CALL PRINT_TIME(NEXT_TIME,IPT,"First SURFACE DATA OUTPUT")

             CALL FATAL_ERROR&
               & ("CAN NOT CONNECT TO EXISTING FILE - THERE IS A TIME GAP",&
               "BETWEEN THE EXISTING OUTPUT AND THE START TIME?")

          END IF

          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first surface file output &
               &must be less than or equal to the end time")

          
          ! SET THE STACK MAX
          STKMAX = NCSF_OUTPUT_STACK

          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! SURFACE NETCDF OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_SF,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)
          
          ! FINISHED WITH FILE
          CALL KILL_FILE(NCF)

       END IF

       IF (RST_ON) THEN
       
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(RST_OUT_INTERVAL),INTERVAL_TIME)

          ! GET THE CURRENT STACK LENGTH
          stklen = NC_START%FTIME%stk_len

          NCF => NC_START

          ! GET THE TIME OF THE FIRST OUTPUT
          IF (STKLEN .GT. 1) THEN

             ! RESET STACK AND TIME AND LET BRACKET FIND IT
             NCF%FTIME%PREV_IO = ZEROTIME
             NCF%FTIME%prev_stkcnt = 0
             CALL UPDATE_FILE_BRACKET(NCF,StartTime,STATUS)
             IF(status /= 0) call Fatal_error&
                  &(" COULD NOT FIND VALID BRACKET FOR START TIME IN THE DATA FILE")

             ! SPECIAL FOR THE RESTART FILE -
             ! MAKE SURE THE PREVIOUS TIME IS THE START TIME FOR
             ! READING HOTSTART
             IF (NCF%FTIME%PREV_IO == StartTime) THEN

                PREV_TIME = NCF%FTIME%PREV_IO
                PREVSTK= NCF%FTIME%PREV_STKCNT
                
                ! SET THE NEXT OUTPUT
                Nextstk = NCF%FTIME%NEXT_STKCNT
                NEXT_TIME = NCF%FTIME%NEXT_IO
             ELSE IF (NCF%FTIME%NEXT_IO == StartTime) THEN

                PREV_TIME= NCF%FTIME%NEXT_IO
                PREVSTK= NCF%FTIME%NEXT_STKCNT
                
                NEXT_TIME= PREV_TIME + INTERVAL_TIME
                Nextstk = PREVSTK + 1
                
             ELSE
                CALL FATAL_ERROR("DID NOT MATCH START TIME TO RESTART FILE")
             END IF


          ELSE IF(STKLEN .eq. 1) THEN
                        
             PREV_TIME = GET_FILE_TIME_NCF(NCF,stklen)
             PREVSTK= STKLEN
             
             ! SET THE NEXT OUTPUT
             Nextstk = stklen + 1
             NEXT_TIME = PREV_TIME + INTERVAL_TIME
             

          ELSE

             CALL FATAL_ERROR("What a mess this is... ",&
                  &"I suggest a strong cocktail before you try and figure this one out!")

          END IF


          ! SET BACK TO CORRECT TIME FOR STARTUP
          NCF%FTIME%PREV_IO= PREV_TIME
          NCF%FTIME%PREV_STKCNT=PREVSTK
          
          NCF%FTIME%NEXT_IO = NEXT_TIME
          NCF%FTIME%NEXT_STKCNT = Nextstk
                   
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first file output &
               &must be less than or equal to the end time")

          ! SET THE STACK MAX
          STKMAX = RST_OUTPUT_STACK

          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! RESTART OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_RST,INTERVAL_TIME,NEXT_TIME&
               &,NextStk,PREV_TIME,PREVSTK,STKLEN,STKMAX)


       END IF

!=================================================
! COLDSTART 
    CASE("coldstart")
!=================================================
       if(DBG_SET(dbg_log)) then 
          WRITE(IPT,*)'!                SET TIME FOR COLDSTART                          !'
          WRITE(IPT,*)'!                                                                !'
       end if


       SELECT CASE(USE_REAL_WORLD_TIME)
       CASE(.TRUE.)
          
          ! GET THE START TIME
          StartTime = READ_DATETIME(START_DATE,DATE_FORMAT,TIMEZONE,status)
          if (status == 0 ) &
               & Call Fatal_Error("Could not read the date string START_DATE: "//trim(START_DATE))

          ! GET THE END TIME
          EndTime = READ_DATETIME(END_DATE,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the date string END_DATE: "&
               &//trim(END_DATE))

          ! SANITY CHECK
          if(StartTime .GT. EndTime) &
               & Call Fatal_Error("Runfile Start_Date exceeds or equal to End_Date")
          
          ! GET THE REFERENCE DATE TIME
          IF(DATE_REFERENCE /= 'default')THEN
            ReferenceDate = READ_DATETIME(DATE_REFERENCE,DATE_FORMAT,TIMEZONE,status)
            if (status == 0 ) &
               & Call Fatal_Error("Could not read the date string DATE_REFERENCE: "//trim(DATE_REFERENCE))
	  ELSE
	    ReferenceDate%MJD   = 0
	    ReferenceDate%MuSOD = 0
	  END IF  

          !CALCULATE THE NUMBER OF STEPS AND IEND
          IINT = 0
          ISTART = 1
          NSTEPS = CALCULATE_NUMBER_OF_TIMESTEPS(StartTime,EndTime)
          IEND = ISTART + NSTEPS

          IF(ASSOCIATED(NC_START)) THEN
             IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
                  & "! SETTING PREV_STCKNT FOR NETCDF COLD START FILE 1"

             CALL SET_STARTUP_FILE_STACK(StartTime)
             
          ELSE
             IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
                  & "! NO NETCDF COLD START FILE"

          END IF

          
          
       CASE(.FALSE.) ! THIS MODEL IS USING IDEALIZED TIME

           ! GET THE START AND END INFORMATION
          CALL IDEAL_TIME_STRING2TIME(START_DATE,BFLAG,StartTIME,IINT)          
          CALL IDEAL_TIME_STRING2TIME(END_DATE,EFLAG,EndTIME,IEND)
          
          ! SANITY CHECK
          IF (BFLAG /= EFLAG) CALL FATAL_ERROR&
               ('IDEALIZED MODEL TIME SPECIFICATION IS INCORRENT',&
               &'BEGIN AND END CAN BE IN EITHER CYCLES OR TIME BUT NOT MIXED',&
               & trim(start_date),trim(end_date) )

          IF (BFLAG == 'time') THEN

             !CALCULATE THE NUMBER OF STEPS             
             IINT = 0
             ISTART = 1
             NSTEPS = CALCULATE_NUMBER_OF_TIMESTEPS(StartTime,EndTime)
             IEND = ISTART + NSTEPS

             IF(ASSOCIATED(NC_START)) THEN
                IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
                     & "! SETTING PREV_STCKNT FOR NETCDF COLD START FI&
                     &LE 2"

             CALL SET_STARTUP_FILE_STACK(StartTime)

             ELSE
                IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
                     & "! NO NETCDF COLD START FILE"
                
             END IF

             
          ELSE IF(BFLAG == 'step') THEN

             ! CALCULATE NSTEPS
             NSTEPS = IEND - IINT 

             ! SANITY CHECK
             IF(NSTEPS .LT. 0) CALL FATAL_ERROR&
                  &('Number of steps can not be less than zero')

             ! CALCULATE THE START TIME
             StartTime = IMDTI * IINT
             
             ! CALCULATE THE END TIME
             EndTime = IMDTI * IEND
             
             ! ADVANCE THE TO THE FIRST TIME STEP OF THE MODEL FROM THE
             ! INITIAL CONDITION
             ISTART = IINT +1


             IF(ASSOCIATED(NC_START)) THEN
                IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
                     & "! SETTING PREV_STCKNT FOR NETCDF COLD START FI&
                     &LE 3"

             CALL SET_STARTUP_FILE_STACK(IINT)

             ELSE
                IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
                     & "! NO NETCDF COLD START FILE"
                
             END IF

          ELSE
             CALL FATAL_ERROR('IDEAL_TIME_STRING2TIME returned invalid flag')

          END IF

       END SELECT
       
       IntTime = StartTime
       ExtTime = StartTime
       RUNFILE_STARTTIME = STARTTIME

#      if defined (WAVE_CURRENT_INTERACTION)
       WaveTime = StartTime
#      endif  

       ! INITIALIZE IEXT
       IEXT = 1

       ! IF THERE IS AN INITIAL CONDITION FILE SET THE STCKNT
       IF(ASSOCIATED(NC_START)) THEN ! ONLY DO THIS IF THERE IS A STARTUP FILE
          CALL SET_STARTUP_FILE_STACK(StartTime)
       END IF

       CALL REPORT_TIME_SETUP

! SET THE OUT PUT TIME FOR THE DIFFERENT FILES
       IF (NCAV_ON) THEN
       
          CALL GET_FIRST_OUTPUT_TIME(TRIM(NCAV_FIRST_OUT),NEXT_TIME)
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(NCAV_OUT_INTERVAL),INTERVAL_TIME)

          NEXTSTK = 0
          PREVSTK = 0 
          PREV_TIME = ZEROTIME
          STKLEN = 0
          STKMAX = NCAV_OUTPUT_STACK


          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TIME AVERAGE OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_AVG,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)


          if (STARTTIME > NEXT_TIME)&
               & Call Fatal_Error("The time of the first averaged output &
               &must be greater than or equal to the start time")
          
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first averaged output &
               &must be less than or equal to the end time")

       END IF

       IF (NC_ON) THEN
       
          CALL GET_FIRST_OUTPUT_TIME(TRIM(NC_FIRST_OUT),NEXT_TIME)
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(NC_OUT_INTERVAL),INTERVAL_TIME)

          NEXTSTK = 0
          PREVSTK = 0 
          PREV_TIME = ZEROTIME
          STKLEN = 0
          STKMAX = NC_OUTPUT_STACK


          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NETCDF OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_DAT,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)



          if (STARTTIME > NEXT_TIME)&
               & Call Fatal_Error("The time of the first file output &
               &must be greater than or equal to the start time")
          
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first file output &
               &must be less than or equal to the end time")

       END IF

       IF (NCSF_ON) THEN
       
          CALL GET_FIRST_OUTPUT_TIME(TRIM(NCSF_FIRST_OUT),NEXT_TIME)
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(NCSF_OUT_INTERVAL),INTERVAL_TIME)

          NEXTSTK = 0
          PREVSTK = 0 
          PREV_TIME = ZEROTIME
          STKLEN = 0
          STKMAX = NCSF_OUTPUT_STACK


          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! SURFACE NETCDF OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_SF,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)



          if (STARTTIME > NEXT_TIME)&
               & Call Fatal_Error("The time of the first surface file output &
               &must be greater than or equal to the start time")
          
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first surface file output &
               &must be less than or equal to the end time")

       END IF

       IF (RST_ON) THEN
       
 
       
          CALL GET_FIRST_OUTPUT_TIME(TRIM(RST_FIRST_OUT),NEXT_TIME)
          
          CALL GET_OUTPUT_FILE_INTERVAL(TRIM(RST_OUT_INTERVAL),INTERVAL_TIME)

          NEXTSTK = 0
          PREVSTK = 0 
          PREV_TIME = ZEROTIME
          STKLEN = 0
          STKMAX = RST_OUTPUT_STACK


          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! RESTART OUTPUT:"
          CALL SET_OUTPUT_FILE_TIME(NC_RST,INTERVAL_TIME,NEXT_TIME&
               &,NEXTSTK,PREV_TIME,PREVSTK,STKLEN,STKMAX)



          if (STARTTIME > NEXT_TIME)&
               & Call Fatal_Error("The time of the first restart output &
               &must be greater than or equal to the start time")
          
          if (ENDTIME < NEXT_TIME)&
               & Call Fatal_Error("The time of the first restart output &
               &must be less than or equal to the end time")

       END IF

!===========================================
    CASE DEFAULT
!===========================================

       CALL FATAL_ERROR("UNKNOWN STARUP TYPE IN RUNFILE")

!===========================================
    END SELECT
!===========================================


!================================================
!
! MAKE SURE THAT OUTPUT INTERVALS ARE AN INTEGER 
! NUMBER OF TIME STEPS
!
!================================================

    IF(RST_ON) THEN
       INTERVAL_TIME = NC_RST%FTIME%INTERVAL
       NEXT_TIME = NC_RST%FTIME%NEXT_IO

       IF(MOD(INTERVAL_TIME,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("RST_OUT_INTERVAL must be an integer number of int&
               &ernal time steps!")
       END IF
       IF(MOD((NEXT_TIME - StartTime),IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("RST_FIRST_OUT must be an integer number of int&
               &ernal time steps from the StartTime!")
       END IF

    END IF

    IF(NC_ON)THEN
       INTERVAL_TIME = NC_DAT%FTIME%INTERVAL
       NEXT_TIME = NC_DAT%FTIME%NEXT_IO
       IF(MOD(INTERVAL_TIME,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("NC_OUT_INTERVAL must be an integer number of int&
               &ernal time steps!")
       END IF
       IF(MOD((NEXT_TIME - StartTime),IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("NC_FIRST_OUT must be an integer number of int&
               &ernal time steps from the StartTime!")
       END IF       
          

    END IF

    IF(NCSF_ON)THEN
       INTERVAL_TIME = NC_SF%FTIME%INTERVAL
       NEXT_TIME = NC_SF%FTIME%NEXT_IO
       IF(MOD(INTERVAL_TIME,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("NCSF_OUT_INTERVAL must be an integer number of int&
               &ernal time steps!")
       END IF
       IF(MOD((NEXT_TIME - StartTime),IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("NCSF_FIRST_OUT must be an integer number of int&
               &ernal time steps from the StartTime!")
       END IF       
          

    END IF

    IF(NCAV_ON)THEN
       INTERVAL_TIME = NC_AVG%FTIME%INTERVAL
       NEXT_TIME = NC_AVG%FTIME%NEXT_IO
       
       IF(MOD(INTERVAL_TIME,IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("NCAV_OUT_INTERVAL must be an integer number of int&
               &ernal time steps!")
       END IF
       IF(MOD((NEXT_TIME - StartTime),IMDTI) /= ZeroTime) THEN
          CALL FATAL_ERROR("NCAV_FIRST_OUT must be an integer number of int&
               &ernal time steps from the StartTime!")
       END IF


       ! CHECK TO MAKE SURE THAT THE RESTART TIMES MATCH THE AVERAGE
       ! OUTPUT TIMES
       IF(RST_ON) THEN
          IF(MOD(NC_RST%FTIME%INTERVAL,INTERVAL_TIME) /= ZeroTime) THEN
             CALL FATAL_ERROR("NCAV_OUT_INTERVAL: The restart file int&
                  &erval must be an integer number of Average intervals!")
       
          END IF

          IF(MOD(NC_RST%FTIME%NEXT_IO-NEXT_TIME,INTERVAL_TIME) /= ZeroTime) THEN
             CALL FATAL_ERROR&
                  &("NCAV_FIRST_OUT: Any Restart dump time must also be a average output time!",&
                  & "MOD(RST_FIRST_OUT - NCAV_FIRST_OUT,NCAV_OUT_INTERVAL)/=0 IS AN ERROR")
       
          END IF


       END IF

    END IF



!================================================
!
! SETUP RECALCULATION TIMES FOR RHO_MEAN
!
!================================================
    IF(RECALCULATE_RHO_MEAN) THEN
       RECALC_RHO_MEAN = STARTTIME

       CALL IDEAL_TIME_STRING2TIME(INTERVAL_RHO_MEAN,FLAG,DELT_RHO_MEAN,NSTEPS)
          
       IF(FLAG == 'step') THEN

          ! SANITY CHECK
          IF(NSTEPS .LE. 0) CALL FATAL_ERROR&
               &('Number of steps for INTERVAL_RHO_MEAN can not be less than zero')
          
          ! CALCULATE THE START TIME
          DELT_RHO_MEAN = IMDTI * NSTEPS
    
          ELSE IF(FLAG == 'time') THEN

          ! SANITY CHECK
          IF(DELT_RHO_MEAN .LE. zerotime) CALL FATAL_ERROR&
               &('INTERVAL_RHO_MEAN can not be LE zero seconds!')
          
       END IF

       IF (RST_ON) THEN

          INTERVAL_TIME = NC_RST%FTIME%INTERVAL
          NEXT_TIME = NC_RST%FTIME%NEXT_IO

          IF(MOD(INTERVAL_TIME,DELT_RHO_MEAN) /= ZeroTime) THEN
             CALL WARNING("MODEL RESTART IS NOT PERFECT WHEN USING RHO MEAN RECALCULATION,",&
                  "TO FIX: SET MOD(RST_OUT_INTERVAL,INTERVAL_RHO_MEAN)==0!")
          END IF


          IF(MOD((NEXT_TIME - StartTime),DELT_RHO_MEAN) /= ZeroTime) THEN
             CALL WARNING("MODEL RESTART IS NOT PERFECT WHEN USING RHO MEAN RECALCULATION,",&
                  "TO FIX: SET MOD(RST_FIRST_OUT-StartTime,INTERVAL_RHO_MEAN)==0!")
          END IF

       END IF

    END IF


  END SUBROUTINE SETUP_TIME
!=========================================================================
  SUBROUTINE REPORT_TIME_SETUP
    IMPLICIT NONE


    SELECT CASE(USE_REAL_WORLD_TIME)
    CASE(.TRUE.)

       if(DBG_SET(dbg_log)) then

          write(IPT,*) "! This case uses real world with specified start and end dates"


          call PRINT_REAL_TIME(StartTime,ipt,"Start Time",TIMEZONE)


          call PRINT_REAL_TIME(EndTime,ipt,"End Time",TIMEZONE)

          call PRINT_TIME(IMDTE,IPT,"External Time STEP")
          write(IPT,*) "! DTE(seconds)       = ",DTE
          call PRINT_TIME(IMDTI,IPT,"Internal Time STEP")
          write(IPT,*) "! DTI(seconds)       = ",DTI

          write(ipt,*) "!============================="
          write(ipt,*) "!  ISTART = ",ISTART
          write(ipt,*) "!  IEND   = ",IEND
          write(ipt,*) "!============================="
          write(ipt,*) "!+++++++++++ FINISED MODEL TIME SETUP ++++++++++++++"
          write(ipt,*) "!==================================================="
       end if
    CASE(.FALSE.)

       if(DBG_SET(dbg_log)) then

          write(IPT,*) "!This is an idealized case with a specified runtime"


          call PRINT_TIME(StartTime,ipt,"Start Time")


          call PRINT_TIME(EndTime,ipt,"End Time")
          
          call PRINT_TIME(IMDTE,IPT,"External Time STEP")
          write(IPT,*) "! DTE(seconds)       = ",DTE
          call PRINT_TIME(IMDTI,IPT,"Internal Time STEP")
          write(IPT,*) "! DTI(seconds)       = ",DTI


          write(ipt,*) "! ============================="
          write(ipt,*) "!   ISTART = ",ISTART
          write(ipt,*) "!   IEND   = ",IEND
          write(ipt,*) "! ============================="
          write(ipt,*) "! +++++++++++ FINISED MODEL TIME SETUP ++++++++++++++"
          write(ipt,*) "! ==================================================="
       end if
    END SELECT

  END SUBROUTINE REPORT_TIME_SETUP



  SUBROUTINE SET_MODEL_TIMESTEP
    IMPLICIT NONE
    if(DBG_SET(dbg_log)) &
         & write(IPT,*) "!========= Setting up model time parameters ==========="

# if defined(SEMI_IMPLICIT)

    IMDTI = SECONDS2TIME(INTSTEP_SECONDS)
    if(DBG_SET(dbg_IO)) write(IPT,*) "! IMDTE (microseconds)= ",IMDTE

    if (IMDTI .LE. ZEROTIME) &
         & Call Fatal_Error("INTSTEP_SECONDS must be greater than 10**-6 seconds!")

    ! SET DTI
    DTI = Seconds(IMDTI)

# else

    ! CONVERT EXTSTEP TO MICROSECONDS AND MAKE IT AN INT
    IMDTE = SECONDS2TIME(EXTSTEP_SECONDS)

    ! Get floating point DTE from IMDTE    
    DTE = SECONDS(IMDTE)

    if(DBG_SET(dbg_IO)) write(IPT,*) "! IMDTE (microseconds)= ",IMDTE
    if(DBG_SET(dbg_IO)) write(IPT,*) "! DTE (seconds)       = ",DTE

    if (IMDTE .LE. ZEROTIME) &
         & Call Fatal_Error("EXTSTEP_SECONDS must be greater than 10**-6 seconds!")

    if (ISPLIT .LE. 0) &
         & Call Fatal_Error("ISPLIT must be greater than zero!")

    ! Get DTI from DTE and ISPLIT
    IMDTI = ISPLIT * IMDTE
    DTI = Seconds(IMDTI)

# endif


  END SUBROUTINE SET_MODEL_TIMESTEP
  
  SUBROUTINE IDEAL_TIME_STRING2TIME(string,flag,ntime,tstep)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: STRING
    TYPE(TIME), INTENT(OUT) :: NTIME
    INTEGER(ITIME), INTENT(OUT) :: tstep
    CHARACTER(LEN=4), INTENT(OUT) :: FLAG
    
    INTEGER :: I

    INTEGER :: NLINE, NCHAR, INTVAL(150), NVAL
    REAL(DP) :: REALVAL(150)
    CHARACTER(LEN=40) :: VARNAME 
    CHARACTER(LEN=80) :: STRINGVAL(150)
    CHARACTER(LEN=7) :: VARTYPE
    LOGICAL :: LOGVAL


    NLINE = -1
    NCHAR = LEN_TRIM(STRING)
    
    CALL GET_VALUE(NLINE,NCHAR,STRING,VARNAME,VARTYPE,LOGVAL&
         &,STRINGVAL,REALVAL,INTVAL,NVAL)
    
    
    SELECT CASE(VARNAME)
    CASE('days')
       IF(VARTYPE == 'float')THEN
          ntime= days2time(realval(1))
          FLAG = "time"
       ELSE
          CALL FATAL_ERROR("IDEAL_TIME_STRING2TIME: Unknown Time or Length of Time:",&
               & TRIM(STRING), "Bad type (check for missing '.')")
       END IF
    CASE('seconds')
       
       IF(VARTYPE == 'float')THEN
          ntime = seconds2time(realval(1))
          FLAG = "time"

       ELSE
          CALL FATAL_ERROR("IDEAL_TIME_STRING2TIME: Unknown Time or Length of Time:",&
               & TRIM(STRING), "Bad type (check for missing '.')")

       END IF
       
    CASE('cycles')
       
       IF(VARTYPE == 'integer')THEN

          tstep = intval(1)
          FLAG = "step"

       ELSE
          CALL FATAL_ERROR("IDEAL_TIME_STRING2TIME: Unknown Time or Length of Time::",&
               & TRIM(STRING), "Bad type (remove '.' ?)")
       END IF
       
    CASE DEFAULT
       
          CALL FATAL_ERROR("IDEAL_TIME_STRING2TIME: Unknown Time or Length of Time::",&
               & TRIM(STRING), "Bad units, can be 'seconds' 'days' or 'cycles'")
       
    END SELECT
    
    
  END SUBROUTINE IDEAL_TIME_STRING2TIME


  FUNCTION CALCULATE_NUMBER_OF_TIMESTEPS(STIME,ETIME) RESULT(TSTEP)
    IMPLICIT NONE
    INTEGER(itime)          :: TSTEP
    TYPE(TIME), INTENT(IN)  :: STIME, ETIME
    real(DP)                :: temp

    
    ! CALCULAT THE NUMBER OF TIME STEPS
    temp = SECONDS(ETIME-STIME) /  SECONDS(IMDTI)
    TSTEP = Ceiling(temp)-1 

  END FUNCTION CALCULATE_NUMBER_OF_TIMESTEPS

  SUBROUTINE SET_LAST_FILE_TIME(NCF,STIME, STEP)
    IMPLICIT NONE
    TYPE(NCFILE), POINTER   :: NCF
    TYPE(TIME), INTENT(OUT)  :: Stime
    INTEGER(ITIME), INTENT(OUT) :: STEP

    INTEGER, TARGET :: FSTEP
    TYPE(TIME)              :: Atime
    integer                 :: status, I
    real(DP)                :: temp
    TYPE(NCVAR),  POINTER   :: VAR
    logical found
    
    IF(.NOT. ASSOCIATED(NCF)) CALL FATAL_ERROR &
         & ("SET_LAST_FILE_TIME: The file object is not assocaited!")
  

    Status = SET_FILE_TIME_TYPE(NCF)
    IF(STATUS /= 0) CALL FATAL_ERROR &
         & ("COULD NOT FIND A VALID TIME VARIABLE IN THE RESTART FILE: ",&
         & TRIM(NCF%FNAME))
    
    
    I = NCF%FTIME%STK_LEN
    
    IF (I .LE. 0) CALL FATAL_ERROR("FILE LENGTH IS LESS THAN ONE - NO VALID RESTART TIMES!")

    STIME  =GET_FILE_TIME_NCF(NCF,I)
    
    NCF%FTIME%PREV_STKCNT = I
    NCF%FTIME%PREV_IO     = STIME
    
    VAR => FIND_VAR(NCF,'iint',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'iint'&
         & IN THE CRASH RESTART FILE OBJECT")
    VAR%scl_int => FSTEP
    CALL NC_READ_VAR(VAR,I)
    
    STEP = FSTEP
    
  END SUBROUTINE SET_LAST_FILE_TIME

  SUBROUTINE SET_STARTUP_FILE_STACK_BY_TIME(STIME,STEP)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN)  :: Stime
    INTEGER(ITIME), INTENT(OUT),OPTIONAL :: STEP

    INTEGER, TARGET :: FSTEP
    TYPE(TIME)              :: Atime
    integer                 :: status, I
    real(DP)                :: temp
    TYPE(NCFILE), POINTER   :: NCF
    TYPE(NCVAR),  POINTER   :: VAR
    logical found
    

     IF(.NOT. ASSOCIATED(NC_START)) CALL FATAL_ERROR &
       & ("The file object NC_START is not assocaited!")
  

    Status = SET_FILE_TIME_TYPE(NC_START)
    IF(STATUS /= 0) CALL FATAL_ERROR &
         & ("COULD NOT FIND A VALID TIME VARIABLE IN THE RESTART FILE: ",&
         & TRIM(NC_START%FNAME))

    I = 0
    DO 
       I = I + 1
       IF (I .GT. NC_START%FTIME%STK_LEN) THEN
          
          IF(DBG_SET(DBG_LOG))&
               & CALL print_real_time(STime,IPT,"Asked for Start Time")
          
          IF(I == 2) THEN
            IF(DBG_SET(DBG_LOG))&
                 & CALL print_real_time(GET_FILE_TIME_NCF(NC_START,1),IPT,"The Only Restart Time Is")
             
          ELSE IF(I==1)THEN
             CALL FATAL_ERROR &
                  & ("Restart file has time dimension equal zero!:",&
                  & TRIM(NC_START%FNAME))
          ELSE
             ATIME  =GET_FILE_TIME_NCF(NC_START,1)
             IF(DBG_SET(DBG_LOG))&
                  & CALL print_real_time(ATIME,IPT,"First Restart Time")
             
             ATIME  =GET_FILE_TIME_NCF(NC_START,I-1)
             IF(DBG_SET(DBG_LOG))&
                  & CALL print_real_time(ATIME,IPT,"Last Restart Time")
          END IF
          
          CALL FATAL_ERROR &
               & ("COULD NOT FIND A MATCHING START TIME IN THE RESTART FILE!:",&
               & "(See time options printed above?)",&
               & TRIM(NC_START%FNAME))
       END IF
          
       ATIME = GET_FILE_TIME_NCF(NC_START,I)
       
       
       !IF (ATIME == STime) THEN
       IF(abs(ATIME -STime)<0.1_SP*IMDTI) THEN
        
          NC_START%FTIME%PREV_STKCNT = I
          !NC_START%FTIME%PREV_IO     = ATIME

          NC_START%FTIME%PREV_IO     = STime

          IF (PRESENT(STEP)) THEN
             VAR => FIND_VAR(NC_START,'iint',FOUND)
             IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'iint'&
                  & IN THE HOTSTART FILE OBJECT")
             VAR%scl_int => FSTEP
             CALL NC_READ_VAR(VAR,I)
             
             STEP = FSTEP
          END IF

          return
       END IF
       
    END DO
    
 
  END SUBROUTINE SET_STARTUP_FILE_STACK_BY_TIME

  SUBROUTINE SET_STARTUP_FILE_STACK_BY_CYCLE(step,ATIME)
    IMPLICIT NONE
    INTEGER(ITIME), INTENT(IN)  :: STEP
    TYPE(TIME), OPTIONAL, INTENT(OUT)     :: Atime
    integer                 :: status, I
    real(DP)                :: temp
    TYPE(NCFILE), POINTER   :: NCF
    TYPE(NCVAR),  POINTER   :: VAR
    INTEGER, TARGET :: FSTEP
    logical found

     IF(.NOT. ASSOCIATED(NC_START)) CALL FATAL_ERROR &
       & ("The file object NC_START is not assocaited!")


    Status = SET_FILE_TIME_TYPE(NC_START)
    IF(STATUS /= 0) CALL FATAL_ERROR &
         & ("COULD NOT FIND A VALID TIME VARIABLE IN THE RESTART FILE: ",&
         & TRIM(NC_START%FNAME))

    I = 0
    DO 
       I = I + 1
       IF (I .GT. NC_START%FTIME%STK_LEN) THEN

          CALL FATAL_ERROR &
               & ("COULD NOT FIND A MATCHING IINT CYCLE NUMBER IN THE RESTART FILE!:",&
               & TRIM(NC_START%FNAME))
       END IF

       VAR => FIND_VAR(NC_START,'iint',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'iint'&
            & IN THE HOTSTART FILE OBJECT")
       VAR%scl_int => FSTEP
       CALL NC_READ_VAR(VAR,I)
       
       IF (FSTEP == STEP) THEN
          
          if (PRESENT(ATIME)) THEN
             ATIME = GET_FILE_TIME_NCF(NC_START,I)
             
             NC_START%FTIME%PREV_STKCNT = I
             NC_START%FTIME%PREV_IO     = ATIME
          END if
          
          return
       END IF
       
    END DO
    
    
  END SUBROUTINE SET_STARTUP_FILE_STACK_BY_CYCLE
 
  SUBROUTINE CHECK_STARTUP_FILE_DIMENSIONS
    IMPLICIT NONE
    integer                 :: status, I
    real(DP)                :: temp
    TYPE(TIME)              :: Atime
    TYPE(NCFILE), POINTER   :: NCF
    TYPE(NCDIM), POINTER   :: DIM
    logical found
    

!!!!!
!!!! I DON'T THINK WE REALLY NEED TO IMPLIMENT THIS - JUST LET IT
!!!! CRASH LATER!
!!!!!

    IF(.NOT. ASSOCIATED(NC_START)) CALL FATAL_ERROR &
         & ("The file object NC_START is not assocaited!")
    


    IF (STARTUP_TYPE .EQ. STARTUP_TYPE_CRASHRESTART .OR. &
         & STARTUP_TYPE .EQ. STARTUP_TYPE_HOTSTART) THEN
    
       DIM => FIND_DIM(NC_START,"nele",found)
       IF(.not. FOUND) CALL FATAL_ERROR&
            & ("START FILE IS MISSING A CRITICAL DIMENSION:",&
            &  TRIM(NC_START%FNAME))
       
       IF (DIM%DIM /= NGL) CALL FATAL_ERROR&
            & ("START FILE DIMENSION 'nele' DOES NOT MATCH:",&
            &  TRIM(NC_START%FNAME))
       
       DIM => FIND_DIM(NC_START,"node",found)
       IF(.not. FOUND) CALL FATAL_ERROR&
            & ("START FILE IS MISSING A CRITICAL DIMENSION:",&
            &  TRIM(NC_START%FNAME))
       
       IF (DIM%DIM /= mgl) CALL FATAL_ERROR&
            & ("START FILE DIMENSION 'node' DOES NOT MATCH THE RUN FILE:",&
            &  TRIM(NC_START%FNAME))
    END IF

    IF(STARTUP_UV_TYPE .eq. STARTUP_TYPE_SETVALUES .OR. &
         & STARTUP_TURB_TYPE .eq. STARTUP_TYPE_SETVALUES .OR.&
         & STARTUP_TS_TYPE .eq. STARTUP_TYPE_SETVALUES) THEN
       
       DIM => FIND_DIM(NC_START,"siglay",found)
       IF(.not. FOUND) CALL FATAL_ERROR&
            & ("RESTART FILE IS MISSING A CRITICAL DIMENSION:",&
            &  TRIM(NC_START%FNAME))
       
       IF (DIM%DIM /= kbm1) CALL FATAL_ERROR&
            & ("RESTART FILE DIMENSION 'siglay' DOES NOT MATCH THE RUN FILE:",&
            &  TRIM(NC_START%FNAME))
       
       DIM => FIND_DIM(NC_START,"siglev",found)
       IF(.not. FOUND) CALL FATAL_ERROR&
            & ("RESTART FILE IS MISSING A CRITICAL DIMENSION:",&
            &  TRIM(NC_START%FNAME))
       
       IF (DIM%DIM /= kbm1) CALL FATAL_ERROR&
            & ("RESTART FILE DIMENSION 'siglay' DOES NOT MATCH THE RUN FILE:",&
            &  TRIM(NC_START%FNAME))
    END IF

  END SUBROUTINE CHECK_STARTUP_FILE_DIMENSIONS

  SUBROUTINE GET_FIRST_OUTPUT_TIME(STRING,FIRST_TIME)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: STRING
    TYPE(TIME), INTENT(OUT) :: FIRST_TIME
    INTEGER ::STATUS

    if (USE_REAL_WORLD_TIME) then
       FIRST_Time = READ_DATETIME(STRING,DATE_FORMAT,TIMEZONE,status)
       if (status == 0 ) Call Fatal_Error &
            ("GET_FIRST_OUTPUT_TIME: Could not read the date string: "//TRIM(STRING))
    else
       ! THIS IS REALLY THE SAME AS WHAT IS NEEDED IN THIS CASE
       ! REALLY WE ARE GETTING THE START TIME NOT THE INTERVAL
       CALL GET_OUTPUT_FILE_INTERVAL(TRIM(STRING),FIRST_TIME)
    end if
    
  END SUBROUTINE GET_FIRST_OUTPUT_TIME



  SUBROUTINE GET_OUTPUT_FILE_INTERVAL(STRING,INTERVAL)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: STRING
    TYPE(TIME), INTENT(OUT) :: INTERVAL

    CHARACTER(LEN=4) :: FLAG
    TYPE(TIME) :: NTIME
    INTEGER(ITIME) :: tstep
    INTEGER :: I

    CALL  IDEAL_TIME_STRING2TIME(string,flag,ntime,tstep)

    IF (FLAG == 'time') THEN ! IF START AND END TIME WERE SPECIFIED

       INTERVAL = ntime
    
    ELSE IF(FLAG == 'step') THEN ! IF START AND END IINT WERE SPECIFIED

       INTERVAL= tstep * IMDTI

    ELSE
       CALL FATAL_ERROR('GET_OUTPUT_FILE_INTERVAL: bad flag value?')
    END IF
    

  END SUBROUTINE GET_OUTPUT_FILE_INTERVAL


  SUBROUTINE SET_OUTPUT_FILE_TIME(NCF,OUT_INTERVAL,NEXT_TIME,NEXT_STK,PREV_TIME,PREV_STK,STKLEN,STACK_MAX)
    IMPLICIT NONE
    TYPE(NCFILE), POINTER      :: NCF
    TYPE(TIME), INTENT(IN)     :: NEXT_TIME, PREV_TIME
    TYPE(TIME), INTENT(IN) :: OUT_INTERVAL
    INTEGER, INTENT(IN)        :: STACK_MAX,STKLEN,NEXT_STK, PREV_STK

    TYPE(TIME)              :: ZEROTIME
    TYPE(NCFTIME), POINTER  :: FTM
    logical found
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"SET_OUTPUT_FILE_TIME: START"
    
    IF(DBG_SET(DBG_IO)) THEN
       WRITE(IPT,*)"=== PRINTING IO INFO FOR: SET_OUTPUT_FILE_TIME ==="
       CALL PRINT_TIME(OUT_INTERVAL,IPT,"OUT_INTERVAL")
       CALL PRINT_TIME(NEXT_TIME,IPT,"NEXT_TIME")
       WRITE(IPT,*)"NEXT_STK=",NEXT_STK
       CALL PRINT_TIME(PREV_TIME,IPT,"PREV_TIME")
       WRITE(IPT,*)"PREV_STK=",PREV_STK
       WRITE(IPT,*)"STKLEN=",STKLEN
       WRITE(IPT,*)"STACK_MAX=",STACK_MAX
       WRITE(IPT,*)"=== END IO INFO FOR: SET_OUTPUT_FILE_TIME ==="
    END IF

    IF (.NOT.ASSOCIATED(NCF)) CALL FATAL_ERROR &
         & ('SET_OUTPUT_FILE_TIME: FILE OBJECT HANDLE IS NOT ASSOCIATED')
       
    IF (.NOT.ASSOCIATED(NCF%FTIME)) NCF%FTIME=>NEW_FTIME()
    FTM => NCF%FTIME
    
    if (OUT_INTERVAL < ZEROTIME) Call Fatal_Error &
         & ("The output interval must be greater than or equal to zero",&
         & "Check the run file INTERVAL")
 
    ! SET CURRENT STACK LENGTH
    FTM%STK_LEN=STKLEN

    ! SET THE MAXIMUM STACK SIZE FOR THE FILE
    FTM%MAX_STKCNT = STACK_MAX
    
    ! SET THE OUPUT INTERVAL TIME
    FTM%INTERVAL = OUT_INTERVAL

    ! SET THE TIME OF THE NEXT OUTPUT
    FTM%NEXT_IO = NEXT_TIME

    !SET THE STK OF THE NEXT OUTPUT
    FTM%NEXT_STKCNT = NEXT_STK

    ! SET THE TIME OF THE PREV OUTPUT
    FTM%PREV_IO = PREV_TIME
    
    !SET THE STK OF THE PREV OUTPUT
    FTM%PREV_STKCNT = PREV_STK
    

    IF(DBG_SET(DBG_IO)) THEN
       WRITE(IPT,*) "! === DUMPING OUTPUT FILE TIMING INFO: ==="
       CALL PRINT_FTIME(FTM)
       WRITE(IPT,*) "! ========================================"

    ELSE IF (DBG_SET(DBG_LOG)) THEN
       
       IF (USE_REAL_WORLD_TIME) then
          CALL PRINT_REAL_TIME(NEXT_TIME,IPT,'First Output Time')
       ELSE
          CALL PRINT_TIME(NEXT_TIME,IPT,'First Output Time')
       END IF
       
       CALL PRINT_TIME(OUT_INTERVAL,IPT,"Output Interval")

    END IF
    
    nullify(FTM)
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"SET_OUTPUT_FILE_TIME: END"
    
  END SUBROUTINE SET_OUTPUT_FILE_TIME
  
END MODULE MOD_SET_TIME
