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

# ifdef VISIT
!!!====================================================================
!!!====================================================================
!!!====================================================================
!!!
!!!         BEGIN VISIT SUBROUTINES
!!!
!!!====================================================================
!!!====================================================================
!!!====================================================================

!NOTES:
! Variables not declared at the header of each subroutine come from
! the other modules in FVCOM.


SUBROUTINE INIT_VISIT
!  USE MOD_UTILS
  USE CONTROL, only : fvcom_version, case_title, PAR, WAIT_FOR_VISIT
  USE LIMS, only : MYID
  USE MOD_VISIT
  USE MOD_UTILS
  USE MOD_NCDIO, only : VISIT_CMD_DUMP
  implicit none
  include "visitfortransiminterface.inc"

  TYPE(VisitMeshType),pointer  :: VMS
! local variables
  integer err,ind
  character(len=100) :: SimName, timestamp
  integer :: LSimName, LTitle

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "BEGIN INIT_VISIT:"

  If(DBG_SET(DBG_LOG)) WRITE(IPT,*)"! INITILIZING VISIT"


  SLAVECALLBACKCOUNT=0
  BROADCASTSTRCOUNT=0
  BROADCASTINTCOUNT=0
 
  if (DBG_SET(DBG_SCL)) then
     WRITE(IPT,*)"=================================="
     WRITE(IPT,*)"       VISIT PARAMETERS"
     WRITE(IPT,*)"=================================="
     WRITE(IPT,*)"VISIT_TWODMESH=",VISIT_TWODMESH
     WRITE(IPT,*)"VISIT_SSHMESH=",VISIT_SSHMESH
     WRITE(IPT,*)"VISIT_BATHYMESH=",VISIT_BATHYMESH
     WRITE(IPT,*)"VISIT_LAYERMESH=",VISIT_LAYERMESH
     WRITE(IPT,*)"VISIT_LEVELMESH=",VISIT_LEVELMESH
     
     WRITE(IPT,*)"VISIT_COMMAND_PROCESS=",VISIT_COMMAND_PROCESS
     WRITE(IPT,*)"VISIT_COMMAND_SUCCESS=",VISIT_COMMAND_SUCCESS
     WRITE(IPT,*)"VISIT_COMMAND_FAILURE=",VISIT_COMMAND_FAILURE
     
  end if


  ! SET MESH VARIABLE TYPES
  VISIT_MESH(VISIT_TWODMESH)%MESH=  VISIT_TWODMESH

  VISIT_MESH(VISIT_SSHMESH)%MESH=  VISIT_SSHMESH

  VISIT_MESH(VISIT_BATHYMESH)%MESH=  VISIT_BATHYMESH

  VISIT_MESH(VISIT_LAYERMESH)%MESH=  VISIT_LAYERMESH

  VISIT_MESH(VISIT_LEVELMESH)%MESH=  VISIT_LEVELMESH


  if (DBG_SET(DBG_SCL)) then
     VMS =>  VISIT_MESH(VISIT_TWODMESH)
     if(associated(VMS))  WRITE(IPT,*)"associated"
     Call PRINT_MESH(VMS)
  end if



! INITIALIZE SOME VISIT VARIABLES

  ! SET THE VISIT RUN FLAG 
  VISITRUNFLAG=1
  IF(WAIT_FOR_VISIT) VisitRunFlag=0 

  !MYID is result of MPI_comm_rank PLUS one!
  VisitParRank=MYID-1

  VisitOneTStep=.FALSE.
  VisitHalt=.false.
  VisitAnimate=.FALSE.
  VISIT_CMD_DUMP=.FALSE.

  VisitStep=VISIT_STEP_INT

  err = visitsetupenv()
  if(dbg_set(dbg_scl)) WRITE(IPT,*) "visitsetupenv=",err


#  if defined (MULTIPROCESSOR)
  IF(PAR) THEN
     err = visitsetparallel(PAR)
     if(dbg_set(dbg_scl)) WRITE(IPT,*)"visitsetparallel(PAR)=",err,": PAR=",PAR
     err= visitsetparallelrank(VisitParRank)
     if(dbg_set(dbg_scl)) WRITE(IPT,*)"visitsetparallelrank=",err,"; Visi&
          &tParRank=",VisitParRank
  END IF
# endif


  if(MSR) then
!     call get_timestamp(timestamp)
     LTitle=Len_Trim(case_title)
     SimName= trim(fvcom_version) ! could include fvcom timestamp BUT
     ! NO COLONS!
     LSimName=Len_Trim(SimName)
     WRITE(IPT,*)trim(SimName),",",LSimName
     WRITE(IPT,*)trim(case_title),',',Ltitle

     err = visitinitializesim(trim(simname),lsimname, trim(case_title), ltitle,&
          & "/no/useful/path", 15,&
          & visit_f77nullstring, visit_f77nullstringlen,&
          & visit_f77nullstring, visit_f77nullstringlen,&
          & visit_f77nullstring, visit_f77nullstringlen)

     WRITE(IPT,*) "MASTER WROTE VISIT SIM FILE: ",trim(simname)
  end if



  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END INIT_VISIT"

  
end SUBROUTINE INIT_VISIT

!!====================================================================
!!====================================================================
!!
!! VISITCHECK decides when check the viewer for connections and commands
!!
!!====================================================================
!!====================================================================


SUBROUTINE VisitCheck
  USE MOD_VISIT
  USE CONTROL, only : IINT, IEXT,IntTime,ExtTime
  implicit none
  


# if defined (SEMI_IMPLICIT)

  visit_cycle=IINT-1
  visit_time_int=IntTime
  visit_time_ext=IntTime

  if (VisitStep .EQ. VISIT_STEP_INT) then
     
     Call VisitSimLoop
     
  elseif (VisitStep .EQ. VISIT_STEP_10XINT) then
     VisitStepCount=VisitStepCount+1
     if (mod(VisitStepCount,10) .EQ. 0) Call VisitSimLoop
               
  elseif (VisitStep .EQ. VISIT_STEP_100XINT) then
     VisitStepCount=VisitStepCount+1
     if (mod(VisitStepCount,100) .EQ. 0) Call VisitSimLoop
          
  end if
# else

  visit_cycle=IINT-1
  visit_time_int=IntTime
  visit_time_ext=ExtTime

  if (VisitStep .EQ. VISIT_STEP_EXT) then
     Call VisitSimLoop
     
  elseif (VisitStep .EQ. VISIT_STEP_INT) then
     
     if (IEXT .EQ. 1) Call VisitSimLoop
     
  elseif (VisitStep .EQ. VISIT_STEP_10XINT) then
     if(IEXT .EQ. 1) then
        VisitStepCount=VisitStepCount+1
        if (mod(VisitStepCount,10) .EQ. 0) Call VisitSimLoop
        
     end if
     
  elseif (VisitStep .EQ. VISIT_STEP_100XINT) then
     if(IEXT .EQ. 1) then
        VisitStepCount=VisitStepCount+1
        if (mod(VisitStepCount,100) .EQ. 0) Call VisitSimLoop
        
     end if
     
  end if
  
# endif

END SUBROUTINE VisitCheck



!!====================================================================
!!====================================================================
!!
!! VisitSimLoop
!!
!! This is the loop that does all the work of talking to visit
!!====================================================================
!!====================================================================


subroutine VisitSimLoop
  USE MOD_VISIT
# if defined (MULTIPROCESSOR)
  USE MPI
# endif
  implicit none
  include "visitfortransiminterface.inc"
  integer :: visitstate, vresult, blocking, loopcount,IERR, ind
  integer, external :: processvisitcommand
  integer, external :: visitbroadcastintfunction

  TYPE(VISITMESHTYPE), POINTER :: vmp ! VISIT MESH POINTER

  TYPE(VISITLAG), POINTER :: lmp ! VISIT LAGRANGIAN POINTER

  TYPE(VisitSphereVel), pointer :: VSV


! VisitRunFlag == 1       Simulation is Running

! VisitRunFlag == 0       Simulation is Stopped


if(dbg_set(dbg_vrb)) &
     &  WRITE(IPT,*)"BEGIN DoVisitSim: PID", VisitParRank



  IF (VISITANIMATE) then
     IF (dbg_set(dbg_log)) WRITE(IPT,*)"Updating Plots:"
     vresult = visittimestepchanged()
     vresult = visitupdateplots()
  END IF
  
  if (VisitOneTStep) then
     if(dbg_set(dbg_vrb)) &
          &  WRITE(IPT,*)"Hit VisitOneTStep at beginning of visitsimloop"

     VisitOneTStep=.false.
     VisitHalt=.true.
  end if

  DO 

     ! LOGIC TO DO ONE TIME STEP ONLY
     
     if(dbg_set(dbg_vrb)) THEN
        WRITE(IPT,*)"Start VISIT SIM LOOP COUNT=", loopcount
        WRITE(IPT,*)"visitRunFlag=",visitrunflag 
     end if
     
     if(VisitRunFlag.eq.1) then
        blocking = 0 ! Simulation running
     else !
        blocking = 1 ! Simulation stopped, waiting for commands
     endif
 
!     !Allow the one step call back to over ride visitrunflag
!     if (VisitOneTStep .eq. .TRUE.) then
!        blocking =0
!        VisitOneTStep=.False.
!        VisitRunFlag=0  ! Make sure FVCOM only goes one step!
!     end if

     if(dbg_set(dbg_vrb)) &
          & WRITE(IPT,*)"blocking=",blocking

     if(MSR) visitstate = visitdetectinput(blocking, -1)
     
     if(dbg_set(dbg_vrb)) then
        IF(MSR) WRITE(IPT,*)"Sending visitstate from master, VISITSTATE=",VISITSTATE
        IF(.not.MSR) WRITE(IPT,*)"Recieving visitstate from master"
     end if

     IERR= visitbroadcastintfunction(visitstate, 0)

if(dbg_set(dbg_vrb)) &
     & WRITE(IPT,*) "Sent: VISITSTATE= ",VISITSTATE


     if (visitstate.lt.0) then
        
        if(DBG_SET(DBG_LOG))then ! equivlent to VisitParRank=0
           WRITE(IPT,*)'VISITSTATE error codes: ',VISITSTATE
           WRITE(IPT,*)"-5: Logic error (fell through all cases)"
           WRITE(IPT,*)'-4: Logic error (no descriptors but blocking)'
           WRITE(IPT,*)'-3: Logic error (a socket was selected but not one we set)'
           WRITE(IPT,*)"-2: Unknown error in select"
           WRITE(IPT,*)"-1: Interrupted by EINTR in select"
        end if
        
        ! Stop the simulation!
        CALL PSTOP
        
     elseif (visitstate.eq.0) then

        if(dbg_set(dbg_vrb)) &
             & WRITE(IPT,*) "VisitRunFlag= ", VisitRunFlag,": PID=",VisitParRank

        if (VisitHalt) then
           if(dbg_set(dbg_vrb)) &
                & WRITE(IPT,*)"Hit VisitHalt while exiting loop with visitsta&
                &te == 0 "

           VisitRunFlag=0
           VisitHalt=.false.
        else
           exit ! return to US_FVCOM Main loop
        end if

        !!call simulate_one_timestep()
     elseif (visitstate.eq.1) then
        VisitRunFlag = 0
        vresult = visitattemptconnection()

        if (MSR) then
           if (vresult.eq.1) then
              if (dbg_set(dbg_log)) WRITE(IPT,*) 'VisIt connected!'
           else
              if (dbg_set(dbg_log)) WRITE(IPT,*)'VisIt did not connect!'
           endif
        end if


        if (vresult.eq.1) VisitAnimate=.TRUE.


     elseif (visitstate.eq.2) then

        if(dbg_set(dbg_vrb)) &
             & WRITE(IPT,*)'Run ProcessVisitCommand'

        IERR = ProcessVisitCommand()
        
        if(IERR.eq.0) then
           if(MSR) WRITE(IPT,*)'ProcessVisitCommand returned disconnect PID=', visitParRank
           vresult = visitdisconnect()
           VisitRunFlag = 1
           VisitAnimate=.FALSE.
           VisitOneTStep=.false.
           VisitHalt=.false.
           VISIT_CMD_DUMP=.FALSE.

!           VisitStep= VISIT_STEP_10XINT
           
           VMP => VISIT_MESH(VISIT_TWODMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           VMP => VISIT_MESH(VISIT_SSHMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           VMP => VISIT_MESH(VISIT_BATHYMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           VMP => VISIT_MESH(VISIT_LAYERMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           VMP => VISIT_MESH(VISIT_LEVELMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           NULLIFY(VMP)

           LMP => VISIT_LAGDATA
           Call DeAllocate_LAG(LMP,IERR)

           NULLIFY(LMP)

# if defined (SPHERICAL)
           VSV =>VisitSphericalVel
           Call DeAllocate_SphereVel(VSV,IERR)

           VSV =>VisitSphericalAVel
           Call DeAllocate_SphereVel(VSV,IERR)

           VSV =>VisitSphericalWindVel
           Call DeAllocate_SphereVel(VSV,IERR)

# if defined (ICE)
           VSV =>VisitSphericalIceVel
           Call DeAllocate_SphereVel(VSV,IERR)
# endif           

           NULLIFY(VSV)
# endif


           exit ! Return to simulation loop
        endif

     endif
     
     if(dbg_set(dbg_vrb)) &
          & WRITE(IPT,*)"END VISIT SIM LOOP COUNT=", loopcount


  end DO
  

  if(dbg_set(dbg_vrb)) &
       & WRITE(IPT,*)"END DoVisitSim"

end subroutine VisitSimLoop

!!=============================================================
!!=============================================================
!!=============================================================

! PARRALLEL HELPER COMMAND FOR DOVISITSIM
integer function ProcessVisitCommand()
  USE MOD_VISIT
  implicit none
  include "visitfortransiminterface.inc"
  integer :: IERR, processcount
  integer :: command, success
  
  if(dbg_set(dbg_sbr)) &
       & WRITE(IPT,*)"BEGIN ProcessVisitCommand"


  IERR=0

  if(MSR) then

     if(dbg_set(dbg_vrb)) &
     & WRITE(IPT,*)"Master Processing Engine command:"

     success =  visitprocessenginecommand()
     
     if(success.eq.1) then

        if(dbg_set(dbg_vrb)) &
             & WRITE(IPT,*)"Master Success Processing Engine command: BROADCAST2SLAVES"

        command = VISIT_COMMAND_SUCCESS
        IERR=1
     else
        if(dbg_set(dbg_vrb)) &
        & WRITE(IPT,*)"Master Failure Processing Engine command: BROADCAST2SLAVES"

        command = VISIT_COMMAND_FAILURE
        IERR=0
     end if
     
     call BroadcastSlaveCommand(command)
     
  else
     
     processcount=0
     DO

        if(dbg_set(dbg_vrb)) &
             & WRITE(IPT,*)"Slave asking for process engine command: processcount=&
             &",processcount

        processcount=processcount+1
        call BroadcastSlaveCommand(command)
        SELECT CASE (command)
           
        CASE (VISIT_COMMAND_PROCESS)

           if(dbg_set(dbg_vrb)) &
                & WRITE(IPT,*)"Slave: process engine command"


           success=visitprocessenginecommand()
           IERR=1

        CASE (VISIT_COMMAND_SUCCESS)
           if(dbg_set(dbg_vrb)) &
                & WRITE(IPT,*)"Slave recieved: visit command success"
           
           IERR =1
           exit
        CASE (VISIT_COMMAND_FAILURE)
           if(dbg_set(dbg_vrb)) &
           & WRITE(IPT,*)"Slave recieved: visit command failure"

           IERR=0
           exit
           
        CASE DEFAULT

           CALL FATAL_ERROR("Slave recieved: !Unrecognized visit command!")

        END SELECT
     END DO

  end if

  if(dbg_set(dbg_vrb)) &
       & WRITE(IPT,*)"END ProcessVisitCommand"

  ProcessVisitCommand = IERR
end function ProcessVisitCommand
!! END PARALLEL HELPER FOR DoVisitSim


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c These functions must be defined to satisfy the visitfortransiminterface lib.
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c---------------------------------------------------------------------------
!c visitcommandcallback
!c---------------------------------------------------------------------------
subroutine visitcommandcallback (cmd, lcmd, intdata, floatdata,&
     & stringdata, lstringdata)
  USE MOD_VISIT
  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8), intent(in) :: cmd, stringdata
  integer, intent(in) ::lcmd, lstringdata, intdata
  real, intent(in) :: floatdata
  integer :: result

  if(dbg_set(dbg_vrb)) &
       & WRITE(IPT,*) "BEGIN VISITCOMMANDCALLBACK: PID=",VisitParRank

  
  if(visitstrcmp(cmd, lcmd, "Halt/Step", 9).eq.0) then
     if (VisitRunFlag .EQ. 1) then
        ! VisitRunFlag = 1      ! if running, halt
        VisitHalt=.TRUE.     ! catch at next entry to visitsimloop

     else  ! visitRunFlag .eq. 0  - stopped
        visitRunFlag =1       ! if stopped go one time step
        VisitOneTStep=.TRUE.
     end if

  elseif(visitstrcmp(cmd, lcmd, "Run", 3).eq.0) then
     VisitRunFlag = 1

  elseif(visitstrcmp(cmd, lcmd, "Dump Restart", 12).eq.0) then
     VISIT_CMD_DUMP = .true.

  elseif(visitstrcmp(cmd, lcmd, "Toggle Step", 11).eq.0) then

# if !defined (SEMI_IMPLICIT)
     SELECT CASE(VisitStep)

     CASE(VISIT_STEP_100XINT)
        VisitStep=VISIT_STEP_EXT
        if(MSR) then
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"     VISIT TIME STEP CHANGED TO EVERY EXT"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
        end if
     CASE(VISIT_STEP_EXT)
        VisitStep=VISIT_STEP_INT
        if (MSR) then
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"     VISIT TIME STEP CHANGED TO ONE IINT"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
        end if
     CASE(VISIT_STEP_INT)
        VisitStep=VISIT_STEP_10XINT
        if(MSR) then
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"     VISIT TIME STEP CHANGED TO 10 X IINT"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
        end if
     CASE(VISIT_STEP_10XINT)
        VisitStep=VISIT_STEP_100XINT
        if(MSR) then
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"     VISIT TIME STEP CHANGED TO 100 X IINT"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
        end if
     CASE DEFAULT
        VisitStep=VISIT_STEP_INT
        if (MSR) then
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)" VISIT TIME STEP CHANGED TO ONE IINT: DEFAULT"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
        end if
     END SELECT
# else
     SELECT CASE(VisitStep)

     CASE(VISIT_STEP_100XINT)
        VisitStep=VISIT_STEP_INT
        if (MSR) then
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"     VISIT TIME STEP CHANGED TO ONE IINT"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
        end if
     CASE(VISIT_STEP_INT)
        VisitStep=VISIT_STEP_10XINT
        if(MSR) then
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"     VISIT TIME STEP CHANGED TO 10 X IINT"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
        end if
     CASE(VISIT_STEP_10XINT)
        VisitStep=VISIT_STEP_100XINT
        if(MSR) then
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"     VISIT TIME STEP CHANGED TO 100 X IINT"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
        end if
     CASE DEFAULT
        VisitStep=VISIT_STEP_INT
        if (MSR) then
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)" VISIT TIME STEP CHANGED TO ONE IINT: DEFAULT"
           WRITE(IPT,*)"================================================"
           WRITE(IPT,*)"================================================"
        end if
     END SELECT

# endif



!  elseif(visitstrcmp(cmd, lcmd, "ext", 3).eq.0) then
!     VisitStep=VISIT_STEP_EXT
!
!  elseif(visitstrcmp(cmd, lcmd, "int", 3).eq.0) then
!     VisitStep=VISIT_STEP_INT
!
!  elseif(visitstrcmp(cmd, lcmd, "10Xint", 6).eq.0) then
!     VisitStep=VISIT_STEP_10XINT
!     VisitStepCount=0
!
!  elseif(visitstrcmp(cmd, lcmd, "100Xint", 7).eq.0) then
!     VisitStep=VISIT_STEP_100XINT
!     VisitStepCount=0
!

  endif


  
  if(dbg_set(dbg_vrb)) &
       & WRITE(IPT,*) "END VISITCOMMANDCALLBACK: PID=",VisitParRank

end subroutine visitcommandcallback

!c---------------------------------------------------------------------------
!c visitbroadcastintfunction
!c---------------------------------------------------------------------------
integer function visitbroadcastintfunction(value, sender)
# if defined (MULTIPROCESSOR)
  USE MPI
  USE MOD_VISIT
# endif
  implicit none
  include "visitfortransiminterface.inc"
  integer, intent(inout) :: value, sender
  integer ::IERR

#  if defined (MULTIPROCESSOR)

  BROADCASTINTCOUNT=BROADCASTINTCOUNT+1
  if(dbg_set(dbg_vrb))&
       & WRITE(IPT,*)"BEGIN VISITBROADCASTINTFUNCTION: COUNT=", BROADCASTINTCOUNT

  CALL MPI_BCAST(value,1,MPI_INTEGER,sender,MPI_FVCOM_GROUP,IERR)

  if(dbg_set(dbg_vrb)) &
       & WRITE(IPT,*)"END VISITBROADCASTINTFUNCTION: COUNT=", BROADCASTINTCOUNT

  ! endif MPI
# endif


  visitbroadcastintfunction = 0


end function visitbroadcastintfunction

!c---------------------------------------------------------------------------
!c visitbroadcaststringfunction
!c---------------------------------------------------------------------------
integer function visitbroadcaststringfunction(str, lstr, sender)
# if defined (MULTIPROCESSOR)
  USE MPI
  USE MOD_VISIT
# endif
  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8), intent(inout):: str
  integer, intent(in):: sender, lstr
  integer :: IERR


# if defined (MULTIPROCESSOR)

  BROADCASTSTRCOUNT= BROADCASTSTRCOUNT +1
  if(dbg_set(dbg_vrb)) &
     & WRITE(IPT,*)"BEGIN VISITBROADCASTSTRINGFUNCTION: COUNT=", BROADCASTSTRCOUNT

  CALL MPI_BCAST(str,lstr,MPI_CHARACTER, sender ,MPI_FVCOM_GROUP,IERR)

  if(dbg_set(dbg_vrb)) &
       & WRITE(IPT,*)"END VISITBROADCASTSTRINGFUNCTION: COUNT=", BROADCASTSTRCOUNT

# endif

  visitbroadcaststringfunction = 0

end function visitbroadcaststringfunction

!c---------------------------------------------------------------------------
!c visitslaveprocesscallback
!c---------------------------------------------------------------------------
subroutine visitslaveprocesscallback()
  USE MOD_VISIT
  implicit none
  include "visitfortransiminterface.inc"
  integer :: command

  command = VISIT_COMMAND_PROCESS
 
  if(dbg_set(dbg_vrb)) &
       & WRITE(IPT,*) "Begin visitslaveprocesscallback: command=",command

  Call BroadcastSlaveCommand(command)
  
  if(dbg_set(dbg_vrb)) &
       WRITE(IPT,*) "End visitslaveprocesscallback"

end subroutine visitslaveprocesscallback


! HELPER FUNCTION FOR BROADCASTING COMMANDS FROM MASTER
subroutine BroadcastSlaveCommand(command)
# if defined (MULTIPROCESSOR)
  USE MPI
  USE MOD_VISIT
# endif
  implicit none
  integer, intent(in) :: command
  integer :: IERR
  
#  if defined (MULTIPROCESSOR)
  
  SLAVECALLBACKCOUNT=SLAVECALLBACKCOUNT+1

  if(dbg_set(dbg_vrb)) &
       & WRITE(IPT,*) "Broadcast from master to slaves:",command,": COUNT=", SLAVECALLBACKCOUNT

  CALL MPI_BCAST(command,1,MPI_INTEGER,0,MPI_FVCOM_GROUP,IERR)


if(dbg_set(dbg_vrb)) &
  WRITE(IPT,*) "Broadcast from master to slaves complete: COUNT=", SLAVECALLBACKCOUNT

# endif

end subroutine BroadcastSlaveCommand



!c---------------------------------------------------------------------------
!c visitgetmetadata
!c---------------------------------------------------------------------------
integer function visitgetmetadata(handle)
  USE ALL_VARS, only : NPROCS, KBM1, KBM2, KB, PAR, VISIT_ALL_VARS
  USE MOD_VISIT

  implicit none
  include "visitfortransiminterface.inc"
  integer, intent(in):: handle
  integer :: err, ind

  real(SP) :: vtime

  ! Mesh Variables
  integer :: tdim2, tdim3, sdim2,sdim3, mesh1,mesh2,mesh3,mesh4,mesh5&
       &, meshtype, meshlag
  character(LEN=20) :: UnitsName, LengthName, WidthName, DepthName,&
       & DomainName
  integer ::  LUnitsName, LLengthName, LWidthName, LDepthName, LDomainName 
  character(LEN=3) :: ch3

  character(LEN=20) :: matname
  character(LEN=20), allocatable :: matnamecopy(:)

  ! Scalar variables
  integer :: scalar

  ! Materials
  integer :: mat1, mat2
  ! Time stuff


  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! TIME AND CYCLE
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  if(dbg_set(dbg_sbrio)) then
     WRITE(IPT,*) "VisitSimCycle= ",Visit_Cycle,": PID=",VisitParRank
     WRITE(IPT,*) "VisitSimTime= ", Visit_Time_Ext,": PID=",VisitParRank
  end if
  vtime = seconds(Visit_Time_Ext)
  err= visitmdsetcycletime(handle, Visit_Cycle, vtime)  

  if(dbg_set(dbg_sbrio)) &
       & WRITE(IPT,*)"visitmdsetcycletime:result=",err

  if(VisitRunFlag.eq.0)then
     err = visitmdsetrunning(handle, VISIT_SIMMODE_RUNNING)
     !     WRITE(IPT,*)"visitmdsetrunning(running):result=",err
  else
     err = visitmdsetrunning(handle, VISIT_SIMMODE_STOPPED)
     !     WRITE(IPT,*)"visitmdsetrunning(stopped):result=",err
  end if
  

  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! MESH META DATA
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! General mesh variables
  meshtype = VISIT_MESHTYPE_UNSTRUCTURED
  tdim2 = 2
  tdim3 = 3

  sdim2 = 2
  sdim3 = 3

# if defined (SPHERICAL)


  UnitsName="Meters"
  LUnitsName=Len_Trim(UnitsName)

  LengthName="X"
  LLengthName=Len_Trim(LengthName)

  WidthName="Y"
  LWidthName=Len_Trim(WidthName)

  DepthName="Z"
  LDepthName=Len_Trim(DepthName)


# else
  
  UnitsName="Meters"
  LUnitsName=Len_Trim(UnitsName)
  
  LengthName="Zonal"
  LLengthName=Len_Trim(LengthName)
  
  WidthName="Meridional"
  LWidthName=Len_Trim(WidthName)
  
  DepthName="Depth"
  LDepthName=Len_Trim(DepthName)
  
# endif


  DomainName="Domain_"
  LDomainName=Len_Trim(DomainName)


# if defined (SPHERICAL)
  ! TWOD_Mesh
  mesh1 = visitmdmeshcreate(handle, "TWOD_Mesh", 9, meshtype, tdim2 &
       &,sdim3, NPROCS)
# else
  ! TWOD_Mesh
  mesh1 = visitmdmeshcreate(handle, "TWOD_Mesh", 9, meshtype, tdim2 &
       &,sdim2, NPROCS)
# endif

  if(mesh1.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh1, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh1, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName) ! Depth not used
     err = visitmdmeshsetblocktitle(handle, mesh1, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh1,DomainName , LDomainName)
  endif

  ! BATHYMETRY_Mesh
  mesh2 = visitmdmeshcreate(handle, "Bathymetry_Mesh", 15, meshtype, tdim2 &
       &,sdim3, NPROCS)

  if(mesh2.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh2, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh2, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName)
     err = visitmdmeshsetblocktitle(handle, mesh2, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh2,DomainName , LDomainName)
  endif

  ! SSH_Mesh
  mesh3 = visitmdmeshcreate(handle, "SSH_Mesh", 8, meshtype, tdim2 &
       &,sdim3, NPROCS)

  if(mesh3.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh3, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh3, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName)
     err = visitmdmeshsetblocktitle(handle, mesh3, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh3,DomainName , LDomainName)
  endif

  ! SIGMALAYER_Mesh
  mesh4 = visitmdmeshcreate(handle, "SigmaLayer_Mesh", 15, meshtype, tdim3 &
       &,sdim3, NPROCS)

  if(mesh4.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh4, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh4, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName)
     err = visitmdmeshsetblocktitle(handle, mesh4, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh4,DomainName , LDomainName)
  endif

  ! SIGMALEVEL_Mesh
  mesh5 = visitmdmeshcreate(handle, "SigmaLevel_Mesh", 15, meshtype, tdim3 &
       &,sdim3, NPROCS)

  if(mesh5.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh5, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh5, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName)
     err = visitmdmeshsetblocktitle(handle, mesh5, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh5,DomainName , LDomainName)
  endif


  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! ADD ONLINE LAGRANGIAN PARTICLE METADATA
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   if (LAG_PARTICLES_ON) then

      meshlag = visitmdmeshcreate(handle, "LAGRANGIAN_MESH", 15,&
           & VISIT_MESHTYPE_POINT,  0, 3, NPROCS)
      if(meshlag.ne.VISIT_INVALID_HANDLE) then
         err = visitmdmeshsetunits(handle, meshlag, UnitsName, LUnitsName)
         err = visitmdmeshsetlabels(handle, meshlag, LengthName,LLengthName&
              &, WidthName, LWidthName, DepthName,LDepthName)
         err = visitmdmeshsetblocktitle(handle, meshlag, "Domains", 7)
         err = visitmdmeshsetblockpiecename(handle, meshlag,DomainName , LDomainName)
      endif

      ! EASTWARD VELOCITY
      scalar = visitmdscalarcreate(handle, "LAG_U", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! NORTHWARD VELOCITY
      scalar = visitmdscalarcreate(handle, "LAG_V", 5, "LAGRANGIAN_MESH",15 , VISIT_VARCENTERING_NODE)

      ! UPWARD VELOCITY
      scalar = visitmdscalarcreate(handle, "LAG_W", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! Z LOCATION (METERS)
      scalar = visitmdscalarcreate(handle, "LAG_Z", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! SOME SCALAR
      scalar = visitmdscalarcreate(handle, "LAG_S1", 6, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! SOME SCALAR
      scalar = visitmdscalarcreate(handle, "LAG_S2", 6, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! SOME SCALAR
      scalar = visitmdscalarcreate(handle, "LAG_S3", 6, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! SOME SCALAR
      scalar = visitmdscalarcreate(handle, "LAG_S4", 6, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! SOME SCALAR
      scalar = visitmdscalarcreate(handle, "LAG_S5", 6, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! PATH LENGTH
      scalar = visitmdscalarcreate(handle, "LAG_D", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)


   end if

  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! ADD SCALAR METADATA
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  
  !  IF (VISIT_OPT .eq. 'basic') THEN
  scalar = visitmdscalarcreate(handle, "H", 1, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "EL", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
  
  scalar = visitmdscalarcreate(handle, "U", 1, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "V", 1, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "WW", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  

  scalar = visitmdscalarcreate(handle, "UA", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "VA", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)

  scalar = visitmdscalarcreate(handle, "T1", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "S1", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "RHO1", 4, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
 
  scalar = visitmdscalarcreate(handle, "VX", 2, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "VY", 2, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)
 
  
#  if defined (MULTIPROCESSOR)
  scalar = visitmdscalarcreate(handle, "EL_PID", 6, "TWOD_Mesh", 9,&
       & VISIT_VARCENTERING_ZONE)


  if(PAR)then
     ! FVCOM NODE INDEX VARIABLES
     scalar = visitmdscalarcreate(handle, "Node_Index_TWOD", 15, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Node_Index_Bathy", 16, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Node_Index_SSH", 14, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Node_Index_SigmaLevel", 21, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Node_Index_SigmaLayer", 21, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     
     ! FVCOM CELL INDEX VARIABLES
     scalar = visitmdscalarcreate(handle, "Cell_Index_TWOD", 15, "TWOD_Mesh", 9, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "Cell_Index_Bathy", 16, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "Cell_Index_SSH", 14, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "Cell_Index_SigmaLevel", 21, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "Cell_Index_SigmaLayer", 21, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_ZONE)
  end if
# endif
  

# if defined (SPHERICAL)

  scalar = visitmdscalarcreate(handle, "XYZ_U", 5, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_V", 5, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_WW", 6, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)


  IF (VISIT_ALL_VARS) THEN 
  scalar = visitmdscalarcreate(handle, "XYZ_UA", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_VA", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_WA", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  END IF

# endif

  
  IF (VISIT_ALL_VARS) THEN 

     ! grid metrics
     scalar = visitmdscalarcreate(handle, "ART", 3, "TWOD_Mesh", 9, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ART1", 4, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "ART2", 4, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)

     ! 2-d Arrays for the general vertical coordinates
     scalar = visitmdscalarcreate(handle, "Z", 1, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "ZZ", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "DZ", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     ! Can not plot DZZ, it has KBM2 node centered values. There is
     ! no mesh for that!
     !scalar = visitmdscalarcreate(handle, "DZZ", 3, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)

     ! 2-d flow variable arrays at elements
# if !defined (SEMI_IMPLICIT)
     scalar = visitmdscalarcreate(handle, "UARK", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "VARK", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
# endif

     scalar = visitmdscalarcreate(handle, "COR", 3, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !!$ DO NOT INCLUDE INTERPOLATED VALUES AT THE CELL CENTERS
     !  scalar = visitmdscalarcreate(handle, "H1", 2, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "D1", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "EL1", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "DT1", 3, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "ET1", 3, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "ELRK1", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "DTFA", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "CC_SPONGE", 9, "TWOD_Mesh", 8, VISIT_VARCENTERING_ZONE)
     ! 2-d flow variable arrays at nodes
     ! INCLUDED IN GENERAL VARIALBES
     scalar = visitmdscalarcreate(handle, "D", 1, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "DT", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "ET", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "EGF", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
# if !defined (SEMI_IMPLICIT)
     scalar = visitmdscalarcreate(handle, "ELRK", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
# endif
     ! SURFACE AND BOTTOM BOUNDARY CONDITIONS
     scalar = visitmdscalarcreate(handle, "CBC", 3, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "SWRAD", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
# if !defined (SEMI_IMPLICIT)
     scalar = visitmdscalarcreate(handle, "WUSURF2", 7, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WVSURF2", 7, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
# endif
     scalar = visitmdscalarcreate(handle, "WUBOT", 5, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WVBOT", 5, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WUSURF", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WVSURF", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     ! Skipping WTSURF: What is it?
     ! GROUND WATER DISCHARGE
     scalar = visitmdscalarcreate(handle, "BFWDIS", 6, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "BFWTMP", 6, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "BFWSLT", 6, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
# if !defined (SEMI_IMPLICIT)
     scalar = visitmdscalarcreate(handle, "BFWDIS2", 7, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
# endif
     ! RIVER DISCHARGE: add later!

     ! SURFACE MET FORCING
     scalar = visitmdscalarcreate(handle, "UUWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "VVWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
# if defined (SPHERICAL)
     scalar = visitmdscalarcreate(handle, "XYZ_UWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "XYZ_VWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "XYZ_WWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
# endif
# if !defined (SEMI_IMPLICIT)
     scalar = visitmdscalarcreate(handle, "QPREC2", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "QEVAP2", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
# endif
     scalar = visitmdscalarcreate(handle, "QEVAP", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "QPREC", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)

     ! 2 D Flow FLuxes
     scalar = visitmdscalarcreate(handle, "PSTX", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "PSTY", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADVUA", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADVVA", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADX2D", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADY2D", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "DRX2D", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "DRY2D", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADVX", 4, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADVY", 4, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)

     ! INTERNAL MODE ARRAYS- ELEMENTS
     scalar = visitmdscalarcreate(handle, "W", 1, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_ZONE)

     scalar = visitmdscalarcreate(handle, "UF", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "VF", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WT", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_ZONE)

     ! INTERNAL MODE ARRAYS- NODES
     scalar = visitmdscalarcreate(handle, "Q2", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "L", 1, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Q2L", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Q2", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)

# if defined (GOTM)
     scalar = visitmdscalarcreate(handle, "TKE", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "TEPS", 4, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
# endif

      !--------------------------------------------------------------
      !--------------------------------------------------------------
# if defined (ICE)

  scalar = visitmdscalarcreate(handle, "AICE", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "VICE", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "STRENGTH", 8, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)

  scalar = visitmdscalarcreate(handle, "UICE2", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "VICE2", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)

# if defined (SPHERICAL)
  scalar = visitmdscalarcreate(handle, "XYZ_UICE", 8, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_VICE", 8, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_WICE", 8, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
# endif

# endif
     !--------------------------------------------------------------
     !--------------------------------------------------------------



     scalar = visitmdscalarcreate(handle, "KM", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "KH", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "KQ", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "AAM", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)

     scalar = visitmdscalarcreate(handle, "TF1", 3, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "SF1", 3, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "WTS", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "WTTS", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)

     ! Baroclinic Pressure Gradients
     scalar = visitmdscalarcreate(handle, "DRHOX", 5, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "DRHOY", 5, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)

     ! Shape Coefficients: add later if needed

     ! Salinity and temperature bottom diffusion stuff

     ! Average arrays ?

     ! WHAT IS VISCOFH


     ! WHAT ABOUT BOUNDARY CONDITION VARIABLES?


  END IF

  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! VISIT MATERIALS
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  ! VISIT Materials are not working yet. Not sure Why? Causes Sig Sev...

  ! Add sigma layers material
  mat1 = visitmdmaterialcreate(handle, "Sigma_Layers", 12, "SigmaLayer_Mesh", 15)
  if(mat1.ne.VISIT_INVALID_HANDLE) then
     
     Do ind=1,KBM2
        
        write(ch3,'(i3)') ind
        
        if(ind<10 .and. ind>=1) then
           matname="Layer_00"//trim(adjustl(ch3))
        else if(ind<100 .and. ind>=10) then
           matname="Layer_0"//trim(adjustl(ch3))
        else if(ind<1000 .and. ind>=100) then
           matname="Layer_"//trim(adjustl(ch3))
        else
           if(MSR) WRITE(IPT,*)"In VisitGetMetaData: Layer Number > 1000"
           if(MSR) WRITE(IPT,*)"Causing error in visit material"
        end if
        
        err = visitmdmaterialadd(handle, mat1, matname, 9)
        
     End Do
  endif


  mat2 = visitmdmaterialcreate(handle, "Sigma_Levels", 12, "SigmaLevel_Mesh", 15)
  if(mat2.ne.VISIT_INVALID_HANDLE) then
     
     Do ind=1,KBM1
        
        write(ch3,'(i3)') ind

        if(ind<10 .and. ind>=1) then
           matname="Level_00"//trim(adjustl(ch3))
        else if(ind<100 .and. ind>=10) then
           matname="Level_0"//trim(adjustl(ch3))
        else if(ind<1000 .and. ind>=100) then
           matname="Level_"//trim(adjustl(ch3))
        else
           if(msr) WRITE(IPT,*)"In VisitGetMetaData: Level Number > 1000"
           if(msr) WRITE(IPT,*)"Causing error in visit material"
        end if

        err = visitmdmaterialadd(handle, mat2, matname, 9)


     End Do
  endif


  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! VISIT CONTROL COMMANDS
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  err = visitmdaddsimcommand(handle, "Halt/Step", 9, VISIT_CMDARG_NONE,1)

  err = visitmdaddsimcommand(handle, "Run", 3, VISIT_CMDARG_NONE,1)

  err = visitmdaddsimcommand(handle, "Dump Restart", 12, VISIT_CMDARG_NONE,1)

  err = visitmdaddsimcommand(handle, "Toggle Step", 11, VISIT_CMDARG_NONE,1)


!  err = visitmdaddsimcommand(handle, "i10", 3, VISIT_CMDARG_NONE,1)
!  err = visitmdaddsimcommand(handle, "i100", 4, VISIT_CMDARG_NONE,1)




  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! FVCOM EXPRESSIONS
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# if defined (SPHERICAL)

  err = visitmdexpressioncreate(handle, "Velocity",8,&
       & "{XYZ_U, XYZ_V, XYZ_WW}", 22, VISIT_VARTYPE_VECTOR) 

  err = visitmdexpressioncreate(handle, "Average_Velocity",16,&
       & "{XYZ_UA, XYZ_VA, XYZ_WA}", 24, VISIT_VARTYPE_VECTOR) 

  if(MSR) WRITE(IPT,*)"Expression create returned: ", err

# if defined (ICE)
  err = visitmdexpressioncreate(handle, "Ice_Velocity",12,&
       & "{XYZ_UICE, XYZ_VICE, XYZ_WICE}", 30, VISIT_VARTYPE_VECTOR) 
# endif

  err = visitmdexpressioncreate(handle, "Wind_Velocity",13,&
       & "{XYZ_UWIND, XYZ_VWIND, XYZ_WWIND}", 33, VISIT_VARTYPE_VECTOR) 

# else  
! if not spherical !!

  err = visitmdexpressioncreate(handle, "Velocity",8, "{U,V,WW}",&
       & 8, VISIT_VARTYPE_VECTOR) 

  err = visitmdexpressioncreate(handle, "Average_Velocity",16, "{UA,VA,0}",&
       & 9, VISIT_VARTYPE_VECTOR) 

# if defined (ICE)
  err = visitmdexpressioncreate(handle, "Ice_Velocity",12,&
       & "{UICE2,VICE2, 0}", 16, VISIT_VARTYPE_VECTOR) 
# endif

  err = visitmdexpressioncreate(handle, "Wind_Velocity",13,&
       & "{UUWIND,VVWIND,0}", 17, VISIT_VARTYPE_VECTOR) 

# endif


  visitgetmetadata = VISIT_OKAY


end function visitgetmetadata

!c---------------------------------------------------------------------------
!c visitgetmesh
!c---------------------------------------------------------------------------
integer function visitgetmesh(handle, domain, name, lname)
  USE MOD_VISIT
  USE ALL_VARS
  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8), intent(in) :: name
  integer, intent(in) :: handle, domain, lname
  integer :: error
  TYPE(VISITMESHTYPE), POINTER :: vmp ! VISIT MESH POINTER

  TYPE(VISITLAG), POINTER :: lmp ! VISIT MESH POINTER

  
  error = VISIT_ERROR

  if (domain .neqv. VisitParRank) then
     if(MSR) WRITE(IPT,*) "Something screwy is going on: visit domain= ",domain, "&
          &; VisitParRank= ",VisitParRank 
     visitgetmesh=VISIT_error
     return
  end if



  if (visitstrcmp(name, lname, "TWOD_Mesh", 9).eq.0) then
     
     vmp => VISIT_MESH(VISIT_TWODMESH)
     Call UPDATE_MESH(VMP,error)
     if (error .ne. 0) then ! CHECK FOR BAD RESULT
        visitgetmesh=VISIT_error
        return
     end if

    
  elseif (visitstrcmp(name, lname, "Bathymetry_Mesh",15).eq.0) then
     
     vmp => VISIT_MESH(VISIT_BATHYMESH)
     Call UPDATE_MESH(VMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=VISIT_error
        return        
     end if
     
  elseif (visitstrcmp(name, lname, "SSH_Mesh",8).eq.0) then
     
     vmp => VISIT_MESH(VISIT_SSHMESH)
     Call UPDATE_MESH(VMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=VISIT_error
        return        
     end if
     
     
  elseif (visitstrcmp(name, lname, "SigmaLayer_Mesh",15).eq.0) then


     vmp => VISIT_MESH(VISIT_LAYERMESH)
     Call UPDATE_MESH(VMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=VISIT_error
        return        
     end if
     
     
  elseif (visitstrcmp(name, lname, "SigmaLevel_Mesh",15).eq.0) then
     
     vmp => VISIT_MESH(VISIT_LEVELMESH)
     Call UPDATE_MESH(VMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=VISIT_error
        return        
     end if

  elseif (visitstrcmp(name, lname, "LAGRANGIAN_MESH", 15).eq.0) then
     
     lmp => VISIT_LAGDATA
     Call UPDATE_LAG(LMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=VISIT_error
        return        
     end if
     
     visitgetmesh = visitmeshpoint(handle, LMP%NDIMS, LMP%NODES, LMP%VX, LMP%VY, LMP%VZ)

     NULLIFY(LMP)
     return
     

  end if
  
   

  error = visitmeshunstructured3(handle,VMP%NDIMS, VMP%Nodes,&
       & VMP%Zones,VMP%GHOSTZONES, VMP%VX,&
       & VMP%VY, VMP%VZ, VMP%LCONN, VMP%NV&
       &,VMP%DATAOWNER)
  
  NULLIFY(VMP)

  visitgetmesh = error
end function visitgetmesh


!c---------------------------------------------------------------------------
!c visitgetscalar
!c---------------------------------------------------------------------------
integer function visitgetscalar(handle, domain, name, lname)
  USE ALL_VARS
  USE MOD_VISIT
  USE BCS
  USE MOD_PAR

      !--------------------------------------------------------------
      !--------------------------------------------------------------
# if defined (ICE)
      use ice_state
      use mod_ice
      use mod_ice2d
# endif
     !--------------------------------------------------------------
     !--------------------------------------------------------------


  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8),intent(in) :: name
  character(LEN=100) :: fnameVVV
  integer, intent(in) :: handle, domain, lname
  integer :: error, I
  integer, dimension(3) :: Vdims
  integer :: Vsize, ind


  integer, allocatable, dimension(:) :: V1di
  real*4, allocatable, dimension(:) :: V1df

  TYPE(VisitSphereVel), pointer :: VSV

  error = VISIT_ERROR


  ! BASIC Components of FVCOM 
  if (visitstrcmp(name, lname, "H",1).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, H(1:MT),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "EL",2).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, EL(1:MT),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "U",1).eq.0) then
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, U(1:NT,1:KBM1),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "V",1).eq.0) then
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, V(1:NT,1:KBM1),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WW",2).eq.0) then
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WW(1:NT,1:KBM1),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "T1",2).eq.0) then
     Vsize=KBM1*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, T1(1:MT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "S1",2).eq.0) then
     Vsize=KBM1*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, S1(1:MT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "RHO1",4).eq.0) then
     Vsize=KBM1*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, RHO1(1:MT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VX",2).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VX(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VY",2).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VY(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

# if defined(MULTIPROCESSOR)
  elseif (visitstrcmp(name, lname, "EL_PID",6).eq.0) then

     allocate(v1di(NT)); v1di=myid
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_TWOD",15).eq.0) then

     allocate(v1di(MT)); v1di=0
     v1di(1:M)=NGID(1:M)
     
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_Bathy",16).eq.0) then

     allocate(v1di(MT)); v1di=0
     v1di(1:M)=NGID(1:M)
     
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_SSH",14).eq.0) then

     allocate(v1di(MT)); v1di=0
     v1di(1:M)=NGID(1:M)
     
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_SigmaLevel",21).eq.0) then

     allocate(v1di(MT*KB)); v1di=0
     Do ind=1,KB
        v1di(((ind-1)*MT+1):((ind-1)*MT+M))=NGID(1:M)
     End Do
     
     Vsize=MT*KB
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_SigmaLayer",21).eq.0) then

     allocate(v1di(MT*KBM1)); v1di=0
     Do ind=1,KBM1
        v1di(((ind-1)*MT+1):((ind-1)*MT+M))=NGID(1:M)
     End Do
     
     Vsize=MT*KBM1
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)


     ! Cell Index arrays!
  elseif (visitstrcmp(name, lname, "Cell_Index_TWOD",15).eq.0) then

     allocate(v1di(NT)); v1di=0
     v1di(1:N)=EGID(1:N)
     
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Cell_Index_Bathy",16).eq.0) then

     allocate(v1di(NT)); v1di=0
     v1di(1:N)=EGID(1:N)
     
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Cell_Index_SSH",14).eq.0) then

     allocate(v1di(NT)); v1di=0
     v1di(1:N)=EGID(1:N)
     
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Cell_Index_SigmaLevel",21).eq.0) then

     allocate(v1di(NT*KBM1)); v1di=0
     Do ind=1,KBM1
        v1di(((ind-1)*NT+1):((ind-1)*NT+N))=EGID(1:N)
     End Do
     
     Vsize=NT*KB
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Cell_Index_SigmaLayer",21).eq.0) then

     allocate(v1di(NT*KBM1)); v1di=0
     Do ind=1,KBM2
        v1di(((ind-1)*NT+1):((ind-1)*NT+N))=EGID(1:N)
     End Do
     
     Vsize=NT*KBM2
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

# endif


# if defined (SPHERICAL)

! Velocity terms
  elseif (visitstrcmp(name, lname, "XYZ_U",5).eq.0) then
     VSV =>VisitSphericalVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereVel(VSV)
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelX,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_V",5).eq.0) then
     VSV =>VisitSphericalVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereVel(VSV)
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelY,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_WW",6).eq.0) then
     VSV =>VisitSphericalVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereVel(VSV)
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelZ,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)


 !! Average velocity terms
  elseif (visitstrcmp(name, lname, "XYZ_UA",6).eq.0) then
     VSV =>VisitSphericalAVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereAVel(VSV)
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelX,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_VA",6).eq.0) then
     VSV =>VisitSphericalAVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereAVel(VSV)
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelY,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_WA",6).eq.0) then
     VSV =>VisitSphericalAVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereAVel(VSV)
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelZ,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

! WIND Velocity terms
  elseif (visitstrcmp(name, lname, "XYZ_UWIND",9).eq.0) then
     VSV =>VisitSphericalWindVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereWindVel(VSV)
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelX,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_VWIND",9).eq.0) then
     VSV =>VisitSphericalWindVel
     if (VSV%Updated_Time .NE. Visit_Time_INT) &
          & Call UpdateSphereWindVel(VSV)
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelY,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_WWIND",9).eq.0) then
     VSV =>VisitSphericalWindVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereWindVel(VSV)
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelZ,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

# if defined (ICE)
! ICE Velocity terms
  elseif (visitstrcmp(name, lname, "XYZ_UICE",8).eq.0) then
     VSV =>VisitSphericalIceVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereICEVel(VSV)
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelX,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_VICE",9).eq.0) then
     VSV =>VisitSphericalIceVel
     if (VSV%Updated_Time .NE. Visit_Time_INT) &
          & Call UpdateSphereIceVel(VSV)
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelY,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_WICE",9).eq.0) then
     VSV =>VisitSphericalIceVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereIceVel(VSV)
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelZ,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

# endif 
! ICE


# endif 
! SPHERICAL



     ! LAGRANGIAN PARTICLE TRACKING STUFF GOES HERE!
  elseif (visitstrcmp(name, lname, "LAG_U",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VU, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "LAG_V",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VV, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "LAG_W",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VW, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "LAG_Z",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VZ, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

! based on number_of_scalars from particle.F

  elseif (visitstrcmp(name, lname, "LAG_S1",6).eq.0) then
     Vsize=Visit_LagData%Nodes
     allocate(v1df(vsize))
     DO I = 1,Vsize
        v1df(I) = Visit_LagData%VS(I)%S(1)
     END DO
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, V1DF, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)
     deallocate(v1df)

  elseif (visitstrcmp(name, lname, "LAG_S2",6).eq.0) then
     Vsize=Visit_LagData%Nodes
     allocate(v1df(vsize))
     DO I = 1,Vsize
        v1df(I) = Visit_LagData%VS(I)%S(2)
     END DO
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, V1DF, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)
     deallocate(v1df)

  elseif (visitstrcmp(name, lname, "LAG_S3",6).eq.0) then
     Vsize=Visit_LagData%Nodes
     allocate(v1df(vsize))
     DO I = 1,Vsize
        v1df(I) = Visit_LagData%VS(I)%S(3)
     END DO
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle,V1DF, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)
     deallocate(v1df)

  elseif (visitstrcmp(name, lname, "LAG_S4",6).eq.0) then
     Vsize=Visit_LagData%Nodes
     allocate(v1df(vsize))
     DO I = 1,Vsize
        v1df(I) = Visit_LagData%VS(I)%S(4)
     END DO
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, V1DF, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)
     deallocate(v1df)

  elseif (visitstrcmp(name, lname, "LAG_S5",6).eq.0) then
     Vsize=Visit_LagData%Nodes
     allocate(v1df(vsize))
     DO I = 1,Vsize
        v1df(I) = Visit_LagData%VS(I)%S(5)
     END DO
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, V1DF, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)
     deallocate(v1df)

  elseif (visitstrcmp(name, lname, "LAG_D",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VD, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


     !=================================================================
     !=================================================================
     !=================================================================
     ! END BASIC VARIABLES   THE REST ARE VARIABLES ACCESSED BY RUNNING 
     ! FVCOM WITH THE VISIT_OPT FLAG SET TO ADVANCED. THESE VARIABLES 
     ! ARE FOR ADVANCDED USERS WHO KNOW THEIR USAGE AND MEANING IN FVCOM


     ! GRID METRICS
  elseif (visitstrcmp(name, lname, "ART",3).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ART(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ART1",4).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ART1(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ART2",4).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ART2(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

     ! 2-d Arrays for the general vertical coordinates
  elseif (visitstrcmp(name, lname, "Z",1).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Z(1:MT,1:KB), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ZZ",2).eq.0) then
     Vsize=KBM1*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ZZ(1:MT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DZ",2).eq.0) then
     Vsize=KBM1*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DZ(1:MT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

     !  elseif (visitstrcmp(name, lname, "DZZ",3).eq.0) then
     !     Vsize=KBM1*MT
     !     Vdims=(/Vsize,1,1/)
     !     error  = visitscalarsetdata(handle, DZZ(1:MT,1:KBM1), Vdims , 1 &
     !          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

     ! 2-d flow variable arrays at elements

  elseif (visitstrcmp(name, lname, "UA",2).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UA(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VA",2).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VA(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


  elseif (visitstrcmp(name, lname, "UARK",4).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UARK(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VARK",4).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UARK(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


  elseif (visitstrcmp(name, lname, "COR",3).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, COR(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DTFA",4).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DTFA(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "CC_SPONGE",9).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, CC_SPONGE(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


     ! 2-d flow variable arrays at nodes
  elseif (visitstrcmp(name, lname, "D",1).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, D(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DT",2).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DT(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ET",2).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ET(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "EGF",3).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, EGF(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ELRK",4).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ELRK(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

     ! SURFACE AND BOTTOM BOUNDARY CONDITIONS

  elseif (visitstrcmp(name, lname, "CBC",3).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, CBC(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "SWRAD",5).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, SWRAD(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WUSURF2",7).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WUSURF2(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WVSURF2",7).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WVSURF2(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


  elseif (visitstrcmp(name, lname, "WUBOT",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WUBOT(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WVBOT",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WVBOT(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WUSURF",6).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WUSURF(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WVSURF",6).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WVSURF(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


     ! Ground water variables!!!
     ! BFWDIS variables are allocated BFWDIS(IBFW) 
         ! => to nodes mesh!
  elseif (visitstrcmp(name, lname, "BFWDIS",6).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, BFWDIS(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "BFWTMP",6).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, BFWTMP(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "BFWSLT",6).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, BFWSLT(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


  elseif (visitstrcmp(name, lname, "BFWDIS2",7).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, BFWDIS2, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

     ! END of Ground water variables!

     ! SURFACE MET FORCING
  elseif (visitstrcmp(name, lname, "UUWIND",6).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UUWIND(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VVWIND",6).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VVWIND(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "QPREC2",6).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, QPREC2(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "QPREC",6).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, QPREC(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "QEVAP2",6).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, QEVAP2(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "QEVAP",6).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, QEVAP(1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "PSTX",4).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, PSTX(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "PSTY",4).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, PSTY(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADVUA",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADVUA(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADVVA",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADVVA(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADX2D",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADX2D(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADY2D",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADY2D(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DRX2D",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DRX2D(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DRY2D",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DRY2D(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADVX",4).eq.0) then
     Vsize=NT*KBM1
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADVX(1:NT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADVY",4).eq.0) then
     Vsize=NT*KBM1
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADVY(1:NT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


     ! INTERNAL MODE ARRAYS- ELEMENTS -- U, V, WW are in the basic
  elseif (visitstrcmp(name, lname, "W",1).eq.0) then
     Vsize=KBM2*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, W(1:NT,2:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "UF",2).eq.0) then
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VF(1:NT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VF",2).eq.0) then
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VF(1:NT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WT",2).eq.0) then
     Vsize=KBM2*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WT(1:NT,2:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


     ! INTERNAL MODE ARRAYS- NODES
  elseif (visitstrcmp(name, lname, "Q2",2).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Q2(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "L",1).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, L(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "Q2L",3).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Q2L(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

# if defined (GOTM)
  elseif (visitstrcmp(name, lname, "TKE",3).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, TKE(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "TEPS",4).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, TEPS(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

# endif

  elseif (visitstrcmp(name, lname, "KM",2).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, KM(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "KH",2).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, KH(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "KQ",2).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, KQ(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "AAM",3).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, AAM(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "TF1",3).eq.0) then
     Vsize=KBM1*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, TF1(1:MT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "SF1",3).eq.0) then
     Vsize=KBM1*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, SF1(1:MT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WTS",3).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WTS(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WTTS",4).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WTTS(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WTS",3).eq.0) then
     Vsize=KB*MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WTS(1:MT,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DRHOX",5).eq.0) then
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DRHOX(1:NT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DRHOY",5).eq.0) then
     Vsize=KBM1*NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DRHOY(1:NT,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


      !--------------------------------------------------------------
      !--------------------------------------------------------------
# if defined (ICE)
  elseif (visitstrcmp(name, lname, "AICE",4).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, AICE(1,1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VICE",4).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, AICE(1,1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "STRENGTH",8).eq.0) then
     Vsize=MT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, STRENGTH(1,1:MT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "UICE2",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UICE2(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VICE2",5).eq.0) then
     Vsize=NT
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VICE2(1:NT), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

# endif
     !--------------------------------------------------------------
     !--------------------------------------------------------------

  endif


  visitgetscalar = error
end function visitgetscalar


!c---------------------------------------------------------------------------
!c visitgetcurve
!c---------------------------------------------------------------------------
integer function visitgetcurve(handle, name, lname)
  implicit none
  include "visitfortransiminterface.inc"
  character*8 name
  integer     handle, lname

  visitgetcurve = VISIT_ERROR
end function visitgetcurve

!c---------------------------------------------------------------------------
!c visitgetdomainlist
!c---------------------------------------------------------------------------
integer function visitgetdomainlist(handle)
  USE MOD_VISIT, only: VisitParRank
  USE LIMS, only : NPROCS
  implicit none
  include "visitfortransiminterface.inc"
  integer handle

  visitgetdomainlist = visitsetdomainlist(handle, NPROCS, VisitParRank,1)

  visitgetdomainlist = VISIT_OKAY
end function visitgetdomainlist



!c---------------------------------------------------------------------------
!c visitgetmaterial
!c---------------------------------------------------------------------------
integer function visitgetmaterial(handle, domain, name, lname)
  USE LIMS
  USE MOD_VISIT
  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8), intent(in) :: name
  integer, intent(in) :: handle, domain, lname
  integer:: rslt, ind, ind2
  character(LEN=3) :: ch3
  character(LEN=20) :: matname
  INTEGER, allocatable :: MatNo(:)
  integer, dimension(3):: dims

  if (visitstrcmp(name, lname, "Sigma_Layers",12).eq.0) then
     
     allocate(MatNo(KBM2))

     dims(1)=NT*KBM2
     dims(2)=1
     dims(3)=1

     rslt = visitMaterialSetDims(handle, dims, 1)
     Do ind =1, KBM2
        write(ch3,'(i3)') ind
        
        if(ind<10 .and. ind>=1) then
           matname="Layer_00"//trim(adjustl(ch3))
        else if(ind<100 .and. ind>=10) then
           matname="Layer_0"//trim(adjustl(ch3))
        else if(ind<1000 .and. ind>=100) then
           matname="Layer_"//trim(adjustl(ch3))
        else
           if(MSR) WRITE(IPT,*)"In GetMaterial: Layer Number > 1000"
           if(MSR) WRITE(IPT,*)"Causing error in visit material"
        end if

        MatNo(ind) = visitMaterialAdd(handle, matname, 9)
        
     End Do
     
     Do ind=1, KBM2
        Do ind2=1,NT
           rslt = visitMaterialAddClean(handle,(ind-1)*NT+ind2,MatNo(ind))
        End Do
     End Do

     Deallocate(MatNo)
     rslt = VISIT_OKAY


  elseif (visitstrcmp(name, lname, "Sigma_Levels",12).eq.0) then
     
     allocate(MatNo(KBM1))

     dims(1)=NT*KBM1
     dims(2)=1
     dims(3)=1
     
     rslt = visitMaterialSetDims(handle, dims, 1)
     
     Do ind =1, KBM1

        write(ch3,'(i3)') ind
        
        if(ind<10 .and. ind>=1) then
           matname="Level_00"//trim(adjustl(ch3))
        else if(ind<100 .and. ind>=10) then
           matname="Level_0"//trim(adjustl(ch3))
        else if(ind<1000 .and. ind>=100) then
           matname="Level_"//trim(adjustl(ch3))
        else
           if(MSR) WRITE(IPT,*)"In GetMaterial: Level Number > 1000"
           if(MSR) WRITE(IPT,*)"Causing error in visit material"
        end if

        MatNo(ind) = visitMaterialAdd(handle, matname, 9)
     End Do
     

     Do ind=1, KBM1
        Do ind2=1,NT
           rslt = visitMaterialAddClean(handle,(ind-1)*NT+ind2,MatNo(ind))
        End Do
     End Do

     Deallocate(MatNo)
     rslt = VISIT_OKAY     


  else
     rslt = VISIT_ERROR
  end if
  


  visitgetmaterial = rslt
end function visitgetmaterial


#else
! if visit is not defined compile a dummy subroutine!

subroutine init_visit
  USE MOD_UTILS
  implicit none

  IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! VISIT SIMULATION MODE IS OFF"

end subroutine init_visit

#endif



























