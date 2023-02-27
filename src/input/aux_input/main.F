PROGRAM AUXILIARY_INPUT
  use mod_input
  use mod_auxiliary
  use mod_setup
  USE MOD_SET_TIME
  USE MOD_REPORT
  IMPLICIT NONE
  INTEGER :: STATUS, I
  CHARACTER(LEN=4) :: BFLAG
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string



  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("Other Input")

# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
# endif

  IF (PAR) CALL FATAL_ERROR &
       & ("THIS PROGRAM IS WAY TO SHORT TO BOTHER WITH WRITING AS A PARALLE CODE.",&
       & "STOP AND RERUN WITH JUST ONE PROCESSOR!")


  CALL NAME_LIST_INITIALIZE ! FROM MAIN CODE

  CALL INITIALIZE_NAMELIST ! FOR THIS PROGRAM

  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ========================="
     WRITE(IPT,*) "!   BEGIN CREATE AUX DATA"
     WRITE(IPT,*) "! ========================="
  END IF


  CALL READ_NAMELIST

  write(ipt,*) "HERE,",timezone 
  CALL CNTRL_PRMTRS
  write(ipt,*) "HERE" 

  CALL CHECK_IO_DIRS

  CALL OPEN_FILES

  CALL READ_COLDSTART_SIGMA
  CLOSE(SIGMAUNIT)

  CALL READ_COLDSTART_GRID(GRIDUNIT,MGL,NGL,NV)
  m = MGL
  mt = MGL
  n = ngl
  nt = ngl


  CALL ALLOCATE_SPACE

  ALLOCATE(X_LCL(0:MGL),Y_LCL(0:MGL))
  CALL READ_COLDSTART_COORDS(GRIDUNIT,MGL,X_LCL,Y_LCL)
  CLOSE(GRIDUNIT)

  CALL COORDINATE_UNITS(X_LCL,Y_LCL)
  CALL SETUP_CENTER_COORDS

  CALL READ_COLDSTART_DEPTH(DEPTHUNIT,MGL,X_LCL,Y_LCL,H)
  CALL N2E2D(H,H1)
  CLOSE(DEPTHUNIT)

  KBM1 = KB - 1
  KBM2 = KB - 2
  CALL Setup_Sigma
  CALL SETUP_SIGMA_DERIVATIVES

  DEALLOCATE(X_LCL,Y_LCL)
!====================================================================
! CREATE YOUR SPECIALIZED INPUT DATA HERE AND SAVE TO A FILE!

  print*,'STOP HERE 1'
  CALL CREATE_HVC(NN_HVC,CC_HVC)
  print*,'STOP HERE 2'

  CALL CREATE_BOTTOM_ROUGHNESS(CC_Z0B)
  print*,'STOP HERE 3'

  CALL DUMP_AUX
  print*,'STOP HERE 4'

!====================================================================
  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN


END PROGRAM AUXILIARY_INPUT
