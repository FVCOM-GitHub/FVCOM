Program fvcom2gnome 

  Use mod_fvcom2gnome
  Use mod_utils, only : initialize_control

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string

  !---------------------------------------------------------------------
  ! initialize basic FVCOM variables (like MSR) 
  !---------------------------------------------------------------------
  call initialize_control("fvcom2gnome")

  !---------------------------------------------------------------------
  ! setup the MPI environment and shutdown if parallel
  !---------------------------------------------------------------------
# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
  IF (PAR) CALL FATAL_ERROR &
       & ("THIS PROGRAM WILL NOT RUN IN PARALLEL", & 
       & "STOP AND RERUN WITH JUST ONE PROCESSOR!")
# endif

  !---------------------------------------------------------------------
  !parse arguments from the command line
  !---------------------------------------------------------------------
  Call Parse_Commandline(CVS_id,CVS_date,CVS_name,CVS_revision)

  !---------------------------------------------------------------------
  !read mesh and number of frames from inputfile 
  !---------------------------------------------------------------------
  Call Setup_Input

  !---------------------------------------------------------------------
  !generate boundary
  !---------------------------------------------------------------------
  Call Get_Edgelist 
  Call Order_Edges
  Call Set_Seg_Neighbors
  Call Get_Seglist
  Call Gen_GNOME_Bndry

  !---------------------------------------------------------------------
  !create output file 
  !---------------------------------------------------------------------
  !Call Setup_Output_NCTOOLS (issues with GNOME, shutdown for now)
  Call Setup_Output

  !---------------------------------------------------------------------
  !loop over frames and dump velocities at vertices to gnome file
  !---------------------------------------------------------------------
  Call Extract

  !---------------------------------------------------------------------
  !Nullify Pointers and Deallocate Memory
  !---------------------------------------------------------------------
  Call Cleanup 

End Program fvcom2gnome
