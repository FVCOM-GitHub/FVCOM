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

module mod_lag
# if defined(LAG_PARTICLE)


  ! MODIFIED BY DAVID FOR THE NEW INPUT: NEEDS TESTING AND EDIT THE HEADER!

  ! ALL TIMES ARE NOW EITHER TYPE(TIME) OR IN DAYS!

  !------------------------------------------------------------------------------
  ! Lagrangian particle tracking using linked lists
  !
  !  Current Capabilities (September 08, 2006)
  !     1.) Multiprocessor Formulation:
  !           Subdomains will only traject particles contained in their interior
  !           Particles passing across interprocessor boundaries are exchanged 
  !           between neighboring processors.
  !     2.) Restart Capability:
  !           The code can be restarted and the particle output file will
  !           act as if the code had run with intermission.  Read
  !           below for details.
  !     3.) General tracking:  particles are trajected at each internal time step
  !         (dti) using a 4 stage RK scheme and the (u,v,sigma) velocity field.
  !
  !     4.) Output File Contents
  !         particle position, (x,y,z, + sigma) as well as an associated scalar
  !         variable at the particle position.  Current options include salinity
  !         temperature, densities, and eddy viscosities but can be modified 
  !         easily to provide user defined fields such as a biological quantitiy
  !         or tracer concentration. 
  !
  !     5.) initial particle position file.  This file is in netcdf format.
  !         It must contain information about number of particles, location
  !         (x,y,z) as well as global parameters including an information
  !         string, time of last output, number of outputs made.  The last
  !         two parameters are used when restarting the model.  To create
  !         an initial position file, start from the create_initial_lagpos.f90
  !         utility supplied in the utilities folder of your FVCOM distribution.
  !         You can use this script to easily convert an ascii particle file
  !         to the FVCOM netcdf format or you can generate a grid of particles
  !         within a user defined box.
  !
  !  Future Capabilities
  !     1.) Horizontal Random Walk.  
  !         To mimic the affects of turbulent eddies, we have use these algorithms. 
  !         They are quite varied, a bit hacky, and have the dangerous property of
  !         perhaps violating local CFL conditions which will mean that particles
  !         passing through domain boundaries will be forever lost.  A particle
  !         must first be detected inside a domain's halo element to know that it should
  !         be sent to the neighboring processor who contains that particular element
  !         in their interior.  If it jumps right into the other processors interior.
  !         It will be lost.  There are ways of dealing with this, of course, but I don't 
  !         suspect it will happen anytime soon.
  !
  !     2.) Vertical random walk. 
  !         This will be easier to implement, but the wide variety of models (read: 
  !         a lot of parameters) will mean a separate inputfile may be necessary
  !         for their specification.  Probably coming soon.
  !    
  !  =======Notes====
  !  Note:  there is no specific array holding the particles because
  !         they may traverse from one processor subdomain to another.
  !         The particles are stored in the linked list 'particle_list'.
  !         Particles may be added or deleted from the list (see linklist.F)
  !         To access a particle property, for example, id, loop through the
  !         linked list as follows:
  !
  ! 1.)  lp  => particle_list%first%next
  ! 2.)  do
  ! 3.)    if(.not. associated(lp) ) exit      
  ! 4.)      write(*,*)'this is the current particles id: ',lp%v%id   
  ! 5.)     lp => lp%next                     
  ! 6.)  end do
  !
  !             A list is comprised of 0 or more nodes
  ! Line 1 sets the pointer to the first node in the list. (lp)
  ! Line 2 starts an infinite loop
  ! Line 3 exits when the end of the list is reached
  ! Line 4 shows example usage:  note lp is the node and lp%v is the particle
  ! Line 5 move to the next node, without this line you have an infinite loop 
  ! Line 6 is the end of the loop
  !
  ! Reference: Object Oriented Programming via Fortran 90/95  , Ed Akin
  !
  ! Note, you must include NETCDF_IO in the Makefile to build the Lagrangian library
  ! This will require you to build the netcdf librarie for your system.
  ! This is a result of a code modernization effort and is part of conspiracy 
  !   against 'old school' programmers.  Well, maybe it is a conspiracy, but
  !   tough luck.
  !
  ! Runtime Parameters (to be placed in file:  casename_run.dat):
  ! Must be capital letters.
  !
  !  LAG_INTERVAL:  output interval in model hours
  !  LAG_INPFILE:   input file name (will read INPDIR directory) 
  !  LAG_RESFILE:   restart file (will dump to executable directory)
  !  LAG_OUTFILE:   output file  (will write to OUTDIR/netcdf directory)
  !  LAG_SCAL:      scalar variable to include with output
  !  LAG_COLD_START logical (T/F) giving startup type
  !
  !  LAG_COLD_START = T
  !       --> a new outputfile will be created
  !       --> initial positions will be output to outputfile 
  !       --> dumps are made every subsequent LAG_INTERVAL to outputfile
  !
  !  LAG_COLD_START = F (assumes input file generated by restart dump)
  !       --> new outputfile will not be created
  !       --> time of last dump is obtained from input file
  !       --> stack # of last dump obtained from input file 
  !       --> next dump made to stack# + 1 at time of last dump + LAG_INTERVAL 
  !       --> initial positions not dumped to file (the assumed initial position
  !           would be the first output in the file.
  !
  !  LAG_SCAL
  !       current choices
  !       s1   (salinity at nodes)
  !       t1   (temperature at nodes)
  !       rho1 (density at nodes)
  !       km   (eddy viscosity)
  !       kh   (vertical eddy diffusivity)
  !       user_defined
  !          --> users must set variable names in get_ud_scalnames
  !          --> users must supply 3D field in add_scalar 
  !          --> recompile after modification!
  !       search the code for "users" to see where to modify
  !
  !  General information on restarting the code.
  !    If you restart the code but don't want to interrupt your particle output 
  !    file:
  !
  !    1.) run the code.  It will dump a particle restart file (designated by
  !        LAG_RESFILE in the casename_run.dat file), every IRESTART iterations.
  !        It will also dump one at the end of the run.
  !    2.) move the LAG_RESFILE to your INPDIR directory.  Rename if you want.
  !        Let's say you name it brainless_coding.nc
  !    3.) modify the following parameters in your runtime file 
  !           RESTART = hot_start  (may already be set to this)
  !           LAG_COLD_START = F
  !           LAG_INPFILE = brainless_coding.nc 
  !    4.) restart the code.
  !
  !    The restart file you used as your LAG_INPFILE contains information about
  !    when the last output was made and the stack# of the output so that the  
  !    Lagrangian tracking module knows when to make the last dump and what
  !    stack# to use in the netcdf file.  By setting LAG_COLD_START=F, the code
  !    knows not to create a new outputfile, but rather use the othe one.
  !
  !    If you do something like restart the code but restart with a different
  !    number of particles but use LAG_COLD_START = F, I cannot be responsible
  !    for what might happen.  In simulations performed during the  extensive 
  !    testing procedure used for every line of FVCOM code, we found that results
  !    can range from crashed code to broken keyboards.
  !
  !------------------------------------------------------------------------------
  use mod_prec
  use linked_list
  use particle_class
  USE MOD_NCTOOLS
  use mod_time
  use mod_utils
  USE MOD_PAR
  USE CONTROL
  implicit none

  private

  public :: set_lag
  public :: lag_update
  public :: output_lag_restart
  public :: add_scalars
  
  !===========================================================
  ! ADDED FOR VISIT
  public :: particle_list
  !===========================================================


  !global variables
  type(link_list) :: particle_list    !linked list holding particles
  type(particle), allocatable, target :: Global_particles(:) ! STORAGE FOR I/O

  ! Some data pointers for IO
  real(sp), pointer,dimension(:) :: Xpnt,Ypnt,Zpnt, TBpnt, TEpnt,&
       & PLpnt, Spnt,latpnt, lonpnt
  real(sp), pointer,dimension(:) :: TMPpnt,SALpnt,KMpnt,KHpnt,DENpnt&
       &,Upnt,Vpnt,Wpnt,Opnt
  INTEGER, pointer,dimension(:)  :: Gpnt

  ! logicals controlling vertical position specification at startup
  logical :: init_with_zloc  = .false.
  logical :: init_with_sigma = .false.

  ! Some time pointers for IO
  INTEGER, pointer,dimension(:) ::TBEGI1,TENDI2
  real(sp), pointer,dimension(:) ::tbegf,tendf
  character(len=80), pointer, dimension(:) :: tbegc,tendc
  TYPE(TIME),POINTER,DIMENSION(:) :: TBEG,TEND


  integer :: nlag_gl = 0              !# global particles
  integer :: nlag = 0                 !# particles in local domain
  integer, allocatable :: ney_list(:) !list of domain neighbor processor ids

  !output control
  TYPE(TIME) :: lag_interval       ! LAG OUTPUT INTERVAL (MuSeconds)
  TYPE(TIME) :: lag_next           ! TIME OF LAST OUTPUT
  INTEGER :: LAG_STACK

  CHARACTER(LEN=80),ALLOCATABLE :: LAG_SCALARS(:)

  !netcdf variables - output file
  INTEGER, PARAMETER :: LAG_OUT_PROC=1 ! SET OUTPUT TO MASTER PROC
  TYPE(NCFILE), POINTER :: NC_LAG
  TYPE(NCFILE), POINTER :: NC_LAGSTART
  TYPE(NCFILE), POINTER :: NC_LAGRESTART


contains !------------------------------------------------------------------!
  ! set_lag             :   read in initial particle locations       !
  !                     :   determine initial host elements          !
  ! lag_update          :   update particle location and write data  !
  ! traject             :   use 4 stage rk scheme to track particle  !
  ! exchange_particles  :   pass particles among processors          !
  ! output_lag          :   gather particles and call to output      !
  ! dump_particles_ncd  :   dump particle information to netcdf file !
  ! collect_particles    :   gather particles to master processor     !
  ! and many more all for a guaranteed low price!!                   !        
  ! -----------------------------------------------------------------!


  !==============================================================================|
  subroutine set_lag 
    !------------------------------------------------------------------------------|
    !  read control parameters from input file
    !  read initial positions from LAG_INPFILE
    !  set remaining particle parameters (sigma location, pathlength, etc) 
    !  report
    !------------------------------------------------------------------------------|
    use mod_set_time
    use all_vars
    use mod_nctools
    implicit none
    integer              :: i,ii,icnt,iscan,ierr
    logical              :: fexist,check,FOUND
    character(len=120)   :: fname,temp
    type(particle),POINTER :: p
    type(particle)       :: psize
    integer              :: status
    integer, allocatable :: tmp(:)
    type(link_node), pointer      :: lp,nlp
    real(sp)             :: ftime

    logical lag_coldstart

    TYPE(TIME) :: tminus,lag_time

    TYPE(NCFILE), pointer :: NCF
    TYPE(NCATT), pointer :: ATT
    TYPE(NCDIM), pointer :: DIM
    TYPE(NCVAR), pointer :: VAR
    character(len=80) :: tzone

    ! IDEAL TIME FLAG
    CHARACTER(LEN=4) :: FLAG
    INTEGER(ITIME) :: IINT_CNT


    !----------------------------------------------------------
    ! initialize particle list (linked list structure)        
    !----------------------------------------------------------
    particle_list = new_list()


    !-------------------------------------------------------------
    !  Read particle tracking parameters from runtime file       
    !-------------------------------------------------------------
    if(.NOT. LAG_PARTICLES_ON)then

       IF(DBG_SET(DBG_LOG)) write(ipt,*)'!  # Lagrangian tracking :  inactive'
       return
    end if
    IF(DBG_SET(DBG_LOG)) write(ipt,*)'!  # Lagrangian tracking :  active'


    !-------------------------------------------------------------
    !  Initialize MPI_PARTICLE datatype       
    !-------------------------------------------------------------
    IF(PAR) CALL CREATE_MPI_PARTICLE


    !-------------------------------------------------------------
    ! Determine which scalars are in the output
    !-------------------------------------------------------------

    CALL SPLIT_STRING(LAG_SCAL_CHOICE,",",LAG_SCALARS)
    IF (size(LAG_SCALARS) .GT. size(psize%s) )CALL FATAL_ERROR&
         &("SET_LAG: You must edit particle.F and mod_lag.F to allow",&
         & "more scalar output terms!",&
         & "You asked for the following output:"//TRIM(LAG_SCAL_CHOICE))

    !-------------------------------------------------------------
    !  Determine run file based time intervals
    !-------------------------------------------------------------

    ! TIME OF FIRST OUTPUT
    IF(USE_REAL_WORLD_TIME) THEN
       lag_next = READ_DATETIME(LAG_FIRST_OUT,DATE_FORMAT,TIMEZONE,status)
       if (.not. status) &
            & Call Fatal_Error("Could not read the date string LAG_FIRST_OUT: "//trim(LAG_FIRST_OUT))          
    ELSE
       CALL IDEAL_TIME_STRING2TIME(LAG_FIRST_OUT,FLAG,lag_next,IINT_CNT)

       IF (FLAG == 'step') lag_next = (IINT_CNT - IINT) * IMDTI + RUNFILE_STARTTIME

    END IF

    ! THE OUTPUT INTERVAL
    CALL IDEAL_TIME_STRING2TIME(LAG_OUT_INTERVAL,FLAG,lag_interval,IINT_CNT)

    IF (FLAG == 'step') THEN
       lag_interval = IINT_CNT * IMDTI
    ELSE
       IF (mod(lag_interval,IMDTI)/= ZeroTime) Call Fatal_error&
            ("The Lagrangian Particle output interval must be an integer number of time steps!",&
            &"But it may be written as seconds,cycles or days.")
    END IF

    ! SANITY CHECK
    IF(lag_interval .LE. ZEROTIME) CALL FATAL_ERROR&
         & ('SET_LAG: PARTICLE OUTPUT INTERVAL IS LE TO ZERO (check your run file!)')

    !-------------------------------------------------------------
    !  Setup particle IO files and time intervals
    !-------------------------------------------------------------

    LAG_STACK=0
    ! FIND STARTUP FILE
    SELECT CASE(STARTUP_TYPE)
    CASE(STARTUP_TYPE_COLDSTART)

       IF(len_trim(LAG_START_FILE)==0) call fatal_error&
            &("LAG_START_FILE is not set in the name list!")

       fname = trim(INPUT_DIR)//trim(LAG_START_FILE)
       ! No need to adjust time parameters

    CASE(STARTUP_TYPE_HOTSTART) 
       IF(len_trim(LAG_START_FILE)==0) call fatal_error&
            &("LAG_START_FILE is not set in the name list!")

       fname = trim(INPUT_DIR)//trim(LAG_START_FILE)
       ! No need to adjust time parameters

    CASE(STARTUP_TYPE_CRASHRESTART)
       ! LOOK FOR THE RESTART FILE FIRST
       IF(len_trim(LAG_RESTART_FILE)==0) call fatal_error&
            &("LAG_RESTART_FILE is not set in the name list!")

       fname = trim(OUTPUT_DIR)//trim(LAG_RESTART_FILE)

       
       DO WHILE ( LAG_NEXT .LT. StartTime )
          LAG_STACK = LAG_STACK+1
          LAG_NEXT = LAG_NEXT + lag_interval
       END DO


    END SELECT


    ! SANITY CHECK FOR TIMES
    IF (LAG_NEXT .LT. StartTime) CALL FATAL_ERROR&
         &('SET_LAG: You specified a Lagragian Particle output before your model starts!')

    IF (LAG_NEXT .GT. EndTime) CALL FATAL_ERROR&
         &('SET_LAG: You specified a Lagragian Particle output before your model starts!')


    IF ( MOD((Lag_next - startTime),IMDTI)/= zerotime) CALL FATAL_ERROR &
         &("The Lag_First_Out value in the run file invalid",&
         & "Output must occur at an integer number of timesteps from the start time.")


    ! OPEN IT: will return an error if file does not exist
    CALL NC_INIT(NCF, fname)
    If(.not. NCF%OPEN) then
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       NC_LAGSTART => NCF
       FILEHEAD => ADD(FILEHEAD,NCF)
    end if

    !-------------------------------------------------------------
    !  Read Initial data
    !-------------------------------------------------------------

    ATT => FIND_ATT(NC_LAGSTART,'info_string',FOUND)

    DIM => FIND_DIM(NC_LAGSTART,'nparticles',FOUND)
    IF(.not. FOUND) DIM => FIND_DIM(NC_LAGSTART,'nlag',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         ("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE DIMENSION: nparticles OR nlag")
    nlag_gl = DIM%DIM


    ! ALLOCATE GLOBAL PARTICLE ARRAY
    CALL NEW_PARTICLES(Global_Particles, nlag_gl)
    ! SET SOME POINTERS TO THIS NEW DATA...    
    Spnt => GLOBAL_PARTICLES%x(3)
    Zpnt => GLOBAL_PARTICLES%zloc

    ! READ IN THE DATA
    IF(has_unlimited(NC_LAGSTART)) THEN

       lag_coldstart=.false.

       DIM => FIND_UNLIMITED(NC_LAGSTART,FOUND)

       lag_stack = 1
       DO WHILE (lag_time /= IntTime)
          lag_stack = lag_stack+1
          IF(lag_stack .GT. dim%dim) call fatal_error&
               ("SET_LAG: COULD NOT FIND VALID START TIME IN LAG RESTART FILE")

          lag_time = GET_FILE_TIME(NC_LAGSTART,lag_stack)

       END DO

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
       ! PARTICLE X LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"lon",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: lon")
       LONpnt => GLOBAL_PARTICLES%x(1)
       CALL NC_CONNECT_PVAR(VAR,LONpnt)
       CALL NC_READ_VAR(VAR,lag_stack)

       ! PARTICLE Y LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"lat",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: lat")
       LATpnt => GLOBAL_PARTICLES%x(2)
       CALL NC_CONNECT_PVAR(VAR,LATpnt)
       CALL NC_READ_VAR(VAR,lag_stack)

       ALLOCATE(Xpnt(nlag_gl))
       ALLOCATE(Ypnt(nlag_gl))
       Xpnt = 0.0_SP
       Ypnt = 0.0_SP

# else
       ! PARTICLE X LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"x",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: x")
       Xpnt => GLOBAL_PARTICLES%x(1)
       CALL NC_CONNECT_PVAR(VAR,Xpnt)
       CALL NC_READ_VAR(VAR,lag_stack)

       ! PARTICLE Y LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"y",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: y")
       Ypnt => GLOBAL_PARTICLES%x(2)
       CALL NC_CONNECT_PVAR(VAR,Ypnt)
       CALL NC_READ_VAR(VAR,lag_stack)

       ALLOCATE(LONpnt(nlag_gl))
       ALLOCATE(LATpnt(nlag_gl))
       LONpnt = 0.0_SP
       LATpnt = 0.0_SP

!---------------------------------------------------------------
# endif
!---------------------------------------------------------------


       ! PARTICLE Z(meters) LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"z",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: z")
       CALL NC_CONNECT_PVAR(VAR,Zpnt)
       CALL NC_READ_VAR(VAR,lag_stack)
       init_with_zloc  = .true.


       ! PARTICLE PATHLENGTH
       VAR => FIND_VAR(NC_LAGSTART,"pathlength",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: pathlength")
       PLpnt => GLOBAL_PARTICLES%pathlength
       CALL NC_CONNECT_PVAR(VAR,PLpnt)
       CALL NC_READ_VAR(VAR,lag_stack)

       ! PARTICLE MARK (0 IN DOMAIN, 1 OUT OF DOMAIN)
!!$       VAR => FIND_VAR(NC_LAGSTART,"mark",FOUND)
!!$       IF(.not. FOUND) CALL FATAL_ERROR&
!!$            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: mark")
!!$       Mpnt => GLOBAL_PARTICLES%mark
!!$       CALL NC_CONNECT_PVAR(VAR,Mpnt)
!!$       CALL NC_READ_VAR(VAR,lag_stack)

    ELSE ! NO TIME DIMENSION IN FILE

       lag_coldstart=.TRUE.

       
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
       ! PARTICLE X LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"lon",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: lon")
       LONpnt => GLOBAL_PARTICLES%x(1)
       CALL NC_CONNECT_PVAR(VAR,LONpnt)
       CALL NC_READ_VAR(VAR)

       ! PARTICLE Y LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"lat",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: lat")
       LATpnt => GLOBAL_PARTICLES%x(2)
       CALL NC_CONNECT_PVAR(VAR,LATpnt)
       CALL NC_READ_VAR(VAR)

       ALLOCATE(Xpnt(nlag_gl))
       ALLOCATE(Ypnt(nlag_gl))
       Xpnt = 0.0_SP
       Ypnt = 0.0_SP


# else

       ! PARTICLE X LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"x",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: x")
       Xpnt => GLOBAL_PARTICLES%x(1)
       CALL NC_CONNECT_PVAR(VAR,Xpnt)
       CALL NC_READ_VAR(VAR)

       ! PARTICLE Y LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"y",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: y")
       Ypnt => GLOBAL_PARTICLES%x(2)
       CALL NC_CONNECT_PVAR(VAR,Ypnt)
       CALL NC_READ_VAR(VAR)


       ALLOCATE(LONpnt(nlag_gl))
       ALLOCATE(LATpnt(nlag_gl))
       LONpnt = 0.0_SP
       LATpnt = 0.0_SP

!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

       ! PARTICLE sigma LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"sigma",FOUND)
!       IF(.not. FOUND) CALL FATAL_ERROR&
!            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: sigma")
       IF(FOUND)THEN
         init_with_sigma = .true. 
         CALL NC_CONNECT_PVAR(VAR,Spnt)
         CALL NC_READ_VAR(VAR)
       ENDIF
       ! PARTICLE Z(meters) LOCATION
       VAR => FIND_VAR(NC_LAGSTART,"z",FOUND)
       IF(FOUND)THEN
         init_with_zloc  = .true. 
         init_with_sigma = .false.  !zloc takes precedence over sigma
!       IF(.not. FOUND) CALL FATAL_ERROR&
!            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: z")
         CALL NC_CONNECT_PVAR(VAR,Zpnt)
         CALL NC_READ_VAR(VAR)
       ENDIF

       !Make sure either sigma or zloc is provided
       IF(.NOT.init_with_sigma .and. .NOT.init_with_zloc)then
        CALL FATAL_ERROR&
          ("SET_LAG: THE PARTICLE INPUT FILE MUST CONTAIN EITHER z OR sigma")
       ENDIF

       ! PARTICLE PATHLENGTH
       VAR => FIND_VAR(NC_LAGSTART,"pathlength",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR&
            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: pathlength")
       PLpnt => GLOBAL_PARTICLES%pathlength
       CALL NC_CONNECT_PVAR(VAR,PLpnt)
       CALL NC_READ_VAR(VAR)

       ! PARTICLE MARK (0 IN DOMAIN, 1 OUT OF DOMAIN)
!!$       VAR => FIND_VAR(NC_LAGSTART,"mark",FOUND)
!!$       IF(.not. FOUND) CALL FATAL_ERROR&
!!$            &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: mark")
!!$       Mpnt => GLOBAL_PARTICLES%mark
!!$       CALL NC_CONNECT_PVAR(VAR,Mpnt)
!!$       CALL NC_READ_VAR(VAR)

    END IF

    ! PARTICLE GROUP DATA
    VAR => FIND_VAR(NC_LAGSTART,"group",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("SET_LAG: THE PARTICLE INPUT FILE DOES NOT CONTAIN THE VARIABLE: group")
    Gpnt => GLOBAL_PARTICLES%group
    CALL NC_CONNECT_PVAR(VAR,Gpnt)
    CALL NC_READ_VAR(VAR)


    ! END TIME FOR PARTICLES
    TEND => GLOBAL_PARTICLES%tend
    VAR => FIND_VAR(NC_LAGSTART,"tend",FOUND)
    IF(FOUND) THEN
       IF(IS_VALID_FLOAT_seconds(VAR,tzone)) THEN

          ! Get adjustment for timezone
          tminus = TIME_ZONE(TZONE,status)
          IF(.not. status) call fatal_error("set lag: could not read timezone?")

          ! Allocate space for data from disks
          ALLOCATE(tendf(nlag_gl))

          CALL NC_CONNECT_PVAR(VAR,tendf)
          CALL NC_READ_VAR(VAR)

          DO I = 1,nlag_gl
             TEND(I) = SECONDS2TIME(tendf(i)) - tminus
          END DO

       ELSEIF(IS_VALID_FLOAT_days(VAR,tzone)) THEN

          ! Get adjustment for timezone
          tminus = TIME_ZONE(TZONE,status)
          IF(.not. status) call fatal_error("set lag: could not read timezone?")

          ! Allocate space for data from disks
          ALLOCATE(tendf(nlag_gl))

          CALL NC_CONNECT_PVAR(VAR,tendf)
          CALL NC_READ_VAR(VAR)

          DO I = 1,nlag_gl
             TEND(I) = DAYS2TIME(tendf(i)) - tminus
          END DO


       ELSEIF(IS_VALID_DATETIME(VAR,tzone)) THEN


          ! Allocate space for data from disks
          ALLOCATE(tendc(nlag_gl))

          CALL NC_CONNECT_PVAR(VAR,tendc)
          CALL NC_READ_VAR(VAR)

          DO I = 1,nlag_gl
             TBEG(I) = READ_DATETIME(tendc(i),DATE_FORMAT,TZONE,status)
             IF(.not.STATUS) CALL FATAL_ERROR("set lag: COULD NOT READ END TIME?")
          END DO

       ELSE
          CALL FATAL_ERROR("SET_LAG: INVALID TEND VARIABLE IN INITIAL &
               &LAG PARTICLE FILE")

       END IF

    ELSE
       CALL FATAL_ERROR("SET_LAG: MISSING VARIABLE TEND IN INITIAL &
               &LAG PARTICLE FILE")
    END IF

    ! END TIME FOR PARTICLES
    TBEG => GLOBAL_PARTICLES%tbeg
    VAR => FIND_VAR(NC_LAGSTART,"tbeg",FOUND)
    IF(FOUND) THEN
       IF(IS_VALID_FLOAT_seconds(VAR,tzone)) THEN

          ! Get adjustment for timezone
          tminus = TIME_ZONE(TZONE,status)
          IF(.not. status) call fatal_error("set lag: could not read timezone?")

          ! Allocate space for data from disks
          ALLOCATE(tbegf(nlag_gl))

          CALL NC_CONNECT_PVAR(VAR,tbegf)
          CALL NC_READ_VAR(VAR)

          DO I = 1,nlag_gl
             TBEG(I) = SECONDS2TIME(tbegf(i)) - tminus
          END DO

       ELSEIF(IS_VALID_FLOAT_days(VAR,tzone)) THEN

          ! Get adjustment for timezone
          tminus = TIME_ZONE(TZONE,status)
          IF(.not. status) call fatal_error("set lag: could not read timezone?")

          ! Allocate space for data from disks
          ALLOCATE(tbegf(nlag_gl))

          CALL NC_CONNECT_PVAR(VAR,tbegf)
          CALL NC_READ_VAR(VAR)

          DO I = 1,nlag_gl
             TBEG(I) = DAYS2TIME(tbegf(i)) - tminus
          END DO


       ELSEIF(IS_VALID_DATETIME(VAR,tzone)) THEN

          ! Allocate space for data from disks
          ALLOCATE(tbegc(nlag_gl))

          CALL NC_CONNECT_PVAR(VAR,tbegc)
          CALL NC_READ_VAR(VAR)

          DO I = 1,nlag_gl
             TBEG(I) = READ_DATETIME(tbegc(i),DATE_FORMAT,TZONE,status)
             IF(.not.STATUS) CALL FATAL_ERROR("set lag: COULD NOT READ BEGIN TIME?")
          END DO

       ELSE
          CALL FATAL_ERROR("SET_LAG: INVALID TBEG VARIABLE IN INITIAL &
               &LAG PARTICLE FILE")

       END IF
    ELSE
       CALL FATAL_ERROR("SET_LAG: MISSING VARIABLE TEND IN INITIAL &
            &LAG PARTICLE FILE")
    END IF


    !------------------------------------------------------------
    !   Adjust special FVCOM setting and add particles to linked list
    !------------------------------------------------------------
!---------------------------------------------------------------
# if !defined(SPHERICAL)
    ! ADJUST X,Y POSITION BY VXMIN AND VYMIN
    GLOBAL_PARTICLES%x(1) = GLOBAL_PARTICLES%x(1) - VXMIN
    GLOBAL_PARTICLES%x(2) = GLOBAL_PARTICLES%x(2) - VYMIN
# endif
!---------------------------------------------------------------

    DO i = 1, nlag_gl
       p => GLOBAL_PARTICLES(I)

       ! SET GLOBAL PARTICLE ID NUMBER
       p%ID = i

       ! Find particle if it is in the local domain
       p%elem = Find_Element_Containing(p%x(1),p%x(2))

       ! IF FOUND INSERT INTO LINKED LIST
       IF (p%elem == 0 .or. p%elem > N) THEN
          p%found = .false.
       ELSE
          p%found = .true.
          p%pid=myid
          call node_insert(particle_list,P)
       END IF

    end do

    !------------------------------------------------------------
    !  modify particle positions in z coordinate [zloc] 
    !  ensure [-h < z > el]
    !  if cold start, assume z = 0 means z = el if el /= 0
    !------------------------------------------------------------

    !get surface elevation and bathymetry at particles (x,y) set [el,h]
    lp  => particle_list%first%next
    do
       if(.not. associated(lp) ) exit
          lp%v%el = INTERP_ANODAL(lp%v%x(1),lp%v%x(2),lp%v%elem,EL)
          lp%v%h = INTERP_ANODAL(lp%v%x(1),lp%v%x(2),lp%v%elem,H)
       lp => lp%next
    end do

    !set z position using sigma (if only sigma supplied in init file)
    if(init_with_sigma)then
      lp  => particle_list%first%next
      do
         if(.not. associated(lp) ) exit
            lp%v%zloc = lp%v%x(3)*(lp%v%h + lp%v%el) + lp%v%el
         lp => lp%next
      end do
    endif

    !now calculate adjust z and to be inside the depth range
    lp  => particle_list%first%next
    do  
       if(.not. associated(lp) ) exit       !end of list, exit
       if(lag_coldstart .and. lp%v%zloc == 0.0_SP) lp%v%zloc = lp%v%el
       if(lp%v%zloc < -lp%v%h) lp%v%zloc = -lp%v%h
       lp => lp%next                              !set object
    end do

    !------------------------------------------------------------
    !  set sigma value of particle positions [x(3)]
    !  (only if z is supplied in init file)
    !------------------------------------------------------------

    if(init_with_zloc)then
      lp  => particle_list%first%next
      do
         if(.not. associated(lp) ) exit       !end of list, exit
         lp%v%x(3) = (lp%v%zloc - lp%v%el)/(lp%v%h + lp%v%el)
         lp => lp%next                              !set object
      end do
    endif

    !------------------------------------------------------------
    !   particle position at last time step [xn]       
    !------------------------------------------------------------
    call shift_pos_list(particle_list)

    !-------------------------------------------------------------
    !  screen reporting of particle processor distribution 
    !-------------------------------------------------------------

    !use mpi_gather to get domain distribution to master 
    allocate(tmp(nprocs))
    nlag = listsize(particle_list) !my # particles
    tmp(1) = nlag

# if defined (MULTIPROCESSOR)
    call mpi_gather(nlag,  1,mpi_integer,tmp,1,mpi_integer,0,MPI_FVCOM_GROUP,ierr)
#  endif

    !check for problems (sum of domain particles > global # particles)
    if(msr)then
       if(sum(tmp) > nlag_gl)then
          write(ipt,*)'some particles are considered to be in multiple domains'
          write(ipt,*)'nlag_gl',nlag_gl 
          write(ipt,*)'sum of particles from all domains',sum(tmp)
          write(ipt,*)'processor       # particles'
          do i=1,nprocs
             write(ipt,*)i,tmp(i)
          end do
          call pstop
       endif
    endif


    ! SETUP AND SAVE RESTART FILE
    IF(RST_ON) THEN
       CALL SETUP_LAGRESTART
       CALL output_lag_restart
    END IF

    ! SET UP AND SAVE DATA FILE
    CALL SETUP_LAGFILE
    call output_lag_file

  end subroutine set_lag

  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

  !==============================================================================|
  subroutine lag_update
    !==============================================================================|
    !  update particle positions, calculate scalar fields and particle velocities  |
    !==============================================================================|

    use all_vars
    implicit none
    integer i,ii,ierr
    type(link_node), pointer    :: lp 

    !---------------------------------------------------
    !  exit if Lagrangian tracking not active          
    !---------------------------------------------------
    if(.not.lag_particles_on) return

    !---------------------------------------------------
    !  update particle position                       
    !---------------------------------------------------
    call traject(dti,IntTime,u,uf,v,vf,wts,wtts,h,el) 

    !--------------------------------------------------
    !  write particle position data to file            
    !--------------------------------------------------
    call output_lag_file

    IF(RST_ON) CALL output_lag_restart
    !--------------------------------------------------
    !  write particle position data to screen         
    !   this will generate a lot of info, but might
    !   be useful for debugging
    !--------------------------------------------------
!    call print_list(particle_list)

    return
  end subroutine lag_update


  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

  !==============================================================================|
  subroutine traject(deltat,now,u1,u2,v1,v2,w1,w2,hl,el) 
    !==============================================================================|
    !  integrate particle position from x0 to xn using velocity fields at time     |
    !  t0 and time tn                                                              |
    !  deltat:   time step between calls to traject (usually dti)                  |
    !  time:     model time in hours
    !  u1/v1/w1: velocity field (u,v,omega) at time t0                             |
    !  u2/v2/w2: velocity field (u,v,omega) at time tn                             |
    !  hl:       bathymetry                                                        | 
    !  el:       free surface elevation at time tn                                 | 
    !==============================================================================|
    use mod_prec
    use lims
    use control, only : iint
    implicit none
    real(sp),                     intent(in)    :: deltat
    TYPE(TIME),                   intent(in)    :: now 
    real(sp), allocatable, intent(in)    :: u1(:,:),u2(:,:),v1(:,:),v2(:,:)
    real(sp), allocatable, intent(in)    :: w1(:,:),w2(:,:)
    real(sp), allocatable, intent(in)    :: hl(:),el(:) 
    !------------------------------------------------------------------------------|
    !  n    : stage counter                                                        |
    !  ns   : number of stages in explicit runga-kutta                             |
    !  chix : stage function evaluation for x-velocity                             | 
    !  pdx  : stage particle x-position                                            |
    !  ul   : stage u velocity                                                     |
    !  eps  : parameter defining depth of dry element                              |
    !  dmax : maximum sigma depth                                                  |
    !  a_rk : erk coefficients (a)                                                 |
    !  b_rk : erk coefficients (b)                                                 |
    !  c_rk : erk_coefficients (c)                                                 |
    !------------------------------------------------------------------------------|
    integer                            :: ns, status,ierr
    integer,  parameter                :: mstage = 4
    real(sp), parameter                :: eps  = 1.0e-5
    real(sp), parameter                :: dmax = -1.0_dp
    real(sp), dimension(3)             :: xtmp,Xm
    real(sp), parameter, dimension(mstage) :: a_rk = (/0.0_dp,0.5_dp,0.5_dp,1.0_dp/) 
    real(sp), parameter, dimension(mstage) :: b_rk = (/1.0_dp/6.0_dp,1.0_dp/3.0_dp, &
         1.0_dp/3.0_dp,1.0_dp/6.0_dp/) 
    real(sp), parameter, dimension(mstage) :: c_rk = (/0.0_dp,0.5_dp,0.5_dp,1.0_dp/) 
    type(link_node), pointer :: p

    ! Model velocity at RK stage 
    real(sp), allocatable       :: ul(:,:),vl(:,:),wl (:,:) 
    !------------------------------------------------------------------------------|

    !------------------------------------------------------------
    !   allocate space for RK stage velocity  
    !------------------------------------------------------------
    ALLOCATE(UL(0:NT,KB),VL(0:NT,KB),WL(0:MT,KB),stat=status)
    IF(STATUS/=0) CALL FATAL_ERROR("TRAJECT: COULD NOT ALLOCATE SPACE")


    !--set particle time step (0 if particle is not released or exceeded T_end)
    p  => particle_list%first%next
    do
       if(.not. associated(p) ) exit  !end of list, exit
       if(now >= p%v%tbeg .and. now <= p%v%tend)then
          p%v%deltat = deltat
       else
          p%v%deltat = 0.0_sp
       endif
       p => p%next                          !set object
    end do



    !--initialize stage functional evaluations
    p  => particle_list%first%next
    do
       if(.not. associated(p) ) exit  !end of list, exit
       p%v%chi = 0.0_sp
       p => p%next                          !set object
    end do


    !--loop over Runge-Kutta stages and calculate stage velocities 
    do ns=1,mstage

       !!particle position at stage n (x,y,sigma)
       IF(NS > 1)THEN
          p  => particle_list%first%next
          do
             if(.not. associated(p) ) exit  !end of list, exit
!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
             XM = a_rk(ns)*p%v%deltat*p%v%chi(:,ns-1)
             p%v%x(1) = p%v%xn(1) + XM(1)/(TPI*COSD(p%v%xn(2)) + 1.0E-6)
             p%v%x(2) = p%v%xn(2) + XM(2)/TPI
             p%v%x(3) = p%v%xn(3) + XM(3)

             IF(p%v%x(1)<0.0_SP) THEN
                p%v%x(1) = p%v%x(1) +360.0_SP
             ELSEIF(p%v%x(1)>360.0_SP) THEN
                p%v%x(1) = p%v%x(1) -360.0_SP
             END IF
             
             IF(p%v%x(2) > 90.0_SP) THEN
                p%v%x(2) = 180.0_SP - p%v%x(2)
             ELSEIF(p%v%x(2) < -90.0_SP) THEN
                p%v%x(2) = -180.0_SP - p%v%x(2)
             END IF

                
# else

             p%v%x = p%v%xn + a_rk(ns)*p%v%deltat*p%v%chi(:,ns-1)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------
             p => p%next                          !set object
          end do
       END IF


       !!determine element location 
       p  => particle_list%first%next
       do
          if(.not. associated(p) ) exit  !end of list, exit
          p%v%elem = Find_Element_Containing(p%v%x(1),p%v%x(2),p%v%elem)
          if (p%v%elem == 0) p%v%found = .false.
          p => p%next  
       end do

       CALL DELETE_NOT_FOUND(particle_list)

       !!shift particles to new domain if they have crossed boundary
       IF(PAR) call exchange_particles


       !!determine el/h at stage position to convert to sigma
       p  => particle_list%first%next
       do
          if(.not. associated(p) ) exit
          if(p%v%elem /= 0) then
             p%v%el = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%elem,EL)
             p%v%h  = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%elem,HL)
          end if
          p => p%next
       end do

       !!adjust sigma position of particle at boundaries 
       p  => particle_list%first%next
       do
          if(.not. associated(p) ) exit  !end of list, exit
          p%v%x(3) = max(p%v%x(3),-(2.0_SP+p%v%x(3))) !mirror bottom 
          p%v%x(3) = min(p%v%x(3),0.0_sp)          !don't penetrate free surf
          p => p%next                         !set object
       end do

       !!calculate velocity field for stage n using c_rk coefficients
       ul  = (1.0_sp-c_rk(ns))*u1 + c_rk(ns)*u2 
       vl  = (1.0_sp-c_rk(ns))*v1 + c_rk(ns)*v2 
       wl  = (1.0_sp-c_rk(ns))*w1 + c_rk(ns)*w2 

       p  => particle_list%first%next
       do
          if(.not. associated(p) ) exit
             p%v%u = INTERP_AZONAL(p%v%x(1),p%v%x(2),p%v%x(3),KBM1,p%v%elem,ul)
             p%v%v = INTERP_AZONAL(p%v%x(1),p%v%x(2),p%v%x(3),KBM1,p%v%elem,vl)
             p%v%w = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%x(3),KB,p%v%elem,wl)
          p => p%next
       end do

       ! THE POSITION, p%v%x has not changed since we last did this?
       ! Do we need to update EL and H?

       !!calculate u,v,w at stage ns
       p  => particle_list%first%next
       do
          if(.not. associated(p) ) exit  !end of list, exit
          p%v%chi(1,ns) = p%v%u
          p%v%chi(2,ns) = p%v%v
          p%v%chi(3,ns) = p%v%w/(p%v%h+p%v%el)         !delta_sigma/deltat = omega/d
          p => p%next                         !set object
       end do

       !!do not allow vertical motion in very shallow water
       p  => particle_list%first%next
       do
          ! THIS SHOULD NEVER HAPPEN, eps IS TOO SMALL?
          if(.not. associated(p) ) exit  !end of list, exit
          if((p%v%h + p%v%el) < eps) p%v%chi(3,:) = 0.0_sp
          p => p%next                         !set object
       end do

    end do !stage loop

    !-sum stage contributions to get updated particle positions-------------------!
    !-store accumulated pathlength
    p  => particle_list%first%next
    do
       if(.not. associated(p) ) exit  !end of list, exit

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
       xtmp = 0.0_SP
       do ns=1,mstage
          XM = p%v%deltat*p%v%chi(:,ns)*b_rk(ns)
          
          XTMP = XTMP + XM
          
          p%v%xn(1) = p%v%xn(1) + XM(1)/(TPI*COSD(p%v%xn(2)) + 1.0E-6)
          p%v%xn(2) = p%v%xn(2) + XM(2)/TPI
          p%v%xn(3) = p%v%xn(3) + XM(3)


          
          IF(p%v%xn(1)<0.0_SP) THEN
             p%v%xn(1) = p%v%xn(1) +360.0_SP
          ELSEIF(p%v%xn(1)>360.0_SP) THEN
             p%v%xn(1) = p%v%xn(1) -360.0_SP
          END IF
          
          IF(p%v%xn(2) > 90.0_SP) THEN
             p%v%xn(2) = 180.0_SP - p%v%xn(2)
          ELSEIF(p%v%xn(2) < -90.0_SP) THEN
             p%v%xn(2) = -180.0_SP - p%v%xn(2)
          END IF



       end do
       p%v%pathlength = p%v%pathlength + sqrt(XTMP(1)**2 + XTMP(2)**2)
       
# else
       xtmp = p%v%xn
       do ns=1,mstage
          p%v%xn = p%v%xn + p%v%deltat*p%v%chi(:,ns)*b_rk(ns)
       end do
       p%v%pathlength = p%v%pathlength + sqrt( (p%v%xn(1)-xtmp(1))**2 + (p%v%xn(2)-xtmp(2))**2)
!---------------------------------------------------------------
# endif
!---------------------------------------------------------------

       p => p%next                         !set object
    end do

    !-adjust sigma locs near top/bottom-------------------------------------------!
    p  => particle_list%first%next
    do
       if(.not. associated(p) ) exit  !end of list, exit
       p%v%xn(3) = max(p%v%xn(3),-(2.0+p%v%xn(3) ))      !mirror bottom 
       p%v%xn(3) = min(p%v%xn(3),0.0_sp)               !don't penetrate free surf
       p => p%next                         !set object
    end do

    !-update x location and find the element -------------------------------------!
    p  => particle_list%first%next
    do
       if(.not. associated(p) ) exit  !end of list, exit
       p%v%x = p%v%xn
       p%v%elem = Find_Element_Containing(p%v%x(1),p%v%x(2),p%v%elem)
       if (p%v%elem == 0) p%v%found = .false.
        p => p%next  
    end do
    
    CALL DELETE_NOT_FOUND(particle_list)

    !--evaluate bathymetry and free surface height at updated particle position----!
    p  => particle_list%first%next
    do
       if(.not. associated(p) ) exit
       p%v%el = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%elem,EL)
       p%v%h = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%elem,HL)
       p%v%zloc = (p%v%x(3)*(p%v%h+p%v%el) + p%v%el) !particle depth (m)
       p => p%next
    end do

    DEALLOCATE(UL,VL,WL)

    return
  end subroutine traject
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

  !==============================================================================|
  !==============================================================================
  ! isintriangle and sameside were moved to mod_utils.F
  !==============================================================================

  !==============================================================================
  ! All FHE routines were moved to mod_utils.F
  !==============================================================================

  !==============================================================================
  ! All interp routines were moved to mod_utils.F
  !==============================================================================

  !==============================================================================|
  ! exchange particles across processor domains                                  |
  !==============================================================================|
  subroutine exchange_particles
    !------------------------------------------------------------------------------|
    use lims,    only : myid ,nprocs, n
    use mod_par, only : el_pid, he_own, he_lst, elid
    use control, only : iint
    use all_vars, only : vx,vy,nv
    implicit none

# if defined (MULTIPROCESSOR)

    integer               ::ireqr(nprocs),ireqs(nprocs)
    !  real(sp), allocatable :: rbuf(:),sbuf(:)
    TYPE(particle), allocatable :: rbuf(:),sbuf(:)
    integer   stat(mpi_status_size),istatr(mpi_status_size,nprocs),ierr,j,n1,n2,ncnt
    integer   i,ifrom,ito,istag,irtag,trcv,tsnd,nvars,lbuf,nmsg,indx,lproc,nsze,rp
    integer   ibuf,recv_buf_size,send_buf_size
    integer   tri_local,tri_global

    integer, allocatable :: nrcv(:)
    integer, allocatable :: nsnd(:)
    type(particle)                :: p
    type(link_node), pointer      :: lp

    logical check
    !------------------------------------------------------------------------------

    !--------------------------------------------
    !find particles in halo
    !set processor id to neighboring processor
    !we will send particle to that processor
    !--------------------------------------------
    lp  => particle_list%first%next
    do
       if(.not. associated(lp) ) exit 
       tri_local = lp%v%elem
       if(tri_local > n)then
          lp%v%pid  = he_own(tri_local-n) 
          lp%v%elem = he_lst(tri_local-n)   
       endif
       lp => lp%next
    end do

    !-------------------------------------------
    !determine # to send to each processor 
    !-------------------------------------------
    allocate(nsnd(nprocs)) ; nsnd = 0
    allocate(nrcv(nprocs)) ; nrcv = 0
    lp  => particle_list%first%next
    do
       if(.not. associated(lp) ) exit  
       if(lp%v%pid /= myid) nsnd(lp%v%pid) = nsnd(lp%v%pid) + 1
       lp => lp%next
    end do


    !-------------------------------------------
    !determine # to recv from each processor 
    !-------------------------------------------
    call recv_count(nsnd,nrcv)

    !-------------------------------------------
    !allocate buffer space      
    !-------------------------------------------
    !  send_buf_size = max(par_type_size*sum(nsnd) , 1)
    !  recv_buf_size = max(par_type_size*sum(nrcv) , 1)
    !  allocate(sbuf(send_buf_size)) ; sbuf = 0.0_sp  
    !  allocate(rbuf(recv_buf_size)) ; rbuf = 0.0_sp  
    call NEW_PARTICLES(sbuf,max(sum(nsnd) , 1))
    call NEW_PARTICLES(rbuf,max(sum(nrcv) , 1))


    !-------------------------------------------
    !Post non-blocking receives from neighbors
    !-------------------------------------------

    trcv = 0
    rp   = 1
    do i=1,nprocs
       if(nrcv(i) > 0)then
          ifrom = i-1
          irtag = i*1000
          trcv  = trcv + 1
          lbuf  = nrcv(i) 
          call mpi_irecv(rbuf(rp),lbuf,mpi_particle,ifrom,irtag,mpi_fvcom_group,ireqr(trcv),ierr)
       end if
       rp = rp + nrcv(i)
    end do

    !-------------------------------------------
    !send data to neighbors
    !-------------------------------------------
    tsnd = 0
    ibuf = 1
    do i=1,nprocs
       lbuf = nsnd(i) 
       if(lbuf == 0 .or. myid == i)cycle

       !load send array
       ncnt = ibuf
       lp   => particle_list%first%next
       do
          if(.not. associated(lp) ) exit  !end of list, exit
          if(lp%v%pid == i)then
             sbuf(ibuf)=lp%v
             ibuf = ibuf +1
          endif
          lp => lp%next                          !set object
       end do

       !send
       tsnd  = tsnd + 1
       ito   = i-1
       istag = myid*1000
       call mpi_isend(sbuf(ncnt),lbuf,mpi_particle,ito,istag,mpi_fvcom_group,ireqs(tsnd),ierr)
    end do

    !--------------------------------------------
    !loop over procs until a message is received
    !--------------------------------------------
    do nmsg = 1,trcv 
       rp = 1
       call mpi_waitany(trcv,ireqr,indx,stat,ierr)
       lproc = stat(mpi_source) +1 
       do i=1,lproc-1
          rp = rp + nrcv(i)
       end do
       do i=1,nrcv(lproc)
          p = rbuf(rp)
          p%elem = elid(p%elem) 
          call node_insert(particle_list,p)
          rp = rp+1
       end do
    end do

    !--------------------------------------------
    !wait for completion of non-blocking sends  
    !--------------------------------------------

    call mpi_waitall(tsnd,ireqs,istatr,ierr)
    call mpi_waitall(trcv,ireqr,istatr,ierr)
    deallocate(sbuf)
    deallocate(rbuf)

    CALL DELETE_NOT_MINE(particle_list,MYID)

!!$    lp  => particle_list%first%next
!!$    do
!!$       if(.not. associated(lp) ) exit  !end of list, exit
!!$       if(lp%v%pid /= myid)then
!!$          call node_delete(particle_list,lp%v,check)
!!$          lp => particle_list%first%next
!!$       else
!!$          lp => lp%next
!!$       endif
!!$    end do


    !--------------------------------------------
    !reset home element
    !--------------------------------------------
!!$    lp  => particle_list%first%next
!!$    do
!!$       if(.not. associated(lp) ) exit  !end of list, exit
!!$
!!$       if(.not. isintriangle(lp%v%elem,lp%v%x(1),lp%v%x(2)))then
!!$          call warning("Exchange particle found particle in the wrong cell?")
!!$       end if
!!$       lp%v%elem = Find_Element_Containing(lp%v%x(1),lp%v%x(2),lp%v%elem)
!!$       if (lp%v%elem == 0) lp%v%found = .false.
!!$       lp => lp%next  
!!$    end do

# endif

  end subroutine exchange_particles


  !==============================================================================|
  ! count particle receives across processor boundaries                          |
  !==============================================================================|
# if defined (MULTIPROCESSOR)
  subroutine recv_count(nsnd,nrcv)
    !------------------------------------------------------------------------------|
    use mod_par, only : nhe, he_own
    use lims,    only : myid ,nprocs, n
    implicit none
    integer, intent(in ) :: nsnd(nprocs)
    integer, intent(out) :: nrcv(nprocs)

    integer               ::ireqr(nprocs),ireqs(nprocs)
    real(sp), allocatable :: rbuf(:),sbuf(:)
    integer   stat(mpi_status_size),istatr(mpi_status_size,nprocs),ierr,j,n1,n2,ncnt
    integer   i,ifrom,ito,istag,irtag,trcv,tsnd,nvars,lbuf,nmsg,indx,lproc,nsze,rp
    integer   ibuf,recv_buf_size,send_buf_size
    integer   tri_local,tri_global

    type(particle)                :: p
    type(particle),allocatable    :: p_arr(:)
    type(link_node), pointer      :: lp
    logical :: first_entry = .true.

    logical check
    !------------------------------------------------------------------------------

    !---------------------------------------------
    !if first entry, determine adjacent processors
    !--------------------------------------------
    if(first_entry)then
       allocate(ney_list(nprocs)) ; ney_list = 0
       do i=1,nhe
          ney_list(he_own(i)) = 1
       end do
       first_entry = .false.
    endif

    !-------------------------------------------
    !Post non-blocking receives from neighbors
    !-------------------------------------------

    trcv = 0
    do i=1,nprocs
       if(ney_list(i) > 0)then
          ifrom = i-1
          irtag = i*1000
          trcv  = trcv + 1
          call mpi_irecv(nrcv(i),1,mpi_integer,ifrom,irtag,mpi_fvcom_group,ireqr(trcv),ierr)
       end if
    end do

    !-------------------------------------------
    !send data to neighbors
    !-------------------------------------------
    tsnd = 0
    do i=1,nprocs
       if(ney_list(i) > 0)then
          tsnd  = tsnd + 1
          ito   = i-1
          istag = myid*1000
          call mpi_isend(nsnd(i),1,mpi_integer,ito,istag,mpi_fvcom_group,ireqs(tsnd),ierr)
       endif
    end do

    !--------------------------------------------
    !loop over procs until a message is received
    !--------------------------------------------
    do nmsg = 1,trcv 
       call mpi_waitany(trcv,ireqr,indx,stat,ierr)
    end do

    !--------------------------------------------
    !wait for completion of non-blocking sends  
    !--------------------------------------------

    call mpi_waitall(tsnd,ireqs,istatr,ierr)


  end subroutine recv_count
# endif

  !==============================================================================|
  ! collect particles to global list in master only 
  !     needed to dump particles to file
  !==============================================================================|
  subroutine collect_particles(PID,RECVID,P_glb)
    USE ALL_VARS, only : vxmin,vymin
    implicit none

    INTEGER, INTENT(IN) :: PID, RECVID
    type(particle), allocatable    :: P_glb(:)
    type(particle), allocatable    :: P_lcl(:)
    integer, allocatable :: tmp(:) 

# if defined (MULTIPROCESSOR)
    integer   stat(mpi_status_size)
# endif
    integer :: ierr,isnd,ircv,ip,nsze,ibuf,i
    type(particle)                :: p
    type(link_node), pointer      :: lp
    INTEGER :: sender,recver,tag


    IF(DBG_SET(DBG_SBR)) write(ipt,*) "start particles_collect: PID:",PID,"; RECVID:",RECVID


    IF(MYID == RECVID .and. .not. allocated(p_glb)) CALL FATAL_ERROR &
         &('particle_collect: The collecting processor must have allocated p_glb argument!')


    !count number of particles in local domain => nlag
    nlag = listsize(particle_list)

    !send to master
    allocate(tmp(nprocs))
    tmp(1) = nlag

# if defined (MULTIPROCESSOR)
    call mpi_allgather(nlag,1,mpi_integer,tmp,1,mpi_integer,mpi_fvcom_group,ierr)
# endif

    !report processor particle distribution to screen
    if(DBG_SET(DBG_IO))then
       write(ipt,*)'   particle/processor distribution'
       write(ipt,*)'processor       # particles'
       do i=1,nprocs
          write(ipt,*)i,tmp(i)
       end do
       write(ipt,*)'max: ',maxval(tmp),' total: ',sum(tmp),' min ',minval(tmp)
    endif

    !loop through other processors, receive particles, add to linked list
    
    !----------------------------------------------------
    !send data to master
    !----------------------------------------------------
    IF(MYID /= RECVID)THEN
       !allocate send space
       call NEW_PARTICLES(P_LCL,nlag)
       
       !load send array
       ibuf = 1
       lp  => particle_list%first%next
       do
          if(.not. associated(lp) ) exit  !end of list, exit
          p_lcl(ibuf) = lp%v
          ibuf=ibuf+1
          lp => lp%next                          !set object
       end do
       
       !send to master processor
       recver = msrid-1
       tag = myid+100
# if defined (MULTIPROCESSOR)
       call mpi_send(p_lcl,nlag,mpi_particle,recver,tag,mpi_fvcom_group,ierr)
# endif
       deallocate(P_lcl)
       
       
    ELSE
       
       !----------------------------------------------------
       !recv data and add to linked list
       !----------------------------------------------------
       
       ! SET P_GLB%ELEM to zero - find which particles are still in domain! 
       P_GLB%ELEM=0

       do ip=1,nprocs
          
          IF(MYID == IP) THEN
             lp  => particle_list%first%next
             do
                if(.not. associated(lp) ) exit  !end of list, exit
                p_glb(lp%v%id) = lp%v
                lp => lp%next                          !set object
             end do
          ELSE
             call NEW_PARTICLES(p_lcl,tmp(ip))
             
             tag = ip+100
             sender = ip -1
# if defined (MULTIPROCESSOR)
             call mpi_recv(p_lcl,tmp(ip),mpi_particle,sender,tag,mpi_fvcom_group,stat,ierr)
# endif
             do i=1,tmp(ip)
                p_glb(p_lcl(i)%id) = p_lcl(i)
             end do
             deallocate(P_lcl)
          END IF
       end do

!---------------------------------------------------------------
# if defined(SPHERICAL)
!---------------------------------------------------------------
       
       IF(USE_PROJ) CALL DEGREES2METERS(LonPnt,LatPnt,PROJECTION_REFERENCE,Xpnt,Ypnt,nlag_gl)


# else
       ! RETURN TO GLOBAL COORDINATES
       WHERE(P_GLB%ELEM /= 0)
          P_glb%x(1) = P_glb%x(1) + vxmin
          P_glb%x(2) = P_glb%x(2) + vymin
       END WHERE

       IF(USE_PROJ) CALL METERS2DEGREES(XPnt,YPnt,PROJECTION_REFERENCE,LonPnt,LatPnt,nlag_gl)


!---------------------------------------------------------------
# endif
!---------------------------------------------------------------


    END IF
 
    deallocate(tmp)
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "end particles_collect"
    
    
  end subroutine collect_particles

  !==============================================================================|
  ! output_lag  
  !==============================================================================|
  subroutine output_lag_file
    USE mod_ncdio, only : update_iodata
    implicit none

    IF(NC_LAG%FTIME%NEXT_IO .GT. IntTime) RETURN 

    NC_LAG%FTIME%NEXT_IO = NC_LAG%FTIME%NEXT_IO + NC_LAG%FTIME%INTERVAL
    NC_LAG%FTIME%NEXT_STKCNT =NC_LAG%FTIME%NEXT_STKCNT + 1   
    
    CALL ADD_SCALARS

    ! ANY PROCESSOR CAN DO OUTPUT - SET IT AT THE TOP OF THE MODULE
    CALL COLLECT_PARTICLES(MYID,LAG_OUT_PROC,GLOBAL_PARTICLES)
    IF(MYID /= LAG_OUT_PROC) RETURN

    CALL UPDATE_IODATA(NC_LAG,IntTime)

    CALL NC_WRITE_FILE(NC_LAG,LOCAL_ONLY=.true.)

  end subroutine output_lag_file

  !==============================================================================|
  ! output_lag_restart 
  !==============================================================================|
  subroutine output_lag_restart
    USE mod_ncdio, only : update_iodata
    implicit none

    IF(NC_LAGRESTART%FTIME%NEXT_IO .GT. IntTime) RETURN 

    NC_LAGRESTART%FTIME%NEXT_IO     = NC_LAGRESTART%FTIME%NEXT_IO + NC_LAGRESTART%FTIME%INTERVAL
    NC_LAGRESTART%FTIME%NEXT_STKCNT = NC_LAGRESTART%FTIME%NEXT_STKCNT + 1   

    ! ANY PROCESSOR CAN DO OUTPUT - SET IT AT THE TOP OF THE MODULE
    CALL COLLECT_PARTICLES(MYID,LAG_OUT_PROC,GLOBAL_PARTICLES)
    IF(MYID /= LAG_OUT_PROC) RETURN

    CALL UPDATE_IODATA(NC_LAGRESTART,IntTime)

    CALL NC_WRITE_FILE(NC_LAGRESTART,LOCAL_ONLY=.true.)

  end subroutine output_lag_restart


  !==============================================================================|
  ! SETUP_LAG FILE 
  !==============================================================================|
  subroutine SETUP_LAGFILE
!    USE MOD_NCDIO, only : DIM_DateStrLen
    use all_vars
    use mod_clock, only : get_timestamp
    implicit none
    character(len=100)    :: timestamp, temp, netcdf_convention,fname
    INTEGER :: I

    TYPE(NCFILE), pointer :: NCF
    TYPE(NCATT), pointer :: ATT
    TYPE(NCVAR), pointer :: VAR

    TYPE(NCDIM), POINTER :: DIM_nlag
    TYPE(NCDIM), POINTER :: DIM_DateStrLenLag
    TYPE(NCDIM), POINTER :: DIM_lagtime

    DIM_nlag    => NC_MAKE_DIM(name='nlag',len=nlag_gl) 
    DIM_lagtime => NC_MAKE_DIM(name='time',len=NF90_UNLIMITED) 
    DIM_DateStrLenLag => NC_MAKE_DIM(name='DateStrLen',len=DateStrLen)
    NCF => NEW_FILE()


    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='title',values="FVCOM PARTICLE TRACKING:"//TRIM(CASE_TITLE)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='institution',values=trim(institution)) 
    NCF => ADD(NCF,ATT)

#  if defined (SPHERICAL)
    ATT => NC_MAKE_ATT(name='CoordinateSystem',values="GeoReferenced" ) 
    NCF => ADD(NCF,ATT)
#endif

    ATT => NC_MAKE_ATT(name='source',values="Particle"//TRIM(FVCOM_VERSION))
    NCF => ADD(NCF,ATT)

    call get_timestamp(temp)
    timestamp = 'model started at: '//trim(temp)

    ATT => NC_MAKE_ATT(name='history',values=trim(timestamp)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='references',values=trim(fvcom_website)) 
    NCF => ADD(NCF,ATT)

    netcdf_convention = 'CF-1.0'
    ATT => NC_MAKE_ATT(name='Conventions',values=trim(netcdf_convention)) 
    NCF => ADD(NCF,ATT)


    ! ADD VARIABLES
    ! time
    VAR => FLOAT_TIME_OBJECT &
         &(USE_MJD=use_real_world_time, &
         & DIM=DIM_LAGTIME)

    NCF  => ADD(NCF,VAR)


    ! Itime
    VAR  => ITIME_OBJECT &
         &(Use_MJD=use_real_world_time, &
         & DIM=DIM_LAGTIME)

    NCF  => ADD(NCF,VAR)

    ! Itime2
    VAR => ITIME2_OBJECT &
         &(Use_MJD=use_real_world_time, &
         & DIM=DIM_LAGTIME)

    NCF => ADD(NCF,VAR)

    IF (use_real_world_time) THEN

       VAR => DATETIME_OBJECT &
            &(DIMSTR=DIM_DateStrLenLag,&
            & DIMTIME=DIM_LAGTIME,&
            TIMEZONE=TIMEZONE)

       NCF  => ADD(NCF,VAR)
    END IF

    ! X
    VAR => NC_MAKE_PVAR(name='x',values=Xpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle x position') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! Y
    VAR => NC_MAKE_PVAR(name='y',values=Ypnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle y position') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! LON
    VAR => NC_MAKE_PVAR(name='lon',values=LONpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! LAT
    VAR => NC_MAKE_PVAR(name='lat',values=LATpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! Z
    VAR => NC_MAKE_PVAR(name='z',values=Zpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle z position') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! Sigma
    VAR => NC_MAKE_PVAR(name='sigma',values=Spnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle sigma position') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! ADD ANY SCALARS SPECIFIED IN THE RUNFILE
    DO I = 1,size(lag_scalars)
       SELECT CASE(LAG_SCALARS(I))
       CASE("none","NONE")
          exit
       CASE("temp","TEMP","temperature","TEMPERATURE")
          
          TMPpnt => GLOBAL_PARTICLES%s(I)

          ! TEMP
          VAR => NC_MAKE_PVAR(name='temp',values=TMPpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='particle temperature') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)
          
       CASE("salinity","SALT","SALINITY","salt")
          
          SALpnt => GLOBAL_PARTICLES%s(I)

          ! SALT
          VAR => NC_MAKE_PVAR(name='salinity',values=SALpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='particle salinity') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)

       CASE("KH","kh")
          
          KHpnt => GLOBAL_PARTICLES%s(I)

          ! Scalar eddy viscosity
          VAR => NC_MAKE_PVAR(name='kh',values=KHpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='Particle Turbulent Eddy Viscosity For Scalars') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)
          

       CASE("KM","km")
          
          KMpnt => GLOBAL_PARTICLES%s(I)

          ! Scalar eddy viscosity
          VAR => NC_MAKE_PVAR(name='km',values=KMpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='Particle Turbulent Eddy Viscosity For Momentum') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)

       CASE("DENSITY","density")
          
          DENpnt => GLOBAL_PARTICLES%s(I)

          ! Scalar eddy viscosity
          VAR => NC_MAKE_PVAR(name='density',values=DENpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='Particle Density') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='kg/m 3') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)

       CASE("U","u")
          
          Upnt => GLOBAL_PARTICLES%s(I)

          ! Scalar eddy viscosity
          VAR => NC_MAKE_PVAR(name='u',values=Upnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Particle Velocity') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)

       CASE("V","v")
          
          Vpnt => GLOBAL_PARTICLES%s(I)

          ! Scalar eddy viscosity
          VAR => NC_MAKE_PVAR(name='v',values=Vpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='Northward Particle Velocity') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)

       CASE("W","w")
          
          Wpnt => GLOBAL_PARTICLES%s(I)

          ! Scalar eddy viscosity
          VAR => NC_MAKE_PVAR(name='w',values=Wpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='Upward Particle Velocity (z)') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)

       CASE("OMEGA","omega")
          
          Opnt => GLOBAL_PARTICLES%s(I)

          ! Scalar eddy viscosity
          VAR => NC_MAKE_PVAR(name='omega',values=Opnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='Particle Sigma Velocity') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='s-1') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)

       END SELECT
    END DO

    ! IF THIS IS A CRASH RESTART RECONNECT OLD OUTPUT
    IF(STARTUP_TYPE == STARTUP_TYPE_CRASHRESTART) THEN

       IF(len_trim(LAG_OUT_FILE)==0) call fatal_error&
            &("LAG_OUT_FILE is not set in the name list!")

       fname = trim(OUTPUT_DIR)//trim(LAG_OUT_FILE)
       NCF%fname=fname
       NCF%CONNECTED = .TRUE.

       NCF%FTIME=>NEW_FTIME()
       NCF%FTIME%NEXT_IO = lag_next
       NCF%FTIME%NEXT_STKCNT =LAG_STACK
       NCF%FTIME%INTERVAL = lag_interval
       
    ELSE
       IF(len_trim(LAG_OUT_FILE)==0) call fatal_error&
            &("LAG_OUT_FILE is not set in the name list!")

       fname = trim(OUTPUT_DIR)//trim(LAG_OUT_FILE)
       NCF%fname=fname
       NCF%FTIME=>NEW_FTIME()


       NCF%FTIME%NEXT_IO = lag_next
       NCF%FTIME%NEXT_STKCNT =0
       NCF%FTIME%INTERVAL = lag_interval

       ! ONLY THE OUTPUT PROC WRITE THE FILES, THE OTHERS BUILD IT AS
       ! A PLACE HOLDER


       IF(MYID == LAG_OUT_PROC) THEN
          CALL NC_WRITE_FILE(NCF,LOCAL_ONLY=.true.)
       ELSE
          NCF%CONNECTED = .TRUE.
       END IF
       
       
    END IF

    NC_LAG => NCF

  end subroutine SETUP_LAGFILE

  !==============================================================================|
  ! SETUP_LAGRESTART FILE 
  !==============================================================================|
  subroutine SETUP_LAGRESTART
!    USE MOD_NCDIO, only : DIM_DateStrLen
    USE MOD_INPUT, only : NC_RST
    USE MOD_SET_TIME, only : GET_FIRST_OUTPUT_TIME
    use all_vars
    use mod_clock, only : get_timestamp
    implicit none
    character(len=100)    :: timestamp, temp, netcdf_convention,fname

    LOGICAL FOUND
    TYPE(NCFILE), pointer :: NCF
    TYPE(NCATT), pointer :: ATT

    TYPE(NCVAR), pointer :: VAR,VAR_T

    TYPE(NCDIM), POINTER :: DIM_nlag
    TYPE(NCDIM), POINTER :: DIM_lagtime
    TYPE(NCDIM), POINTER :: DIM_DateStrLenLag

    ! IF THIS IS A CRASH RESTART RECONNECT OLD OUTPUT
    IF(STARTUP_TYPE == STARTUP_TYPE_CRASHRESTART) THEN

       NC_LAGRESTART => NC_LAGSTART

       DIM_lagtime => FIND_DIM(NC_LAGRESTART,'time',FOUND)
       IF(.NOT. FOUND) CALL FATAL_ERROR&
            & ('CAN NOT CRASH RESTART PARTICLES FROM A FILE THAT DOES N&
            &OT HAVE THE TIME DIMENSION?')

       NC_LAGRESTART%FTIME%NEXT_STKCNT = DIM_lagtime%DIM
       NC_LAGRESTART%FTIME%NEXT_IO     = NC_RST%FTIME%NEXT_IO
       NC_LAGRESTART%FTIME%INTERVAL    = NC_RST%FTIME%INTERVAL
       
       RETURN
    END IF

    ! MAKE A NEW RESTART FILE
    DIM_nlag    => NC_MAKE_DIM(name='nlag',len=nlag_gl) 
    DIM_lagtime => NC_MAKE_DIM(name='time',len=NF90_UNLIMITED) 
    DIM_DateStrLenLag => NC_MAKE_DIM(name='DateStrLen',len=DateStrLen)

    NCF => NEW_FILE()


    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='title',values="FVCOM PARTICLE TRACKING:"//TRIM(CASE_TITLE)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='institution',values=trim(institution)) 
    NCF => ADD(NCF,ATT)

#  if defined (SPHERICAL)
    ATT => NC_MAKE_ATT(name='CoordinateSystem',values="GeoReferenced" ) 
    NCF => ADD(NCF,ATT)
# else
    ATT => NC_MAKE_ATT(name='CoordinateSystem',values="Cartesian" ) 
    NCF => ADD(NCF,ATT)
#endif

    ATT => NC_MAKE_ATT(name='source',values="Particle"//TRIM(FVCOM_VERSION))
    NCF => ADD(NCF,ATT)

    call get_timestamp(temp)
    timestamp = 'model started at: '//trim(temp)

    ATT => NC_MAKE_ATT(name='history',values=trim(timestamp)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='references',values=trim(fvcom_website)) 
    NCF => ADD(NCF,ATT)

    netcdf_convention = 'CF-1.0'
    ATT => NC_MAKE_ATT(name='Conventions',values=trim(netcdf_convention)) 
    NCF => ADD(NCF,ATT)

        ! ADD VARIABLES
    ! time
    VAR => FLOAT_TIME_OBJECT &
         &(USE_MJD=use_real_world_time, &
         & DIM=DIM_LAGTIME)

    NCF  => ADD(NCF,VAR)


    ! Itime
    VAR  => ITIME_OBJECT &
         &(Use_MJD=use_real_world_time, &
         & DIM=DIM_LAGTIME)

    NCF  => ADD(NCF,VAR)

    ! Itime2
    VAR => ITIME2_OBJECT &
         &(Use_MJD=use_real_world_time, &
         & DIM=DIM_LAGTIME)

    NCF => ADD(NCF,VAR)

    IF (use_real_world_time) THEN

       VAR => DATETIME_OBJECT &
            &(DIMSTR=DIM_DateStrLenLag,&
            & DIMTIME=DIM_LAGTIME,&
            TIMEZONE=TIMEZONE)

       NCF  => ADD(NCF,VAR)
    END IF

    ! X
    VAR => NC_MAKE_PVAR(name='x',values=Xpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle x position') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! Y
    VAR => NC_MAKE_PVAR(name='y',values=Ypnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle y position') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! LON
    VAR => NC_MAKE_PVAR(name='lon',values=LONpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! LAT
    VAR => NC_MAKE_PVAR(name='lat',values=LATpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)



    ! Z
    VAR => NC_MAKE_PVAR(name='z',values=Zpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle z position') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! PathLength
    VAR => NC_MAKE_PVAR(name='pathlength',values=PLpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle integrated path length') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

!!$    ! Particle Mark
!!$    VAR => NC_MAKE_PVAR(name='mark',values=Mpnt,DIM1=DIM_nlag,DIM2=DIM_lagtime)
!!$
!!$    ATT  => NC_MAKE_ATT(name='long_name',values='particle marker (0=in domain)') 
!!$    VAR  => ADD(VAR,ATT)
!!$
!!$    ATT  => NC_MAKE_ATT(name='units',values='-') 
!!$    VAR  => ADD(VAR,ATT)
!!$
!!$    NCF  => ADD(NCF,VAR)

    ! Group
    VAR => NC_MAKE_PVAR(name='group',values=Gpnt,DIM1=DIM_nlag)

    ATT  => NC_MAKE_ATT(name='long_name',values='particle group') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='-') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(.not.ASSOCIATED(NC_LAGSTART)) CALL FATAL_ERROR&
         & ("MOD_LAG: SETUP_LAGRESTART; NC_LAGSTART IS NOT ASSOCIATED?")

    ! TBEG - use the same type as the startup file!
    VAR_T => FIND_VAR(NC_LAGSTART,"tbeg",FOUND)
    VAR => COPY_VAR(VAR_T)

    NCF  => ADD(NCF,VAR)

    ! TEND - use the same type as the startup file!
    VAR_T => FIND_VAR(NC_LAGSTART,"tend",FOUND)
    VAR => COPY_VAR(VAR_T)

    NCF  => ADD(NCF,VAR)



    ! SETUP FILE NAME AND TIME
    IF(len_trim(LAG_RESTART_FILE)==0) call fatal_error&
         &("LAG_RESTART_FILE is not set in the name list!")
    fname = trim(OUTPUT_DIR)//trim(LAG_RESTART_FILE)
    NCF%fname=fname
    NCF%FTIME=>NEW_FTIME()

    NCF%FTIME%NEXT_STKCNT =0
    CALL GET_FIRST_OUTPUT_TIME(TRIM(RST_FIRST_OUT),NCF%FTIME%NEXT_IO)
    NCF%FTIME%INTERVAL = NC_RST%FTIME%INTERVAL    

    ! ONLY THE OUTPUT PROC WRITE THE FILES, THE OTHERS BUILD IT AS
    ! A PLACE HOLDER
    IF(MYID == LAG_OUT_PROC) THEN
       CALL NC_WRITE_FILE(NCF,LOCAL_ONLY=.true.)
    ELSE
       NCF%CONNECTED = .TRUE.
    END IF
    

    ! Move the pointer
    NC_LAGRESTART => NCF

    
  end subroutine SETUP_LAGRESTART



  subroutine add_scalars
    USE ALL_VARS
    implicit none
    integer:: i
    type(link_node), pointer :: p

    p  => particle_list%first%next
    outer : do
       if(.not. associated(p) ) exit  !end of list, exit
       
       
       DO I = 1,size(lag_scalars)
          SELECT CASE(LAG_SCALARS(I))
          CASE("none","NONE")
             exit outer
          CASE("temp","TEMP","temperature","TEMPERATURE")

             p%v%s(I) = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%x(3),KBM1,p%v%elem,T1)
             
          CASE("salinity","SALT","SALINITY","salt")

             p%v%s(I) = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%x(3),KBM1,p%v%elem,S1)
                          
          CASE("KH","kh")

             p%v%s(I) = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%x(3),KB,p%v%elem,KH)
                          
          CASE("KM","km")

             p%v%s(I) = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%x(3),KB,p%v%elem,KM)

          CASE("DENSITY","density")

             p%v%s(I) = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%x(3),KBM1,p%v%elem,rho1)

          CASE("U","u")

             p%v%s(I) = INTERP_AZONAL(p%v%x(1),p%v%x(2),p%v%x(3),KBM1,p%v%elem,U)

          CASE("V","v")

             p%v%s(I) = INTERP_AZONAL(p%v%x(1),p%v%x(2),p%v%x(3),KBM1,p%v%elem,V)

          CASE("OMEGA","omega")

             p%v%s(I) = INTERP_AZONAL(p%v%x(1),p%v%x(2),p%v%x(3),KBM1,p%v%elem,W)

          CASE("W","w")

             p%v%s(I) = INTERP_ANODAL(p%v%x(1),p%v%x(2),p%v%x(3),KB,p%v%elem,wts)
             
          END SELECT
       END DO
       
       p => p%next                          !set object
    end do outer
    

  end subroutine add_scalars

# endif
END module mod_lag
