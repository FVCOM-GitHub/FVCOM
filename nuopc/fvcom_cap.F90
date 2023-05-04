!>
!! @mainpage FVCOM NUOPC Cap
!! @author Jianhua Qi (jqi@umassd.edu) basded on Saeed Moghimi (moghimis@gmail.com)
!! @date 05/25/2020 Original documentation
!------------------------------------------------------

#define MULTIPROCESSOR
#define SPHERICAL
#define WET_DRY
#define TWO_D_MODEL___disabled
#define WAVE_ROLLER___disabled

module fvcom_cap

  !-----------------------------------------------------------------------------
  ! FVCOM Component.
  !-----------------------------------------------------------------------------

  use mpi
  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS        => SetServices,       &
    model_label_SetClock    => label_SetClock,    &
    model_label_CheckImport => label_CheckImport, &    
    model_label_Advance     => label_Advance,     &
    model_label_Finalize    => label_Finalize

!JQI  USE MESSENGER, ONLY  : MPI_COMM_FVCOM,UPDATER

  use NUOPC_FVCOM, only : NUOPC_FVCOM_init
  use NUOPC_FVCOM, only : NUOPC_FVCOM_run
  use NUOPC_FVCOM, only : NUOPC_FVCOM_Final
  use LIMS       , only : MGL,MT,M,NGL,NT,N,KB,KBM1,MYID,NPROCS    
  use ALL_VARS   , only : NTVE,NBVE,EL,U,V,VX,YC
  use ALL_VARS   , only : WAVESTRX_2D, WAVESTRY_2D,WAVESTRX_3D,WAVESTRY_3D,nv
  use ALL_VARS   , only : U_STOKES_2D, V_STOKES_2D,U_STOKES_3D,V_STOKES_3D
  use ALL_VARS   , only : PAR, INTSTEP_SECONDS
  use MOD_PAR    , only : NBN,BN_MLT,BN_LOC,BNC,NODE_MATCH
  use MOD_PAR    , only : NC,EC,AEXCHANGE
  use MOD_PAR    , only : NMAP,MPI_F,MPI_FVCOM_GROUP,NGID_X,ACOLLECT
  use LIMS       , only : MSRID
  use mod_prec   , only : SP

  use mod_wd     , only : ISWETC,ISWETCT
  use mod_spherical, only : TPI,DEG2RAD,DELTUY

  use fvcom_mod, only: FVCOM_WHS,  FVCOM_WLEN, FVCOM_WDIR

  use fvcom_mod, only: NUOPC4WAV

  use fvcom_mod, only: meshdata
  use fvcom_mod, only: create_parallel_esmf_mesh_from_meshdata
  use fvcom_mod, only: extract_parallel_data_from_mesh
  
  implicit none

  private
  public SetServices
  public fvcom_name

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: unit
    logical           :: assoc    ! is the farrayPtr associated with internal data
    logical           :: connected
    real(ESMF_KIND_R8), dimension(:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter :: fldsMax = 100
  integer :: fldsToFvcom_num = 0
  type (fld_list_type) :: fldsToFvcom(fldsMax)
  integer :: fldsFrFvcom_num = 0
  type (fld_list_type) :: fldsFrFvcom(fldsMax)

  type(meshdata),save  :: mdataIn, mdataOut
  
  character(len=2048):: info
  integer :: dbrc     ! temporary debug rc value

  !to test field halo update.
!JQI  type (ESMF_RouteHandle), save :: ATM_HaloRouteHandel
  
  !real(ESMF_KIND_R8)      :: WaveCouplingIntervalSec, WindCouplingIntervalSec  !in seconds
  !type(ESMF_TimeInterval) :: WaveCouplingInterval, WindCouplingInterval

!JQI2021  logical, save            :: first_exchange = .true.
!JQI2021  integer, save            :: iunit_log = 10000


!JQI2021  real,parameter :: wave_force_limmit = 0.05

! The FVCOM domain name
  character(len=20) :: fvcom_name
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  !> NUOPC SetService method is the only public entry point.
  !! SetServices registers all of the user-provided subroutines
  !! in the module with the NUOPC layer.
  !!
  subroutine SetServices(model, rc)
  
    implicit none
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    character(len=*),parameter   :: subname='(fvcom_cap:SetServices)'
    ! Local variables
    integer                      :: num,i
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    !Assume no need to change clock settings
    ! attach specializing method(s)
!    call NUOPC_CompSpecialize(model, specLabel=model_label_SetClock, &
!      specRoutine=SetClock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    !comment out for now to avoid over writing NOUPC check import
!    call NUOPC_CompSpecialize(model, specLabel=model_label_CheckImport, &
!      specPhaseLabel="RunPhase1", specRoutine=CheckImport, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
      specRoutine=FVCOM_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!!!JQI    call FVCOM_FieldsSetup()
!
!    do num = 1,fldsToFvcom_num
!        print *,  "fldsToFvcom_num  ", fldsToFvcom_num
!        print *,  "fldsToFvcom(num)%shortname  ", fldsToFvcom(num)%shortname
!        print *,  "fldsToFvcom(num)%stdname  ", fldsToFvcom(num)%stdname
!     write(info,*) subname,'fldsToFvcom(num)%stdname  ', fldsToFvcom(num)%stdname
!!     write(info,*) subname,"fldsToFvcom(num)%shortname  ", fldsToFvcom(num)%shortname
!   end do
!
!   do num = 1,fldsFrFvcom_num
!     !print *,  "fldsFrFvcom_num  ", fldsFrFvcom_num
!     !print *,  "fldsFrFvcom(num)%shortname  ", fldsFrFvcom(num)%shortname
!     !print *,  "fldsFrFvcom(num)%stdname  ", fldsFrFvcom(num)%stdname
!     write(info,*) subname,'fldsFrFvcom(num)%stdname  ', fldsFrFvcom(num)%stdname
!   end do
!
    
    write(info,*) subname,' --- fvcom SetServices completed --- '
    !print *,      subname,' --- fvcom SetServices completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
  end subroutine SetServices
  
  !-----------------------------------------------------------------------------
  !> First initialize subroutine called by NUOPC.  The purpose
  !! is to set which version of the Initialize Phase Definition (IPD)
  !! to use.
  !!
  !! For this FVCOM cap, we are using IPDv01.
  !!
  !! @param model an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code

  subroutine InitializeP1(model, importState, exportState, clock, rc)
    
    implicit none
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    type(ESMF_VM)                :: vm
    integer                      :: esmf_comm,fvcom_comm,ierr
    character(len=*),parameter   :: subname='(fvcom_cap:InitializeP1)'

    rc = ESMF_SUCCESS

    write(info,*) subname,' --- initialization phase 1 begin --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    ! details Get current ESMF VM.
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! details Get MPI_communicator from ESMF VM.
    call ESMF_VMGet(vm, mpiCommunicator=esmf_comm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    print*, 'MPI_1 ',esmf_comm
    call MPI_Comm_dup(esmf_comm, fvcom_comm, ierr)
!    print*, 'MPI_2 ',esmf_comm,fvcom_comm
    fvcom_comm = esmf_comm
!    print*, 'MPI_3 ',esmf_comm,fvcom_comm
    ! Duplicate the MPI communicator not to interfere with ESMF communications.
    ! The duplicate MPI communicator can be used in any MPI call in the user
    ! code. Here the MPI_Barrier() routine is called.
    call MPI_Barrier(fvcom_comm, ierr)
    !Initialize fvcom before setting up fields
    
    !NUOPC4MET = .true.
    !NUOPC4WAV = .true.
    
    call NUOPC_FVCOM_init(fvcom_comm,fvcom_name)

    call FVCOM_FieldsSetup()

    !WTIMINC = 3*3600.0
    !RSTIMINC = adc_cpl_int +  adc_cpl_num / adc_cpl_den
    !print *,   'WTIMINC > ',WTIMINC,'  RSTIMINC > ',RSTIMINC

    call FVCOM_AdvertiseFields(importState, fldsToFvcom_num, fldsToFvcom, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call FVCOM_AdvertiseFields(exportState, fldsFrFvcom_num, fldsFrFvcom, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) subname,' --- initialization phase 1 completed --- '
    !print *,      subname,' --- initialization phase 1 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine InitializeP1


  !> Advertises a set of fields in an ESMF_State object by calling
  !! NUOPC_Advertise in a loop.
  !!
  !! @param state the ESMF_State object in which to advertise the fields
  !! @param nfields number of fields
  !! @param field_defs an array of fld_list_type listing the fields to advertise
  !! @param rc return code
  subroutine FVCOM_AdvertiseFields(state, nfields, field_defs, rc)
  
    implicit none
    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(fvcom_cap:FVCOM_AdvertiseFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields
      !print *, 'Advertise: '//trim(field_defs(i)%stdname)//'---'//trim(field_defs(i)%shortname)
      call ESMF_LogWrite('Advertise: '//trim(field_defs(i)%stdname), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call NUOPC_Advertise(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    !print *,      subname,' --- IN   --- '
    write(info,*) subname,' --- completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine FVCOM_AdvertiseFields
  !
  !----------------------------------------------------------------------------------
  subroutine FVCOM_FieldsSetup
  
    implicit none
    integer                    :: rc,k
    character*(ESMF_MAXSTR)    :: tmpName, tmpShortName
    character(len=*),parameter :: subname='(fvcom_cap:FVCOM_FieldsSetup)'
    logical :: hasEntry

    !--------- import fields to Sea FVCOM -------------
    !TODO: Consider moving these lines to driver to avoid doing it in both CAPS

    write(tmpName,'(a)') 'significant_wave_height'
    write(tmpShortName,'(a)') 'wavhss'

    hasEntry = NUOPC_FieldDictionaryHasEntry(trim(tmpName),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if(.not. hasEntry) call NUOPC_FieldDictionaryAddEntry(trim(tmpName), "mx", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=trim(tmpName), shortname=trim(tmpShortName))

    write(tmpName,'(a)') 'average_wave_length'
    write(tmpShortName,'(a)') 'wavlen'

    hasEntry = NUOPC_FieldDictionaryHasEntry(trim(tmpName),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if(.not. hasEntry) call NUOPC_FieldDictionaryAddEntry(trim(tmpName), "mx", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=trim(tmpName), shortname=trim(tmpShortName))

    write(tmpName,'(a)') 'average_wave_direction'
    write(tmpShortName,'(a)') 'wavdir'

    hasEntry = NUOPC_FieldDictionaryHasEntry(trim(tmpName),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if(.not. hasEntry) call NUOPC_FieldDictionaryAddEntry(trim(tmpName), "mx", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=trim(tmpName), shortname=trim(tmpShortName))

    !--------- import fields from atm to Fvcom -------------
!    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=
!    "air_pressure_at_sea_level", shortname= "pmsl" )
!    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=
!    "inst_merid_wind_height10m", shortname= "imwh10m" )
!    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom,
!    stdname="inst_zonal_wind_height10m" , shortname= "izwh10m" )

    !--------- export fields from Sea Fvcom -------------
    call fld_list_add(num=fldsFrFvcom_num, fldlist=fldsFrFvcom, &
         stdname="sea_surface_height_above_sea_level", shortname= "seahgt" )
    call fld_list_add(num=fldsFrFvcom_num, fldlist=fldsFrFvcom, &
         stdname="surface_eastward_sea_water_velocity", shortname= "uucurr" )
    call fld_list_add(num=fldsFrFvcom_num, fldlist=fldsFrFvcom, &
         stdname="surface_northward_sea_water_velocity",shortname= "vvcurr" )

!NEMS hycod standard names>>>
!https://esgf.esrl.noaa.gov/projects/couplednems/coupling_fields
!ocn_current_zonal
!ocncz
!m s-1	Ocean current X component.	 	 	 	
!ocn_current_merid
!ocncm
!m s-1	Ocean current Y component.	 	

   !
    write(info,*) subname,' --- Passed--- '
    !print *,      subname,' --- Passed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine FVCOM_FieldsSetup

  !---------------------------------------------------------------------------------

  subroutine fld_list_add(num, fldlist, stdname, data, shortname, unit)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    implicit none
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    real(ESMF_KIND_R8), dimension(:), optional, target :: data
    character(len=*),    intent(in),optional :: shortname
    character(len=*),    intent(in),optional :: unit

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(fvcom_cap:fld_list_add)'

    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=rc)
      return
    endif

    fldlist(num)%stdname = trim(stdname)
    if (present(shortname)) then
      fldlist(num)%shortname = trim(shortname)
    else
      fldlist(num)%shortname = trim(stdname)
    endif

    if (present(data)) then
      fldlist(num)%assoc     = .true.
      fldlist(num)%farrayPtr => data
    else
      fldlist(num)%assoc     = .false.
    endif

    if (present(unit)) then
      fldlist(num)%unit      = unit
    endif


    write(info,*) subname,' --- Passed--- '
    !print *,      subname,' --- Passed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine fld_list_add

  !-----------------------------------------------------------------------------
  !> Called by NUOPC to realize import and export fields.

  !! The fields to import and export are stored in the fldsToFvcom and fldsFrFvcom
  !! arrays, respectively.  Each field entry includes the standard name,
  !! information about whether the field's grid will be provided by the cap,
  !! and optionally a pointer to the field's data array.  Currently, all fields
  !! are defined on the same mesh defined by the cap.
  !! The fields are created by calling fvcom_cap::fvcom_XXXXXXXXXXXXXXXXXXX.
  !!
  !! @param model an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code

  subroutine InitializeP2(model, importState, exportState, clock, rc)
    
    implicit none
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field)        :: field
    !Saeed added
    type(meshdata)               :: mdata
    type(ESMF_Mesh)              :: ModelMesh,meshIn,meshOut
    type(ESMF_VM)                :: vm
    integer                      :: localPet, petCount
    character(len=*),parameter   :: subname='(fvcom_cap:RealizeFieldsProvidingGrid)'

    rc = ESMF_SUCCESS

    !> \details Get current ESMF VM.
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Get query local pet information for handeling global node information
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    ! call ESMF_VMPrint(vm, rc=rc)

    !! Assign VM to mesh data type.
    mdata%vm = vm
    !print *,localPet,"< LOCAL pet, ADC ..1.............................................. >> "
    ! create a Mesh object for Fields
    !call extract_parallel_data_from_mesh(ROOTDIR, mdata, localPet)
    call extract_parallel_data_from_mesh(mdata, localPet)
    !    print *,"ADC ..2.............................................. >> "
    call create_parallel_esmf_mesh_from_meshdata(mdata,ModelMesh )
    !    print *,"ADC ..3.............................................. >> "
    !
    
    call ESMF_MeshWrite(ModelMesh, filename="fvcom_mesh.nc", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !print *,"FVCOM >> "
    !print *,"NumNd", mdata%NumNd
    !print *,"NumOwnedNd", mdata%NumOwnedNd
    !print *,"NumEl", mdata%NumEl
    !print *,"NumND_per_El", mdata%NumND_per_El


    meshIn  = ModelMesh ! for now out same as in
    meshOut = meshIn

    mdataIn  = mdata
    mdataOut = mdata

!    print *,"FVCOM ..4.............................................. >> "
    call FVCOM_RealizeFields(importState, meshIn , mdata, fldsToFvcom_num, fldsToFvcom, "FVCOM import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    !
    call FVCOM_RealizeFields(exportState, meshOut, mdata, fldsFrFvcom_num, fldsFrFvcom, "FVCOM export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) subname,' --- initialization phase 2 completed --- '
!    print *,      subname,' --- initialization phase 2 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
  end subroutine InitializeP2
  

  !> Adds a set of fields to an ESMF_State object.  Each field is wrapped
  !! in an ESMF_Field object.  Memory is either allocated by ESMF or
  !! an existing FVCOM pointer is referenced.
  !!
  !! @param state the ESMF_State object to add fields to
  !! @param grid the ESMF_Grid object on which to define the fields
  !! @param nfields number of fields
  !! @param field_defs array of fld_list_type indicating the fields to add
  !! @param tag used to output to the log
  !! @param rc return code
  subroutine FVCOM_RealizeFields(state, mesh, mdata, nfields, field_defs, tag, rc)

    implicit none
    type(ESMF_State), intent(inout)    :: state
    type(ESMF_Mesh), intent(in)        :: mesh
    type(meshdata)                     :: mdata
    integer, intent(in)                :: nfields
    type(fld_list_type), intent(inout) :: field_defs(:)
    character(len=*), intent(in)       :: tag
    integer, intent(inout)             :: rc


    type(ESMF_Field)                   :: field
    type(ESMF_DistGrid)                :: nodeDistgrid
    type(ESMF_Array)                   :: array
    type(ESMF_RouteHandle)             :: haloHandle
    integer                            :: i
    character(len=*),parameter         :: subname='(fvcom_cap:FVCOM_RealizeFields)'

    rc = ESMF_SUCCESS
    
    ! Get node DistGrid from the Mesh.
    call ESMF_MeshGet(mesh, nodalDistgrid=nodeDistgrid, rc=rc)
    
    ! Create an ESMF Array with a halo region from a node DistGrid.

    do i = 1, nfields

      array = ESMF_ArrayCreate(nodeDistgrid, typekind=ESMF_TYPEKIND_R8,  &
              haloSeqIndexList=mdata%owned_to_present_halo_nodes, name=field_defs(i)%shortname,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      field = ESMF_FieldCreate(name=field_defs(i)%shortname, mesh=mesh, array=array, &
              meshLoc=ESMF_MESHLOC_NODE, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
	  
      ! Create the RouteHandle for the halo communication.
      call ESMF_FieldHaloStore(field,routehandle=haloHandle, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
	
      ! Do the halo communication.
!!JQI      call ESMF_FieldHalo(field, routehandle=haloHandle, rc=rc)  

        if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

          call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
            ESMF_LOGMSG_INFO, &
            line=__LINE__, &
            file=__FILE__, &
            rc=dbrc)

          !print *,      subname,' --- Connected --- '
          field_defs(i)%connected = .true.
        else
          call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
            ESMF_LOGMSG_INFO, &
            line=__LINE__, &
            file=__FILE__, &
            rc=dbrc)
          ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
          !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
          ! remove a not connected Field from State
          call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          !print *,      subname," Field ", field_defs(i)%stdname ,' --- Not-Connected --- '
          field_defs(i)%connected = .false.
        end if
    end do

    ! After its last use the RouteHandle can be released.
!    call ESMF_FieldHaloRelease(haloHandle, rc=rc)
	
    ! The Field can now be destroyed.
!    call ESMF_FieldDestroy(field, rc=rc)
	
    ! The Array can now be destroyed.
!    call ESMF_ArrayDestroy(array, rc=rc)
	
    write(info,*) subname,' --- OUT--- '
    !print *,      subname,' --- OUT --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
  end subroutine FVCOM_RealizeFields
  !-----------------------------------------------------------------------------

!  subroutine SetClock_mine_not_active(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc

!    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: ADCTimeStep

!    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
!    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    ! initialize internal clock
    ! - on entry, the component clock is a copy of the parent clock
    ! - the parent clock is on the slow timescale atm timesteps
    ! - reset the component clock to have a timeStep that is for fvcom-wav of the parent
    !   -> timesteps
    
    !call ESMF_TimeIntervalSet(ADCTimeStep, s=     adc_cpl_int, sN=adc_cpl_num, sD=adc_cpl_den, rc=rc) ! 5 minute steps
    !TODO: use nint !!?
!    call ESMF_TimeIntervalSet(ADCTimeStep, s= adc_cpl_int , rc=rc) ! 5 minute steps
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    call NUOPC_CompSetClock(model, clock, ADCTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
      
!    print *, "ADC Timeinterval1 = "
!    call ESMF_TimeIntervalPrint(ADCTimeStep, options="string", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!        return  ! bail out    

!  end subroutine
  !-----------------------------------------------------------------------------
  ! From CICE model uses same clock as parent gridComp
!  subroutine SetClock_not_active(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc
!    
!    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: ADCTimeStep, timestep
!    character(len=*),parameter  :: subname='(fvcom_cap:SetClock)'

!    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
!    call ESMF_GridCompGet(model, clock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    !call ESMF_TimeIntervalSet(ADCTimeStep, s=     adc_cpl_int, sN=adc_cpl_num, sD=adc_cpl_den, rc=rc) ! 5 minute steps
    ! tcraig: dt is the cice thermodynamic timestep in seconds
!   call ESMF_TimeIntervalSet(timestep, s=adc_cpl_int, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    call ESMF_ClockSet(clock, timestep=timestep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
      
    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
!    call ESMF_TimeIntervalSet(ADCTimeStep, s=adc_cpl_int, rc=rc) 
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!    call NUOPC_CompSetClock(model, clock, ADCTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
!  end subroutine


  !-----------------------------------------------------------------------------

  !> Called by NUOPC to advance the FVCOM model a single timestep >>>>>
  !! <<<<<<<<<  TODO: check! this is not what we want!!!.
  !!
  !! This subroutine copies field data out of the cap import state and into the
  !! model internal arrays.  Then it calls FVCOM_Run to make a NN timesteps.
  !! Finally, it copies the updated arrays into the cap export state.
  !!
  !! @param model an ESMF_GridComp object
  !! @param rc return code
  !-----------------------------------------------------------------------------
  subroutine ModelAdvance(model, rc)
  
    implicit none
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_State)           :: importState, exportState
    type(ESMF_Time)            :: currTime
    type(ESMF_TimeInterval)    :: timeStep
    character(len=*),parameter :: subname='(fvcom_cap:ModelAdvance)'
    !tmp vector
!    real(ESMF_KIND_R8), pointer:: tmp(:)
    real(ESMF_KIND_R8), ALLOCATABLE, TARGET :: tmp(:)

    ! exports
    real(ESMF_KIND_R8), pointer:: dataPtr_zeta(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_velx(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_vely(:)

    !imports
    real(ESMF_KIND_R8), pointer:: dataPtr_wavhs(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_wavlen(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_wavdir(:)

!    real(ESMF_KIND_R8), pointer:: dataPtr_pmsl(:)
!    real(ESMF_KIND_R8), pointer:: dataPtr_imwh10m(:)
!    real(ESMF_KIND_R8), pointer:: dataPtr_izwh10m(:)

    type(ESMF_StateItem_Flag)  :: itemType
    type(ESMF_Mesh)            :: mesh
    type(ESMF_Field)           :: lfield
    character(len=128)         :: fldname,timeStr
    integer                    :: i1,num
    integer                    :: ITIME_BGN_FVCOM, ITIME_END_FVCOM
    integer                    :: nCplFVCOM
    real(ESMF_KIND_R8)         :: timeStepAbs
    ! local variables for Get methods
    integer :: YY, MM, DD, H, M, S
    integer :: ss,ssN,ssD
    integer :: j, k,ierr
    character*(ESMF_MAXSTR)    :: tmpShortName
    logical :: wave_forcing, meteo_forcing, surge_forcing

    type(ESMF_Time) :: BeforeCaribbeanTime,AfterCaribbeanTime
    
    real(ESMF_KIND_R8), allocatable :: unode(:),vnode(:)
    
    real(sp), allocatable :: TMP_WHS(:),TMP_WLEN(:),TMP_WDIR(:)
    real(sp), allocatable :: TMP1_WHS(:),TMP1_WLEN(:),TMP1_WDIR(:)
    real(sp), allocatable :: WHS(:),WLEN(:),WDIR(:)

    rc = ESMF_SUCCESS
    dbrc = ESMF_SUCCESS
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.

!    if(msr)then
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing FVCOM from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!    end if

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    if(msr)then
    call ESMF_TimePrint(currTime + timeStep, &
      preString="------------FVCOM---------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!    end if

    call ESMF_TimeGet(currTime, yy=YY, mm=MM, dd=DD, h=H, m=M, s=S, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info, *)  "FVCOM currTime = ", YY, "/", MM, "/", DD," ", H, ":", M, ":", S
!JQI    call allMessage(1,info)

    call ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !global values for fvcom time step number
!JQI    ITIME_BGN_FVCOM = ITHS + 1   !if hot start then ITHS>0
!JQI    ITIME_END_FVCOM = NT         !NT is set in read_input.F

    call ESMF_TimeIntervalGet(timeStep, s=ss,sN=ssN,sD=ssD,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    timeStepAbs = real(ss) + real(ssN/ssD)
!JQI    if (mod(timeStepAbs,DTDP) .EQ. 0) then
!JQI      nCplFVCOM = nint(timeStepAbs/DTDP)
!JQI    else
!JQI      nCplFVCOM = nint(timeStepAbs/DTDP)+1
!JQI    endif

    if (mod(timeStepAbs,IntStep_Seconds) .EQ. 0) then
      nCplFVCOM = nint(timeStepAbs/IntStep_Seconds)
    else
      STOP
!JQI      nCplFVCOM = nint(timeStepAbs/DTDP)+1
    endif


     !print *, '  nCplFVCOM = ', nCplFVCOM
     !print *, '  FVCOMCouplingTimeInterval = ', timeStepAbs
!    print*, "STOP HERE?????????"
    !-----------------------------------------
    !   IMPORT
    !-----------------------------------------
    !Get and fill imported fields
    wave_forcing= .true.
   
    do num = 1,fldsToFvcom_num

      if (fldsToFvcom(num)%shortname == 'wavhss') &
         wave_forcing = wave_forcing .and. fldsToFvcom(num)%connected
      if (fldsToFvcom(num)%shortname == 'wavlen') &
         wave_forcing = wave_forcing .and. fldsToFvcom(num)%connected
      if (fldsToFvcom(num)%shortname == 'wavdir') &
         wave_forcing = wave_forcing .and. fldsToFvcom(num)%connected
!      print*,"WAVE PARAMETER IS ",num, wave_forcing, fldsToFvcom(num)%connected
    end do

!    print*,"WAVE FORCING IS ",wave_forcing, fldsToFvcom(1)%connected
!    print*,"SHORT NAME ",fldsToFvcom(1)%shortname,fldsToFvcom(2)%shortname,fldsToFvcom(3)%shortname
    if (wave_forcing) then
    
      ! Wave time step
!JQI        RSTIMINC = nint(timeStepAbs)  !TODO: Get it from coupler based on the time slots.
!JQI                                      !TODO: This implemetation wotks for one time slot for wind and wave right now.
!JQI                                      !TODO: 
!JQI        !RSTIMINC = fvcom_cpl_int +  fvcom_cpl_num / fvcom_cpl_den
!JQI        !print *, ' in cap   ....> RSTIMINC > ', RSTIMINC

      !-----------------------------------------
      ! <<<<< RECEIVE and UN-PACK WAVHSS

      call State_getFldPtr_(ST=importState,fldname='wavhss',fldptr=dataPtr_wavhs, &
                            dump=.true.,timeStr=timeStr,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
!      dataPtr_sxx3(:,k) =  dataPtr_tmp(:)
!      print*,'3D size = ', size(dataPtr_tmp),size(FVCOM_SXX3,1),M,MT
        
      !-----------------------------------------
      ! <<<<< RECEIVE and UN-PACK WAVLEN
      call State_getFldPtr(ST=importState,fldname='wavlen',fldptr=dataPtr_wavlen,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !-----------------------------------------
      ! <<<<< RECEIVE and UN-PACK WAVDIR
      call State_getFldPtr(ST=importState,fldname='wavdir',fldptr=dataPtr_wavdir,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

!JQI        write(800+myid,*) size(dataPtr_sxx2)
!JQI        do i1=1,mdataOut%NumOwnedNd_NOHALO
 !JQI         write(800+myid,*) i1, mdataOut%owned_to_present_nodes(i1),dataPtr_sxx2(i1)
!JQI	end do  

      !print *, 'size dataPtr_sxy > ', size(dataPtr_sxy)
      !print *, 'in cap maxval(RSNX2)', maxval(RSNX2)
      !print *, 'in cap maxval(RSNY2)', maxval(RSNY2)
      ! Allocate arrays for radiation stresses.

      IF(.NOT.ALLOCATED(FVCOM_WHS))  ALLOCATE(FVCOM_WHS(0:MT,1:2))
      IF(.NOT.ALLOCATED(FVCOM_WLEN)) ALLOCATE(FVCOM_WLEN(0:MT,1:2))
      IF(.NOT.ALLOCATED(FVCOM_WDIR)) ALLOCATE(FVCOM_WDIR(0:MT,1:2))
      FVCOM_WHS  = 0.0 
      FVCOM_WLEN = 0.0
      FVCOM_WDIR = 0.0

      ! update last rad str
      FVCOM_WHS(:,1)  = FVCOM_WHS(:,2)
      FVCOM_WLEN(:,1) = FVCOM_WLEN(:,2)
      FVCOM_WDIR(:,1) = FVCOM_WDIR(:,2)

      ! Fill owned nodes from imported data to model variable
      ! devide by water density to convert from N.m-2 to m2s-2
	
!      print*, "SIZE =",size(dataPtr_sxx2),mdataIn%NumOwnedNd_NoHalo,  &
!                                          mdataIn%NumOwnedNd_Halo
	
      do i1 = 1, mdataIn%NumOwnedNd_NoHalo
        FVCOM_WHS(mdataIn%owned_to_present_nodes(i1),1)  = dataPtr_wavhs(i1)
        FVCOM_WLEN(mdataIn%owned_to_present_nodes(i1),1) = dataPtr_wavlen(i1)
        FVCOM_WDIR(mdataIn%owned_to_present_nodes(i1),1) = dataPtr_wavdir(i1)
      end do
      do i1 = 1, mdataIn%NumOwnedNd_Halo
        FVCOM_WHS(mdataIn%owned_to_present_halo_nodes(i1),1)  = dataPtr_wavhs(i1+mdataIn%NumOwnedNd_NoHalo)
        FVCOM_WLEN(mdataIn%owned_to_present_halo_nodes(i1),1) = dataPtr_wavlen(i1+mdataIn%NumOwnedNd_NoHalo)
        FVCOM_WDIR(mdataIn%owned_to_present_halo_nodes(i1),1) = dataPtr_wavdir(i1+mdataIn%NumOwnedNd_NoHalo)
      end do

      IF(PAR)THEN
!        CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,    &
!                   FVCOM_SXX2(:,1),FVCOM_SYY2(:,1))
       
        IF(.NOT. ALLOCATED(TMP_WHS))  ALLOCATE(TMP_WHS(0:MGL))
        IF(.NOT. ALLOCATED(TMP_WLEN)) ALLOCATE(TMP_WLEN(0:MGL))
        IF(.NOT. ALLOCATED(TMP_WDIR)) ALLOCATE(TMP_WDIR(0:MGL))
        IF(.NOT. ALLOCATED(TMP1_WHS)) ALLOCATE(TMP1_WHS(0:MT))
        IF(.NOT. ALLOCATED(TMP1_WLEN)) ALLOCATE(TMP1_WLEN(0:MT))
        IF(.NOT. ALLOCATED(TMP1_WDIR)) ALLOCATE(TMP1_WDIR(0:MT))

        TMP1_WHS(:)  = FVCOM_WHS(:,1)
        TMP1_WLEN(:) = FVCOM_WLEN(:,1)
        TMP1_WDIR(:) = FVCOM_WDIR(:,1)

        CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,TMP1_WHS,TMP_WHS)
        CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,TMP1_WLEN,TMP_WLEN)
        CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,TMP1_WDIR,TMP_WDIR)

        CALL MPI_BCAST(TMP_WHS,MGL,MPI_F,0,MPI_FVCOM_GROUP,IERR)
        CALL MPI_BCAST(TMP_WLEN,MGL,MPI_F,0,MPI_FVCOM_GROUP,IERR)
        CALL MPI_BCAST(TMP_WDIR,MGL,MPI_F,0,MPI_FVCOM_GROUP,IERR)

        do i1 = 1, mt
        FVCOM_WHS(i1,1)    = TMP_WHS(NGID_X(i1))
        FVCOM_WLEN(i1,1)   = TMP_WLEN(NGID_X(i1))
        FVCOM_WDIR(i1,1)   = TMP_WDIR(NGID_X(i1))
        end do

       ! DEALLOCATE(TMP_SXX2,TMP_SYY2)
       ! DEALLOCATE(TMP1_SXX2,TMP1_SYY2)
      END IF
      FVCOM_WHS(:,2)  = FVCOM_WHS(:,1)
      FVCOM_WLEN(:,2) = FVCOM_WLEN(:,1)
      FVCOM_WDIR(:,2) = FVCOM_WDIR(:,1)

      !print *, 'Hard Coded >>>>>  SXX >>>>>> '
      !mask
      where(abs(FVCOM_WHS) .gt. 1e6)  FVCOM_WHS  =  1e-10
      where(abs(FVCOM_WLEN).gt. 1e6)  FVCOM_WLEN =  1e-10
      where(abs(FVCOM_WDIR).gt. 1e6)  FVCOM_WDIR =  1e-10

      !max values
      where(FVCOM_WHS  >  1e3)  FVCOM_WHS  =  1e3
      where(FVCOM_WLEN >  1e3)  FVCOM_WLEN =  1e3
      where(FVCOM_WDIR >  1e3)  FVCOM_WDIR =  1e3
      where(FVCOM_WHS  < -1e3)  FVCOM_WHS  = -1e3
      where(FVCOM_WLEN < -1e3)  FVCOM_WLEN = -1e3
      where(FVCOM_WDIR < -1e3)  FVCOM_WDIR = -1e3

      if(.not.allocated(whs))  allocate(whs(0:MT))
      if(.not.allocated(wlen)) allocate(wlen(0:MT))
      if(.not.allocated(wdir)) allocate(wdir(0:MT))

      whs(0:MT)  = FVCOM_WHS(0:MT,1)
      wlen(0:MT) = FVCOM_WLEN(0:MT,1)
      wdir(0:MT) = FVCOM_WDIR(0:MT,1)

      call radiation_stress(whs,wlen,wdir)

        CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,WHS,TMP_WHS)
        CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,WLEN,TMP_WLEN)
        CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,WDIR,TMP_WDIR)

      do i1=1,mgl
        write(100+myid,*) i1,tmp_whs(i1)
      end do
      do i1=1,mgl
        write(200+myid,*) i1,tmp_wlen(i1)
      end do
      do i1=1,mgl
        write(300+myid,*) i1,tmp_wdir(i1)
      end do


!JQI!    NOTE: FVCOM accepts wave-driven stresses "in units of velocity squared
!JQI!    (consistent with the units of gravity).  Stress in these units is obtained
!JQI!    by dividing stress in units of force/area by the reference density of water."


!
!
!From WW3 source : w3iogomd.ftn  line: 1756
!constants.ftn:34:!      DWAT      Real  Global   Density of water               (kg/m3)
!constants.ftn:63:      REAL, PARAMETER         :: DWAT   = 1000.
!
!w3iogomd.ftn:1753:      SXX    = SXX * DWAT * GRAV
!w3iogomd.ftn:1754:      SYY    = SYY * DWAT * GRAV
!w3iogomd.ftn:1755:      SXY    = SXY * DWAT * GRAV
!
!    Rad.Str info from netcdf header
!		 sxx:long_name = "Radiation stress component Sxx" ;
!		 sxx:standard_name = "radiation_stress_component_sxx" ;
!		 sxx:globwave_name = "significant_wave_height" ;
!		 sxx:units = "N m-1" ;
!		 sxx:_FillValue = 9.96921e+36f ;
!		 sxx:scale_factor = 1.f ;
!		 sxx:add_offset = 0.f ;
!		 sxx:valid_min = -3000 ;
!	 	 sxx:valid_max = 3000 ;
 
!    Therefore we need to divide sxx/rho to change its unit to m3s-2
!    in force calculation we do d(sxx)/dx therefore the final force 
!    unit will be m2s-2 which is the correct one.   

      !print *, 'in cap maxval(FVCOM_SXX)', maxval(FVCOM_SXX)

!JQI        InterpoWeight = 0.0  !avoid time interpolation

      ! Calculate wave forces
!JQI        call ComputeWaveDrivenForces

      !iunit_log = iunit_log + 1
      !open(unit = iunit_log, ACTION = "write", STATUS ="replace" )
      !write(iunit_log, *) RSNX2, RSNY2
      !close(iunit_log)

!JQI        write(info,*) subname,' --- wave data exchange OK / wave feilds are all connected --- / Model advances '
!JQI        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
      !print *, info

      !RSNX2 = 0.0001
      !RSNY2 = 0.0001
    
      !RSNX2 = 0.0
      !RSNY2 = 0.0

      ! initialize time reach to Caribean Islands
      !call ESMF_TimeSet(BeforeCaribbeanTime, yy=2008, mm=9, dd=6 , h=12, m=0, s=0, rc=rc)
      !call ESMF_TimeSet(AfterCaribbeanTime , yy=2008, mm=9, dd=12, h=12, m=0, s=0, rc=rc)

!!!!#ifdef NO_COMPILE00000
        !-----------------------
   !     if ((currTime > BeforeCaribbeanTime) .and. (currTime < AfterCaribbeanTime)) then
!JQI            write(info,*) subname, 'in cap after maxval(RSNX2)', maxval(RSNX2)
!JQI            call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

!JQI            write(info,*) subname, 'in cap after maxval(RSNY2)', maxval(RSNY2)
!JQI            call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

            !print *, 'Hard Coded >>>>>>>>>>>  where(abs(RSNX2).gt. wave_force_limmit) RSNX2 =  wave_force_limmit'
            !write(info,*) subname,'Hard Coded >>>>>>>>>>>  where(abs(RSNX2).gt. wave_force_limmit) RSNX2 =  wave_force_limmit'
            !call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

            !where(RSNX2.gt. wave_force_limmit) RSNX2 =  wave_force_limmit
            !where(RSNY2.gt. wave_force_limmit) RSNY2 =  wave_force_limmit

            !where(RSNX2.le. (-1.0 * wave_force_limmit)) RSNX2 =  -1.0 * wave_force_limmit
            !where(RSNY2.le. (-1.0 * wave_force_limmit)) RSNY2 =  -1.0 * wave_force_limmit

!!!!#endif

    !        endif

    else
      NUOPC4WAV = .false.
      write(info,*) subname,' --- no wave forcing exchange / waves are not all connected --- '
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
      !print *, info
      !stop
    end if        
    !-----------------------------------------
    !   IMPORT from ATM
    !-----------------------------------------
    meteo_forcing= .false.
!JQI   do num = 1,fldsToFvcom_num
!JQI      if (fldsToFvcom(num)%shortname == 'pmsl')    meteo_forcing = meteo_forcing .and. fldsToFvcom(num)%connected
!JQI      if (fldsToFvcom(num)%shortname == 'imwh10m') meteo_forcing = meteo_forcing .and. fldsToFvcom(num)%connected
!JQI      if (fldsToFvcom(num)%shortname == 'izwh10m') meteo_forcing = meteo_forcing .and. fldsToFvcom(num)%connected
!JQI   end do
    
    if ( meteo_forcing) then
!JQI        !NWS = 39   ! over write NWS option to be sure we incldue wind forcing
        
!JQI        !WTIMINC = wrf_int 
!JQI        !WTIMINC = wrf_int +  wrf_num / wrf_den
!JQI        WTIMINC =  nint(timeStepAbs)  !TODO: Get it from coupler based on the time slots.
                                      !TODO: This implemetation wotks for one time slot for wind and wave right now.
                                      !TODO: 
        !Get and fill imported fields
        ! <<<<< RECEIVE and UN-PACK pmsl
!JQI        call State_getFldPtr_(ST=importState,fldname='pmsl',fldptr=dataPtr_pmsl,rc=rc,dump=.true.,timeStr=timeStr)
!JQI        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!JQI          line=__LINE__, &
!JQI          file=__FILE__)) &
!JQI          return  ! bail out

        !print *, 'size > 5 dataPtr_pmsl>', size(dataPtr_pmsl) !, dataPtr_pmsl(5:)
        
        !-----------------------------------------
        ! <<<<< RECEIVE and UN-PACK imwh10m    V-Y wind comp
!JQI        call State_getFldPtr_(ST=importState,fldname='imwh10m',fldptr=dataPtr_imwh10m,rc=rc,dump=.false.,timeStr=timeStr)
!JQI        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!JQI          line=__LINE__, &
!JQI          file=__FILE__)) &
!JQI          return  ! bail out
        !-----------------------------------------
        ! <<<<< RECEIVE and UN-PACK izwh10m    U-X  wind comp
!JQI        call State_getFldPtr_(ST=importState,fldname='izwh10m',fldptr=dataPtr_izwh10m,rc=rc,dump=.false.,timeStr=timeStr)
 !JQI       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!JQI          line=__LINE__, &
!JQI          file=__FILE__)) &
!JQI          return  ! bail out
        
!JQI        write(info,*) subname,' --- meteo forcing exchange OK / atm feilds are all connected --- / Model advances '
!JQI        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        !print *, info

        ! Allocate arrays for radiation stresses.
 !JQI       IF(.NOT.ALLOCATED(WVNX1)) ALLOCATE(WVNX1(1:NP))
!JQI        IF(.NOT.ALLOCATED(WVNY1)) ALLOCATE(WVNY1(1:NP))
!JQI        IF(.NOT.ALLOCATED(PRN1) ) ALLOCATE(PRN1 (1:NP))

 !JQI       IF(.NOT.ALLOCATED(WVNX2)) ALLOCATE(WVNX2(1:NP))
!JQI        IF(.NOT.ALLOCATED(WVNY2)) ALLOCATE(WVNY2(1:NP))
!JQI        IF(.NOT.ALLOCATED(PRN2) ) ALLOCATE(PRN2 (1:NP))

        !print *, 'maxval(WVNX2)', maxval(WVNX2)

 !JQI       WVNX1 = WVNX2   
!JQI        WVNY1 = WVNY2  
!JQI        PRN1  = PRN2

        !call UPDATER( dataPtr_izwh10m(:), dataPtr_imwh10m(:), dataPtr_pmsl(:),3)
       
        ! Fill owned nodes from imported data to model variable
        !TODO: unit check
!JQI        do i1 = 1, mdataOut%NumOwnedNd, 1
!JQI            WVNX2(mdataOut%owned_to_present_nodes(i1)) =  dataPtr_izwh10m(i1) !* 0.0  !zonal is u-comp or x-comp
!JQI        end do

!JQI        do i1 = 1, mdataOut%NumOwnedNd, 1
!JQI            WVNY2(mdataOut%owned_to_present_nodes(i1)) =  dataPtr_imwh10m(i1) !* 0.0  !Meridionalis v-comp or y-comp
!JQI        end do
        
!JQI        do i1 = 1, mdataOut%NumOwnedNd, 1
!JQI            PRN2(mdataOut%owned_to_present_nodes(i1) ) = dataPtr_pmsl(i1) / (1025 * 9.81)    !convert Pascal to mH2O
        
          !if ( abs(dataPtr_pmsl(i1) ).gt. 1e11)  then
          !  STOP '  dataPtr_pmsl > mask '     
          !end if
!JQI        end do
         
        ! Ghost nodes update 
!JQI        call UPDATER( WVNX1(:), WVNY1(:), PRN1(:),3)
!JQI        call UPDATER( WVNX2(:), WVNY2(:), PRN2(:),3)

!JQI        if (first_exchange .and. sum(PRN1) .le. 1.0) then
!JQI          WVNX2 = 1e-10
!JQI          WVNY2 = 1e-10
!JQI          WVNX1 = 1e-10
!JQI          WVNY1 = 1e-10

!JQI          PRN2 = 10.0  !hard coded to handel the 1st exchange zeros problem :TODO! Need to resolve this!
!JQI          PRN1 = PRN2
!JQI          first_exchange = .false.
!JQI        end if  
    
        !if (sum(PRN1) .eq. 0.0 ) then
        !  PRN1 = 10000.0
        !end if  
    
        !if (sum(WVNX2) .eq. 0.0 ) then
        !  WVNX2 = 8.0
        !end if  
    
        !if (sum(WVNX1) .eq. 0.0 ) then
        !  WVNX1 = 8.0
        !end if  
        

        !if (sum(WVNY2) .eq. 0.0 ) then
        !  WVNY2 = -8.0
        !end if  
    
        !if (sum(WVNY1) .eq. 0.0 ) then
        !  WVNY1 = -8.0
        !end if  


        !where(abs(PRN1).gt. 1e11)  PRN1 =  1e4
        !where(abs(PRN2).gt. 1e11)  PRN2 =  1e4

        
        ! where(abs(WVNX1).gt. 1e6)  WVNX1 =  8.0
        !where(abs(WVNX2).gt. 1e6)  WVNX2 =  8.0

        !where(abs(WVNY1).gt. 1e6)  WVNY1 =  -8.0
        !where(abs(WVNY2).gt. 1e6)  WVNY2 =  -8.0

            
          !PRN2 = 10000.0
          !PRN1 = 10000.0
          !WVNX2 =  8.0
          !WVNX1 =  8.0
          !WVNY2 = -8.0
          !WVNY1 = -8.0       
            
    
        !where(dataPtr_pmsl .gt. 1e20)  dataPtr_pmsl =  10e4
        !where(dataPtr_pmsl .lt. 8e4 )  dataPtr_pmsl =  10e4        
        
        !print *, 'size(dataPtr_pmsl) > ',  size(dataPtr_pmsl)
        !print *, 'size(PRN2)         > ',  size(PRN2)


        !
        !where((WVNX2).gt. 20)  WVNX2 =  20
        !where((WVNY2).gt. 20)  WVNY2 =  20
        !where((PRN2) .gt. 1e20)   PRN2 =  1e4                

        !where((WVNX1).gt. 20)  WVNX1 =  20
        !where((WVNY1).gt. 20)  WVNY1 =  20
        !where((PRN1).gt.  1e20)   PRN1 =  1e4  
        
        !PRN1 =  1e4  
        !PRN2 =  1e4  
        
        ! Ghost nodes update for meteo infos
        !call UPDATER( WVNX2(:), WVNY2(:), PRN2(:),3)

        !print *, 'in cap before maxval(WVNX2)', maxval(WVNX2)
        !print *, 'in cap before maxval(WVNY2)', maxval(WVNY2)
         !WVNX2 = 10.0
         !WVNY2 = 10.0
         !PRN2  = 10000.0

         !WVNX2 = 0.0
         !WVNY2 = 0.0
         !PRN2  = 10000.0

        !print *, 'in cap after maxval(WVNX2)', maxval(WVNX2)
        !print *, 'in cap after maxval(WVNY2)', maxval(WVNY2)    


!JQI    else
!JQI        NUOPC4MET = .false.
!JQI        write(info,*) subname,' --- no meteo forcing exchange / atm feilds are not all connected --- '
!JQI        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        !print *, info
        !stop
    endif

    surge_forcing= .true.
    do num = 1,fldsFrFvcom_num
      if (fldsFrFvcom(num)%shortname == 'seahgt') &
          surge_forcing = surge_forcing .and. fldsFrFvcom(num)%connected
      if (fldsFrFvcom(num)%shortname == 'uucurr') &
          surge_forcing = surge_forcing .and. fldsFrFvcom(num)%connected
      if (fldsFrFvcom(num)%shortname == 'vvcurr') &
          surge_forcing = surge_forcing .and. fldsFrFvcom(num)%connected
    end do

    !------------------------------------------
    !---------------  RUN  --------------------
    !------------------------------------------
    write(info,*) subname,' --- run phase 2 called --- '
    !print *,      subname,' --- run phase 2 called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_FVCOM_Run(nCplFVCOM)

    write(info,*) subname,' --- run phase 3 called --- '         ! nCplFVCOM = ',nCplFVCOM
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
    !-------------------------------------------

!    surge_forcing = .true.

    if (surge_forcing) then
    
      allocate(unode(mdataOut%NumOwnedNd)); unode = 0.0
      allocate(vnode(mdataOut%NumOwnedNd)); vnode = 0.0
      
      do i1 = 1,mdataOut%NumOwnedNd
        unode(i1) = 0.0
        vnode(i1) = 0.0
	do j = 1, ntve(i1)
	  unode(i1) = unode(i1)+u(nbve(i1,j),1)
	  vnode(i1) = vnode(i1)+v(nbve(i1,j),1)
	end do
	unode(i1) = unode(i1)/ntve(i1)  
	vnode(i1) = vnode(i1)/ntve(i1)  
      end do
      
      !-----------------------------------------
      !   EXPORT
      !-----------------------------------------
      allocate (dataPtr_zeta(mdataOut%NumOwnedNd))

      !pack and send exported fields
      
      allocate (tmp(mdataOut%NumOwnedNd))

      ! >>>>> PACK and send ZETA
      call State_getFldPtr_(ST=exportState,fldname='seahgt',fldptr=dataPtr_zeta, &
!JQI        rc=rc,dump=.true.,timeStr=timeStr)
        rc=rc,dump=.false.,timeStr=timeStr)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

!      print*,"INSIDE FVCOM CAP ",size(dataPtr_zeta)
      !fill only owned nodes for tmp vector
      do i1 = 1, mdataOut%NumOwnedNd_NoHalo
        tmp(i1) = EL(mdataOut%owned_to_present_nodes(i1))
      end do
      do i1 = 1, mdataOut%NumOwnedNd_Halo
        tmp(i1+mdataOut%NumOwnedNd_NoHalo) = EL(mdataOut%owned_to_present_halo_nodes(i1))
      end do

      !assign to field
      dataPtr_zeta = tmp
      
      deallocate(tmp)
      
      allocate(tmp(mdataOut%NumOwnedNd))
      allocate(dataPtr_velx(mdataOut%NumOwnedNd))
      !----------------------------------------
      ! >>>>> PACK and send VELX
      call State_getFldPtr(ST=exportState,fldname='uucurr',fldptr=dataPtr_velx,rc=rc)
!JQI      call State_getFldPtr_(ST=exportState,fldname='velx',fldptr=dataPtr_velx, &
!JQI        rc=rc,dump=.true.,timeStr=timeStr)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill elements for tmp vector
      do i1 = 1, mdataOut%NumOwnedNd_NoHalo
        tmp(i1) = UNODE(mdataOut%owned_to_present_nodes(i1))
      end do
      do i1 = 1, mdataOut%NumOwnedNd_Halo
        tmp(i1+mdataOut%NumOwnedNd_NoHalo) = UNODE(mdataOut%owned_to_present_halo_nodes(i1))
      end do
      !assign to field
      dataPtr_velx = tmp
      deallocate(tmp)
      
      allocate(tmp(mdataOut%NumOwnedNd))
      allocate(dataPtr_vely(mdataOut%NumOwnedNd))
      !----------------------------------------
      ! >>>>> PACK and send VELY
      call State_getFldPtr(ST=exportState,fldname='vvcurr',fldptr=dataPtr_vely,rc=rc)
!JQI      call State_getFldPtr_(ST=exportState,fldname='vely',fldptr=dataPtr_vely, &
!JQI        rc=rc,dump=.true.,timeStr=timeStr)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill elements for tmp vector
      do i1 = 1, mdataOut%NumOwnedNd_NoHalo
        tmp(i1) = VNODE(mdataOut%owned_to_present_nodes(i1))
      end do
      do i1 = 1, mdataOut%NumOwnedNd_Halo
        tmp(i1+mdataOut%NumOwnedNd_NoHalo) = VNODE(mdataOut%owned_to_present_halo_nodes(i1))
      end do
      !assign to field
      dataPtr_vely = tmp
      deallocate(tmp)

      deallocate(unode,vnode)
    else
      write(info,*) subname,' --- no surge forcing for wave. 1way coupled WW3 -> FVCOM  ---'
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
      !print *, info
    end if

!    deallocate(unode)
!    deallocate(vnode)
      
  end subroutine ModelAdvance

!-----------------------------------------------------------
  !> Retrieve a pointer to a field's data array from inside an ESMF_State object.
  !!
  !! @param ST the ESMF_State object
  !! @param fldname name of the fields
  !! @param fldptr pointer to 1D array
  !! @param rc return code
  subroutine State_GetFldPtr_(ST, fldname, fldptr, rc, dump,timeStr)
    
    implicit none
    type(ESMF_State), intent(in) :: ST
    type(ESMF_RouteHandle)       :: haloHandle
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:)
    integer, intent(out), optional :: rc
    logical, intent(in), optional  :: dump
    character(len=128),intent(inout), optional :: timeStr
    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(fvcom_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (present(rc)) rc = lrc

    call ESMF_FieldHaloStore(lfield, routehandle = haloHandle, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !Halo update
    call ESMF_FieldHalo(lfield, routehandle = haloHandle, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    !!
    !!TODO: this should not be here. It should finalize once
!!JQI    call ESMF_FieldHaloRelease (routehandle = haloHandle, rc=lrc)
!!JQI    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 
    !if (present(rc)) rc = lrc
    !write(info,*) ' --- ATM  halo routehandel in work >>>>>  ---'
    !call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    if (dump) then
      if (.not. present(timeStr)) timeStr="_"
      call ESMF_FieldWrite(lfield, &
        fileName='field_ocn_'//trim(fldname)//trim(timeStr)//'.nc', &
        rc=rc,overwrite=.true.)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if
  end subroutine State_GetFldPtr_

  !> Retrieve a pointer to a field's data array from inside an ESMF_State object.
  !!
  !! @param ST the ESMF_State object
  !! @param fldname name of the fields
  !! @param fldptr pointer to 1D array
  !! @param rc return code
  subroutine State_GetFldPtr(ST, fldname, fldptr, rc)
    
    implicit none
    type(ESMF_State), intent(in) :: ST
    type(ESMF_RouteHandle)       :: haloHandle
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:)
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(FVCOM:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldHaloStore(lfield, routehandle = haloHandle, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !Halo update
    call ESMF_FieldHalo(lfield,routehandle = haloHandle, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    !!
    !!TODO: this should not be here. It should finalize once
    call ESMF_FieldHaloRelease (routehandle = haloHandle, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    if (present(rc)) rc = lrc
  end subroutine State_GetFldPtr

  subroutine CheckImport_not_comp(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc
    
    ! This is the routine that enforces correct time stamps on import Fields
    
    ! local variables
    type(ESMF_Clock)        :: driverClock
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_State)        :: importState
    type(ESMF_Field)        :: field
    logical                 :: atCorrectTime

    rc = ESMF_SUCCESS
    return
    
!    ! query Component for the driverClock
!    call NUOPC_ModelGet(model, driverClock=driverClock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
!    ! get the start time and current time out of the clock
!    call ESMF_ClockGet(driverClock, startTime=startTime, &
!      currTime=currTime, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
   
  end subroutine

  !-----------------------------------------------------------------------------
  


  !-----------------------------------------------------------------------------
  !> Called by NUOPC at the end of the run to clean up.  The cap does
  !! this simply by calling FVCOM_Final.
  !!
  !! @param model the ESMF_GridComp object
  !! @param rc return code
  subroutine FVCOM_model_finalize(model, rc)

    implicit none  
    !input arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    !local variables
    type(ESMF_Clock)     :: clock
    type(ESMF_Time)      :: currTime
    character(len=*),parameter  :: subname='(fvcom_cap:fvcom_model_finalize)'

    rc = ESMF_SUCCESS

    write(info,*) subname,' --- finalize called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    call NUOPC_FVCOM_Final()

    write(info,*) subname,' --- finalize completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine FVCOM_model_finalize

!==============================================================================|
   SUBROUTINE RADIATION_STRESS(WHS,WLEN,WDIR)
!------------------------------------------------------------------------------|
   
   USE MOD_PREC
   USE ALL_VARS, ONLY : ISBCE,DLTXC,DLTYC,IEC,IENODE,NE,RAMP,DEG2RAD,GRAV_N,D1,ZZ,D,DZ,DZ1,H1
   IMPLICIT NONE

   REAL(SP), DIMENSION(:) :: WHS,WLEN,WDIR
   REAL(SP), ALLOCATABLE     :: SXX(:,:),SXY(:,:),SYY(:,:)
   REAL(SP), DIMENSION(0:NT,KB) :: PSXXPX,PSXYPX,PSXYPY,PSYYPY
   REAL(SP), DIMENSION(0:NT,KB) :: PSPXPZ,PSPYPZ

   REAL(SP), DIMENSION(0:MT) :: WAVE_NUMBER,WAVE_NUMBER_X,WAVE_NUMBER_Y,SIN_DIR,COS_DIR
   REAL(SP), DIMENSION(0:MT) :: WAVE_ENERGY,KD,WAVE_C
   REAL(SP), DIMENSION(0:MT) :: O_WAVE_NUMBER
   REAL(SP), DIMENSION(0:MT) :: O_COSH,O_SINH,O_2SINH
   REAL(SP) :: EXFLUX

   INTEGER  :: I,K,IA,IB,J1,J2
   REAL(SP) :: FSS,FCS,FSC,FCC
   REAL(SP) :: CFF1,CFF2,CFF3,CFF4,CFF5,CFF6,FAC2,sum3dsxx
   REAL(SP) :: SXXIJ,SXYIJ,SYYIJ,DIJ

   REAL(SP) :: XTMP,XTMP1
   REAL(SP), ALLOCATABLE :: TPZDIST(:,:)


   real(sp), allocatable :: sxx_gl(:),sxy_gl(:),syy_gl(:)
   real(sp), allocatable :: sxx_tmp(:),sxy_tmp(:),syy_tmp(:)
   real(sp), allocatable :: U_STOKES_3D_TMP(:,:), V_STOKES_3D_TMP(:,:)

!==============================================================================|
   REAL(SP), PARAMETER :: KDMAX = 3.0_SP ! Based on MELLOR(2015). For old code: KDMAX = 5.0_SP
   REAL(SP), PARAMETER :: eps1 = 1E-14_SP
   REAL(SP), PARAMETER :: WAVE_LENGTH_MIN = 0.01_SP
   REAL(SP), PARAMETER :: PI = 3.1415926
!---------------Jianzhong----------------------------
   IF(.NOT.ALLOCATED(SXX)) ALLOCATE(SXX(0:MT,KB))
   IF(.NOT.ALLOCATED(SXY)) ALLOCATE(SXY(0:MT,KB))
   IF(.NOT.ALLOCATED(SYY)) ALLOCATE(SYY(0:MT,KB))
   IF(.NOT.ALLOCATED(U_STOKES_3D_TMP)) ALLOCATE(U_STOKES_3D_TMP(0:MT,KB))
   IF(.NOT.ALLOCATED(V_STOKES_3D_TMP)) ALLOCATE(V_STOKES_3D_TMP(0:MT,KB))
!----------------------------------------------------
   ALLOCATE(TPZDIST(0:NT,KB));       TPZDIST     = 0.0_SP

   WAVE_NUMBER   = 0.0_SP   ;WAVE_NUMBER_X = 0.0_SP   ;WAVE_NUMBER_Y = 0.0_SP
   O_COSH        = 0.0_SP   ;O_SINH        = 0.0_SP   ;O_2SINH       = 0.0_SP
   O_WAVE_NUMBER = 0.0_SP

!
!  Compute wave numbers and wave energy.
!
   DO I=1,MT
    WAVE_NUMBER(I) = 2.0_SP*PI/MAX(WLEN(I),WAVE_LENGTH_MIN)
   END DO 
   O_WAVE_NUMBER = 1.0_SP/WAVE_NUMBER
!JQI   SIN_DIR       = SIN(WDIR*DEG2RAD)
!JQI   COS_DIR       = COS(WDIR*DEG2RAD)
   SIN_DIR       = SIN(WDIR)
   COS_DIR       = COS(WDIR)
   WAVE_NUMBER_X = WAVE_NUMBER*COS_DIR
   WAVE_NUMBER_Y = WAVE_NUMBER*SIN_DIR
   WAVE_ENERGY   = 0.0625_SP*GRAV_N*WHS*WHS
!
!  Compute wave celerity and phase velocity.
!
   DO I=1,MT
!     KD(I) = MIN(WAVE_NUMBER(I)*D(I)+eps1,KDMAX)
!JQI error     KD(I) = WAVE_NUMBER(I)*D1(I)+eps1
     KD(I) = WAVE_NUMBER(I)*D(I)+eps1
   END DO 
   
   WHERE(KD <= KDMAX) 
    WAVE_C = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))

    O_COSH  = 1.0_SP/COSH(KD)
    O_SINH  = 1.0_SP/SINH(KD)
    O_2SINH = 1.0_SP/SINH(2.0_SP*KD)
   ELSEWHERE
    WAVE_C = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))
   END WHERE

#if defined(WAVE_ROLLER)
   OROLLER = 0.0_SP;GAMW = 0.0_SP; ROLLA = 0.0_SP
   DO I=1,MT
     GAMW(I) = MIN(D(I)/(HSC1(I)+eps1),5.0_SP) 
     DO K=1,KBM1
        OROLLER(I)=OROLLER(I)+D(I)*DZ(I,K)*(1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4))
     END DO
     OROLLER(I)=1.0_SP/(OROLLER(I)+eps1)
     ROLLA(I)=0.0424*HSC1(I)*QB1(I)*WLEN(I)
   END DO
#endif
   
!----------INITIALIZE STRESS ARRAY ----------------------------------------------!

   SXX    = 0.0_SP   ;SXY    = 0.0_SP   ;SYY    = 0.0_SP
!   SXXA   = 0.0_SP   ;SXYA   = 0.0_SP   ;SYYA   = 0.0_SP
   PSXXPX = 0.0_SP   ;PSXYPX = 0.0_SP   ;PSXYPY = 0.0_SP   ;PSYYPY = 0.0_SP
   PSPXPZ = 0.0_SP   ;PSPYPZ = 0.0_SP

   DO I=1,M
     sum3dsxx=0
    IF(KD(I) <= KDMAX)THEN 
     DO K=1,KBM1
       FAC2 = 1.0_SP+ZZ(I,K)
       FCC  = COSH(KD(I)*FAC2)*O_COSH(I)
       FCS  = COSH(KD(I)*FAC2)*O_SINH(I)
       FSC  = SINH(KD(I)*FAC2)*O_COSH(I)
       FSS  = SINH(KD(I)*FAC2)*O_SINH(I)

       CFF1 = WAVE_NUMBER(I)*WAVE_ENERGY(I)
!       CFF4 = CFF1*FCS*FSS
       CFF4 = CFF1*FSC*FSS
       CFF6 = CFF1*FCS*FCC
       CFF5 = CFF1*FCS*FCC*O_WAVE_NUMBER(I)*O_WAVE_NUMBER(I)
#if defined(WAVE_ROLLER)
       CFF3 = 1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4)
       CFF3 = CFF3*OROLLER(I)*ROLLA(I)/(WLEN(I)+eps1)*WAVE_C(I)**2
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)-CFF4 + &
                  + CFF6 + CFF3*COS_DIR(I)*COS_DIR(I)
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)-CFF4 + &
                  + CFF6 + CFF3*SIN_DIR(I)*SIN_DIR(I)
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)      + &
                  CFF3*SIN_DIR(I)*COS_DIR(I)
#else
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)+CFF6-CFF4
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)+CFF6-CFF4
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)
#endif       
       TPZDIST(I,K) = FCC*FSS
     END DO  
    ELSE
     DO K=1,KBM1
       FAC2 = ZZ(I,K)
       FCC  = EXP(KD(I)*FAC2)
       FCS  = FCC
       FSC  = FCC
       FSS  = FCC

       CFF1 = WAVE_NUMBER(I)*WAVE_ENERGY(I)
!       CFF4 = CFF1*FCS*FSS
       CFF4 = CFF1*FSC*FSS
       CFF6 = CFF1*FCS*FCC
       CFF5 = CFF1*FCS*FCC*O_WAVE_NUMBER(I)*O_WAVE_NUMBER(I)
#if defined(WAVE_ROLLER)
       CFF3 = 1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4)
       CFF3 = CFF3*OROLLER(I)*ROLLA(I)/(WLEN(I)+eps1)*WAVE_C(I)**2
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)-CFF4 + &
                  + CFF6 + CFF3*COS_DIR(I)*COS_DIR(I)
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)-CFF4 + &
                  + CFF6 + CFF3*SIN_DIR(I)*SIN_DIR(I)
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)      + &
                  CFF3*SIN_DIR(I)*COS_DIR(I)
#else
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)+CFF6-CFF4
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)+CFF6-CFF4
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)
#endif       
       TPZDIST(I,K) = FCC*FSS
     END DO  
    END IF
   END DO  

#if defined (MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,SXX,SXY,SYY)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SXX,SXY,SYY)   !Jianzhong
#endif

!JQI       allocate(sxx_gl(0:MGL)); sxx_gl = 0.0_sp
!JQI       allocate(sxy_gl(0:MGL)); sxy_gl = 0.0_sp
!JQI       allocate(syy_gl(0:MGL)); syy_gl = 0.0_sp
!JQI       allocate(sxx_tmp(0:MT)); sxx_tmp = 0.0_sp
!JQI       allocate(sxy_tmp(0:MT)); sxy_tmp = 0.0_sp
!JQI       allocate(syy_tmp(0:MT)); syy_tmp = 0.0_sp

 !JQI      sxx_tmp(:) = sxx(:,1)
!JQI       sxy_tmp(:) = sxy(:,1)
!JQI       syy_tmp(:) = syy(:,1)

!JQI       CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,sxx_tmp,sxx_gl)
!JQI       CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,sxy_tmp,sxy_gl)
!JQI       CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,syy_tmp,syy_gl)
!JQI       do i=1,mgl
!JQI         write(700+myid,*) i,sxx_gl(i),sxy_gl(i),syy_gl(i)
!JQI       end do

!JQI       deallocate(sxx_gl)
!JQI       deallocate(sxy_gl)
!JQI       deallocate(syy_gl)
!JQI       deallocate(sxx_tmp)
!JQI       deallocate(sxy_tmp)
!JQI       deallocate(syy_tmp)
!JQI

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

#    if defined (WET_DRY)
     IF(ISWETC(IA) == 1 .OR. ISWETC(IB) == 1)THEN
#    endif

     DO K=1,KBM1
       SXXIJ=0.5_SP*(SXX(J1,K)+SXX(J2,K))
       SXYIJ=0.5_SP*(SXY(J1,K)+SXY(J2,K))
       SYYIJ=0.5_SP*(SYY(J1,K)+SYY(J2,K))
       DIJ = 0.5_SP*(D(J1)*DZ(J1,K)+D(J2)*DZ(J2,K))

#      if defined (SPHERICAL)
       !for spherical coordinator and domain across 360^o          
       XTMP  = VX(J2)*TPI-VX(J1)*TPI
       XTMP1 = VX(J2)-VX(J1)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF

       EXFLUX       = DIJ*SXXIJ*DLTYC(I)
       PSXXPX(IA,K) = PSXXPX(IA,K) - EXFLUX
       PSXXPX(IB,K) = PSXXPX(IB,K) + EXFLUX

       EXFLUX       = DIJ*SXYIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSXYPY(IA,K) = PSXYPY(IA,K) + EXFLUX
       EXFLUX       = DIJ*SXYIJ*XTMP*COS(DEG2RAD*YC(IB))
       PSXYPY(IB,K) = PSXYPY(IB,K) - EXFLUX

       EXFLUX     = DIJ*SXYIJ*DLTYC(I)
       PSXYPX(IA,K) = PSXYPX(IA,K) - EXFLUX
       PSXYPX(IB,K) = PSXYPX(IB,K) + EXFLUX

       EXFLUX       = DIJ*SYYIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSYYPY(IA,K) = PSYYPY(IA,K) + EXFLUX
       EXFLUX       = DIJ*SYYIJ*XTMP*COS(DEG2RAD*YC(IB))
       PSYYPY(IB,K) = PSYYPY(IB,K) - EXFLUX

#      else
       EXFLUX       = DIJ*SXXIJ*DLTYC(I)
       PSXXPX(IA,K) = PSXXPX(IA,K) - EXFLUX
       PSXXPX(IB,K) = PSXXPX(IB,K) + EXFLUX

       EXFLUX       = DIJ*SXYIJ*DLTXC(I)
       PSXYPY(IA,K) = PSXYPY(IA,K) + EXFLUX
       PSXYPY(IB,K) = PSXYPY(IB,K) - EXFLUX

       EXFLUX       = DIJ*SXYIJ*DLTYC(I)
       PSXYPX(IA,K) = PSXYPX(IA,K) - EXFLUX
       PSXYPX(IB,K) = PSXYPX(IB,K) + EXFLUX

       EXFLUX       = DIJ*SYYIJ*DLTXC(I)
       PSYYPY(IA,K) = PSYYPY(IA,K) + EXFLUX
       PSYYPY(IB,K) = PSYYPY(IB,K) - EXFLUX
#      endif     
     END DO
#    if defined (WET_DRY)
     END IF
#    endif
   END DO

   WAVESTRX_3D = 0.0_SP
   WAVESTRY_3D = 0.0_SP
   
   CALL RADIATION_STRESS_Z(WAVE_ENERGY,KD,KDMAX,PSPXPZ,PSPYPZ)

   WAVESTRX_3D = PSXXPX + PSXYPY - PSPXPZ
   WAVESTRY_3D = PSXYPX + PSYYPY - PSPYPZ

!qxu set rediation stress limit to 200 Pa 01/19/2021
     WAVESTRX_3D = max(min(WAVESTRX_3D,200.0_SP),-200.0_SP)
     WAVESTRY_3D = max(min(WAVESTRY_3D,200.0_SP),-200.0_SP)
!qxu}

#  if defined (WET_DRY)
   DO I = 1,NT
     WAVESTRX_3D(I,:) = WAVESTRX_3D(I,:)*ISWETC(I)
     WAVESTRY_3D(I,:) = WAVESTRY_3D(I,:)*ISWETC(I)

!JQI
     IF(H1(I) <= 0.0_sp)THEN
       WAVESTRX_3D(I,:) = 0.0_SP
       WAVESTRY_3D(I,:) = 0.0_SP
     END IF
!JQI    
     IF(ISBCE(I) == 2)THEN
       WAVESTRX_3D(I,:) = 0.0_SP
       WAVESTRY_3D(I,:) = 0.0_SP
     END IF
   END DO  
#  endif
#  if !defined (TWO_D_MODEL)
   WAVESTRX_2D(:) = 0.0_SP; WAVESTRY_2D(:) = 0.0_SP
   DO I = 1,NT
     DO K=1,KBM1
        WAVESTRX_2D(I) = WAVESTRX_2D(I)+WAVESTRX_3D(I,K)
        WAVESTRY_2D(I) = WAVESTRY_2D(I)+WAVESTRY_3D(I,K)
     END DO
   END DO
#  else
   CALL RADIATION_STRESS_2D
#  endif
   WAVESTRX_2D = WAVESTRX_2D*RAMP
   WAVESTRY_2D = WAVESTRY_2D*RAMP
   WAVESTRX_3D = WAVESTRX_3D*RAMP
   WAVESTRY_3D = WAVESTRY_3D*RAMP

#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WAVESTRX_3D,WAVESTRY_3D) 
#  endif
!Calculate stokes velocity
   U_STOKES_3D_TMP = 0.0_SP; U_STOKES_3D = 0.0_SP; U_STOKES_2D = 0.0_SP
   V_STOKES_3D_TMP = 0.0_SP; V_STOKES_3D = 0.0_SP; V_STOKES_2D = 0.0_SP
   DO I=1,M
    IF(KD(I) <= KDMAX)THEN
     DO K=1,KBM1
        FAC2 = 1.0_SP+ZZ(I,K)
# if defined (WAVE_ROLLER)
        CFF2=2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1))
# else
        CFF2=2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*WAVE_ENERGY(I)
# endif
        U_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_X(I)
        V_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_Y(I)
     ENDDO
    ELSE
     DO K=1,KBM1
        FAC2 = ZZ(I,K)
# if defined (WAVE_ROLLER)
        CFF2=2/WAVE_C(I)*EXP(KD(I)*FAC2)*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1))
# else
        CFF2=2/WAVE_C(I)*EXP(KD(I)*FAC2)*WAVE_ENERGY(I)
# endif
        U_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_X(I)
        V_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_Y(I)
     ENDDO
    END IF 
   ENDDO
   DO I=1,NT
     DO K=1,KBM1
       U_STOKES_3D(I,K)=(U_STOKES_3D_TMP(NV(I,1),K)+U_STOKES_3D_TMP(NV(I,2),K)+U_STOKES_3D_TMP(NV(I,3),K))/3.0_SP
       V_STOKES_3D(I,K)=(V_STOKES_3D_TMP(NV(I,1),K)+V_STOKES_3D_TMP(NV(I,2),K)+V_STOKES_3D_TMP(NV(I,3),K))/3.0_SP
       U_STOKES_2D(I)=U_STOKES_2D(I)+U_STOKES_3D(I,K)*DZ1(I,K)
       V_STOKES_2D(I)=V_STOKES_2D(I)+V_STOKES_3D(I,K)*DZ1(I,K)
     ENDDO
   ENDDO
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U_STOKES_3D,V_STOKES_3D) 
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U_STOKES_2D,V_STOKES_2D) 
#  endif

   RETURN
   END SUBROUTINE RADIATION_STRESS  
!==============================================================================|

   SUBROUTINE RADIATION_STRESS_Z(WAVE_ENERGY,KD,KDMAX,PSPXPZ,PSPYPZ) 

!==============================================================================|
   USE MOD_PREC
   USE ALL_VARS, ONLY : Z,Z1,VX,VY,N2E2D

   IMPLICIT NONE
   REAL(SP), INTENT(IN)  :: WAVE_ENERGY(0:MT),KD(0:MT),KDMAX
   REAL(SP), INTENT(OUT) :: PSPXPZ(0:NT,KB),PSPYPZ(0:NT,KB)
   REAL(SP)              :: SPX(KB),SPY(KB)

   REAL(SP) :: WAVE_ENERGY1(0:NT), KD1(0:NT)
   REAL(SP), DIMENSION(0:MT) :: O_COSH,O_SINH
   INTEGER  :: I,K,J,J1,J2,I1,I2,I3
   REAL(SP) :: FSS1,FCS1,FSC1,FCC1
   REAL(SP) :: CFF1,CFF2,FAC1,FAC2,FAC3
   REAL(SP) :: WEIJ,KDIJ,DIJ,SIJ
#  if defined (SPHERICAL)
   REAL(SP) :: XTMP,XTMP1
#  endif
!==============================================================================|

!----------INITIALIZE ARRAYS---------------------------------------------------!
    CALL N2E2D(WAVE_ENERGY,WAVE_ENERGY1)
    CALL N2E2D(KD,KD1)
    PSPXPZ  = 0.0_SP   ;PSPYPZ  = 0.0_SP  
    O_COSH  = 1.0_SP/COSH(KD)
    O_SINH  = 1.0_SP/SINH(KD)

   DO I = 1, N
     SPX = 0.0_SP; SPY = 0.0_SP
#    if defined (WET_DRY)
     IF(ISWETCT(I)*ISWETC(I) == 1)THEN
#    endif
       I1=NV(I,1);I2=NV(I,2);I3=NV(I,3)
       IF(KD1(I) <= KDMAX)THEN
       
       DO K=1,KBM1
! Calculate some coefficients
         FAC1 = 1.0_SP+Z(I1,K) 
         FAC2 = 1.0_SP+Z(I2,K)
         FAC3 = 1.0_SP+Z(I3,K)
         FCC1 = (COSH(KD(I1)*FAC1)*O_COSH(I1)+COSH(KD(I2)*FAC2)*O_COSH(I2)+COSH(KD(I3)*FAC3)*O_COSH(I3))/3
         FCS1 = (COSH(KD(I1)*FAC1)*O_SINH(I1)+COSH(KD(I2)*FAC2)*O_SINH(I2)+COSH(KD(I3)*FAC3)*O_SINH(I3))/3
         FSC1 = (SINH(KD(I1)*FAC1)*O_COSH(I1)+SINH(KD(I2)*FAC2)*O_COSH(I2)+SINH(KD(I3)*FAC3)*O_COSH(I3))/3
         FSS1 = (SINH(KD(I1)*FAC1)*O_SINH(I1)+SINH(KD(I2)*FAC2)*O_SINH(I2)+SINH(KD(I3)*FAC3)*O_SINH(I3))/3
         CFF1=(FCC1-FSS1)*(FSS1*0.5_SP)
         CFF2=(FCC1-FSS1)*(FCS1*(1+Z1(I,K))*WAVE_ENERGY1(I)-WAVE_ENERGY1(I)*FSS1/TANH(KD1(I)))
         DO J = 1, 3
           J1=J+1-INT((J+1)/4)*3
           J2=J+2-INT((J+2)/4)*3
           WEIJ=0.5_SP*(WAVE_ENERGY(NV(I,J1))+WAVE_ENERGY(NV(I,J2)))*CFF1
           KDIJ=0.5_SP*(KD(NV(I,J1))+KD(NV(I,J2)))*CFF2
           SIJ=WEIJ+KDIJ
#          if defined (SPHERICAL)
           SPX(K)=SPX(K)-DELTUY(I,J)*SIJ
#          else
           SPX(K)=SPX(K)-(VY(NV(I,J2))-VY(NV(I,J1)))*SIJ
#          endif

#          if defined (SPHERICAL)
           XTMP  = VX(NV(I,J2))*TPI-VX(NV(I,J1))*TPI
           XTMP1 = VX(NV(I,J2))-VX(NV(I,J1))
           IF(XTMP1 >  180.0_SP)THEN
             XTMP = -360.0_SP*TPI+XTMP
           ELSE IF(XTMP1 < -180.0_SP)THEN
             XTMP =  360.0_SP*TPI+XTMP
           END IF  

           SPY(K)=SPY(K)+XTMP*COS(DEG2RAD*YC(I))*SIJ
#          else
           SPY(K)=SPY(K)+(VX(NV(I,J2))-VX(NV(I,J1)))*SIJ
#          endif
         END DO
       END DO
       
       ELSE
       
       DO K=1,KBM1
! Calculate some coefficients
         FAC1 = Z(I1,K) 
         FAC2 = Z(I2,K)
         FAC3 = Z(I3,K)
         FCC1 = (EXP(KD(I1)*FAC1)+EXP(KD(I2)*FAC2)+EXP(KD(I3)*FAC3))/3.0_SP
	 FCS1 = FCC1
         FSC1 = FCC1
         FSS1 = FCC1
	 
         CFF1=(FCC1-FSS1)*(FSS1*0.5_SP)
         CFF2=(FCC1-FSS1)*(FCS1*(1+Z1(I,K))*WAVE_ENERGY1(I)-WAVE_ENERGY1(I)*FSS1)
         DO J = 1, 3
           J1=J+1-INT((J+1)/4)*3
           J2=J+2-INT((J+2)/4)*3
           WEIJ=0.5_SP*(WAVE_ENERGY(NV(I,J1))+WAVE_ENERGY(NV(I,J2)))*CFF1
           KDIJ=0.5_SP*(KD(NV(I,J1))+KD(NV(I,J2)))*CFF2
           SIJ=WEIJ+KDIJ
#          if defined (SPHERICAL)
           SPX(K)=SPX(K)-DELTUY(I,J)*SIJ
#          else
           SPX(K)=SPX(K)-(VY(NV(I,J2))-VY(NV(I,J1)))*SIJ
#          endif

#          if defined (SPHERICAL)
           XTMP  = VX(NV(I,J2))*TPI-VX(NV(I,J1))*TPI
           XTMP1 = VX(NV(I,J2))-VX(NV(I,J1))
           IF(XTMP1 >  180.0_SP)THEN
             XTMP = -360.0_SP*TPI+XTMP
           ELSE IF(XTMP1 < -180.0_SP)THEN
             XTMP =  360.0_SP*TPI+XTMP
           END IF  

           SPY(K)=SPY(K)+XTMP*COS(DEG2RAD*YC(I))*SIJ
#          else
           SPY(K)=SPY(K)+(VX(NV(I,J2))-VX(NV(I,J1)))*SIJ
#          endif
         END DO
       END DO
       
       END IF
       
       DO K = 1,KBM1
         PSPXPZ(I,K) = SPX(K)-SPX(K+1) 
	 PSPYPZ(I,K) = SPY(K)-SPY(K+1)
       END DO
#    if defined (WET_DRY)
     END IF
#    endif
   END DO

   RETURN
   END SUBROUTINE RADIATION_STRESS_Z
!==============================================================================|

end module
