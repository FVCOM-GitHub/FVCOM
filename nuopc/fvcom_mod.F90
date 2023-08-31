!
!! @mainpage FVCOM NUOPC Cap
!! @author Jianhua Qi based on Saeed Moghimi (moghimis@gmail.com)
!! @date 02/25/2020 Original documentation
!------------------------------------------------------

module fvcom_mod

  !-----------------------------------------------------------------------------
  ! FVCOM mesh utility
  !-----------------------------------------------------------------------------
  use mpi
  use ESMF
  use NUOPC

  use MOD_DRIVER, only : fvcom_pet_num
  use LIMS   ,    only : MGL,MT,M,NGL,NT,N,MYID,MSRID,NPROCS    !np,ne,nm,slam,sfea
  use ALL_VARS,   only : NV,VX,VY,H,LAT,LON,NVG,LATC,LONC
  use MOD_PAR,    only : NGID_X,EGID_X,NLID,NMAP,ACOLLECT
  use MOD_PAR,    only : NGID, EL_PID, NHE, HE_LST, EGID, ELID_X,HE_OWN,NLID_X,NLID
  use CONTROL,    only : MPI_FVCOM_GROUP
!  use GLOBAL , only: IMAP_EL_LG,NODES_LG
!  use GLOBAL , only: ETA2, UU2, VV2  ! Export water level and velocity fileds to wave model
!  USE GLOBAL,  ONLY: RSNX2, RSNY2    ! Import wave 2D forces from wave model
!  use SIZES  , only: ROOTDIR

  implicit none
  type meshdata
    !> \details vm is an ESMF_VM object.  ESMF_VM is just an ESMF virtual machine class,
    !! which we will use to get the data about the local PE and PE count.
    type(ESMF_VM)                      :: vm
    !> \details This array contains the node coordinates of the mesh. For
    !! example, in a 2D mesh, the \c jth coordinate of the \c nth node
    !!oords
    real(ESMF_KIND_R8), allocatable    :: NdCoords(:)
    !> \details This array contains the elevation of different nodes of the mesh
    real(ESMF_KIND_R8), allocatable    :: bathymetry(:)
    !> \details Number of nodes present in the current PE. This is different from the
    !! number of nodes owned by this PE (cf. NumOwnedNd)
    integer(ESMF_KIND_I4)              :: NumNd
    !> \details Number of nodes owned by this PE. This is different from the number of
    !! nodes present in the current PE (cf. NumNd)
    integer(ESMF_KIND_I4)              :: NumOwnedNd
    !> \details Number of elements in the current PE. This includes ghost elements and
    !! owned elements. However, we do not bother to distinguish between owned
    !! element and present element (as we did for the nodes).
    integer(ESMF_KIND_I4)              :: NumEl
    !> \details Number of nodes of each element, which is simply three in 2D ADCIRC.
    integer(ESMF_KIND_I4)              :: NumND_per_El
    !> \details Global node numbers of the nodes which are present in the current PE.
    integer(ESMF_KIND_I4), allocatable :: NdIDs(:)
    !> \details Global element numbers which are present in the current PE.
    integer(ESMF_KIND_I4), allocatable :: ElIDs(:)
    !> \details The element connectivity array, for the present elements in the current PE.
    !! The node numbers are the local numbers of the present nodes. All the element
    !! connectivities are arranged in this one-dimensional array.
    integer(ESMF_KIND_I4), allocatable :: ElConnect(:)
    !> \details The number of the PE's which own each of the nodes present this PE.
    !! This number is zero-based.
    integer(ESMF_KIND_I4), allocatable :: NdOwners(:)
    !> \details An array containing the element types, which are all triangles in our
    !! application.
    integer(ESMF_KIND_I4), allocatable :: ElTypes(:)
    !> \details This array contains the element coordinates of the mesh. 
    real(ESMF_KIND_R8), allocatable    :: ElCoords(:)
    !> \details This is an array, which maps the indices of the owned nodes to the indices of the present
    !! nodes. For example, assume we are on <tt>PE = 1</tt>, and we have four nodes present, and the
    !! first and third nodes belong to <tt>PE = 0</tt>. So we have:
    !! \code
    !! NumNd = 4
    !! NumOwnedNd = 2
    !! NdOwners = (/0, 1, 0, 1/)
    !! NdIDs = (/2, 3, 5, 6/)
    !! owned_to_present = (/2, 4/)    <-- Because the first node owned by this PE is actually
    !!                                    the second node present on this PE, and so on.
    !! \endcode
    integer(ESMF_KIND_I4), allocatable :: owned_to_present_nodes(:)
    integer(ESMF_KIND_I4), allocatable :: owned_to_present_halo_nodes(:)

    !! number of nodes owned by this PE
    integer(ESMF_KIND_I4)              :: NumOwnedNd_NOHALO
    integer(ESMF_KIND_I4)              :: NumOwnedNd_HALO
  end type meshdata

  ! reading data time management info WW3 <-----> FVCOM exchange
  integer :: fvcom_cpl_int,fvcom_cpl_num,fvcom_cpl_den

  REAL,ALLOCATABLE,TARGET :: FVCOM_SXX2(:,:)
  REAL,ALLOCATABLE,TARGET :: FVCOM_SXY2(:,:)
  REAL,ALLOCATABLE,TARGET :: FVCOM_SYY2(:,:)

  REAL,ALLOCATABLE,TARGET :: FVCOM_WHS(:,:)
  REAL,ALLOCATABLE,TARGET :: FVCOM_WLEN(:,:)
  REAL,ALLOCATABLE,TARGET :: FVCOM_WDIR(:,:)

  REAL,ALLOCATABLE,TARGET :: FVCOM_ZETA(:,:)
  REAL,ALLOCATABLE,TARGET :: FVCOM_VELX(:,:)
  REAL,ALLOCATABLE,TARGET :: FVCOM_VELY(:,:)

  REAL,ALLOCATABLE,TARGET :: FVCOM_WVNX(:,:)
  REAL,ALLOCATABLE,TARGET :: FVCOM_WVNY(:,:)
  REAL,ALLOCATABLE,TARGET :: FVCOM_PRN(:,:)

  logical :: NUOPC4WAV,NUOPC4MET

  ! module name
  character(*), parameter :: modName = "(fvcom_mod)"

  !-----------------------------------------------------------------------------
  contains
  !> \author Ali Samii - 2016
  !! See: https://github.com/samiiali
  !> @details Using the data available in <tt> fort.14, fort.18, partmesh.txt</tt> files
  !! this function extracts the scalars and arrays required for construction of a
  !! meshdata object.
  !! After calling this fucntion, one can call create_parallel_esmf_mesh_from_meshdata()
  !! or create_masked_esmf_mesh_from_data() to create an ESMF_Mesh.
  !! @param vm This is an ESMF_VM object, which will be used to obtain the \c localPE
  !! and \c peCount of the \c MPI_Communicator.
  !! @param global_fort14_dir This is the directory path (relative to the executable
  !! or an absolute path) which contains the global \c fort.14 file (not the fort.14
  !! after decomposition).
  !! @param the_data This is the output meshdata object.
  !!

  !> \details As the name of this function suggests, this funciton creates a parallel
  !! ESMF_Mesh from meshdata object. This function should be called collectively by
  !! all PEs for the parallel mesh to be created. The function, extract_parallel_data_from_mesh()
  !! should be called prior to calling this function.
  !! \param the_data This the input meshdata object.
  !! \param out_esmf_mesh This is the ouput ESMF_Mesh object.
  subroutine create_parallel_esmf_mesh_from_meshdata(the_data, out_esmf_mesh)
    implicit none
    type(ESMF_Mesh), intent(out) :: out_esmf_mesh
    type(meshdata), intent(in)   :: the_data
    integer, parameter           :: dim1=2, spacedim=2, NumND_per_El=3
    type(ESMF_Distgrid)          :: nodeDistgrid, elementDistgrid
    integer                      :: rc

    ! create node distgrid
    nodeDistgrid = ESMF_DistgridCreate(the_data%NdIDs, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    ! create element distgrid
    elementDistgrid = ESMF_DistgridCreate(the_data%ElIDs, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    ! create mesh
    out_esmf_mesh=ESMF_MeshCreate(parametricDim=dim1, spatialDim=spacedim, &
        nodeIDs=abs(the_data%NdIDs), &
        nodeCoords=the_data%NdCoords, &
        nodeOwners=the_data%NdOwners, &
        nodalDistgrid=nodeDistgrid, &
        elementIDs=abs(the_data%ElIDs), &
        elementTypes=the_data%ElTypes, &
        elementConn=the_data%ElConnect, &
        elementCoords=the_data%ElCoords, &
        elementDistgrid=elementDistgrid, &
        coordSys=ESMF_COORDSYS_SPH_DEG, &
        rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

  end subroutine create_parallel_esmf_mesh_from_meshdata
  !
  !> \author Ali Samii - 2016
  !! See: https://github.com/samiiali
  !> @details Using the data available in <tt> fort.14, fort.18, partmesh.txt</tt> files
  !! this function extracts the scalars and arrays required for construction of a
  !! meshdata object.
  !! After calling this fucntion, one can call create_parallel_esmf_mesh_from_meshdata()
  !! or create_masked_esmf_mesh_from_data() to create an ESMF_Mesh.
  !! @param vm This is an ESMF_VM object, which will be used to obtain the \c localPE
  !! and \c peCount of the \c MPI_Communicator.
  !! @param global_fort14_dir This is the directory path (relative to the executable
  !! or an absolute path) which contains the global \c fort.14 file (not the fort.14
  !! after decomposition).
  !! @param the_data This is the output meshdata object.
  !!
  subroutine extract_parallel_data_from_mesh(the_data,localPet)
    implicit none

    type(meshdata), intent(inout) :: the_data
!JQI        character(len=*), intent(in)          :: global_fort14_dir
    integer, intent(in)           :: localPet
    character(len=6)              :: PE_ID, garbage1
        
    character(len=200)            :: fort14_filename, fort18_filename, partmesh_filename
    integer                       :: i1, i11, j1,j2, i_num, petCount, num_global_nodes, rc,ierr
    integer, allocatable          :: local_node_numbers(:), local_elem_numbers(:), node_owner(:),node_owner_gl(:)
    integer, parameter            :: dim1=2, NumND_per_El=3

    integer(ESMF_KIND_I4), allocatable :: tmp_owned_to_present_nodes(:)
    integer(ESMF_KIND_I4), allocatable :: tmp_owned_to_present_halo_nodes(:)
!        write(PE_ID, "(A,I4.4)") "PE", localPet
!        fort14_filename = TRIM(global_fort14_dir)//'/'//PE_ID//"/fort.14"
!        fort18_filename = TRIM(global_fort14_dir)//'/'//PE_ID//"/fort.18"
!        partmesh_filename = TRIM(global_fort14_dir)//'/'//'partmesh.txt'

!        open(unit=23514, file=fort14_filename, form='FORMATTED', status='OLD', action='READ')
!        open(unit=23518, file=fort18_filename, form='FORMATTED', status='OLD', action='READ')
!        open(unit=235100, file=partmesh_filename, form='FORMATTED', status='OLD', action='READ')

!        read(unit=23514, fmt=*)
!        read(unit=23514, fmt=*) the_data%NumEl, the_data%NumNd

    the_data%NumEl = NT 
    the_data%NumNd = MT
    the_data%NumOwnedND = M

    print*, 'NT,N,MT,M: ', NT, N, MT, M
    allocate(the_data%NdIDs(the_data%NumNd))
    allocate(local_node_numbers(the_data%NumNd))
    allocate(the_data%ElIDs(the_data%NumEl))
    allocate(local_elem_numbers(the_data%NumEl))
    allocate(the_data%NdCoords(dim1*the_data%NumNd))
    allocate(the_data%bathymetry(the_data%NumNd))
    allocate(the_data%ElConnect(NumND_per_El*the_data%NumEl))
    allocate(the_data%NdOwners(the_data%NumNd))
    allocate(the_data%ElTypes(the_data%NumEl))
    allocate(the_data%ElCoords(dim1*the_data%NumEl))

!        read(unit=23518, fmt=*)
!        read(unit=23518, fmt=*)
!        read(unit=23518, fmt=*) local_elem_numbers
!        the_data%ElIDs = abs(local_elem_numbers)
!        read(unit=23518, fmt=*) garbage1, num_global_nodes, garbage2, garbage3
!        read(unit=23518, fmt=*) local_node_numbers
!        the_data%NumOwnedND = 0
!        do i1 = 1, the_data%NumNd, 1
!            if (local_node_numbers(i1) > 0) then
!                the_data%NumOwnedND = the_data%NumOwnedND + 1
!            end if
!        end do
    do i1 = 1,the_data%NumEl
      local_elem_numbers(i1) = EGID_X(i1)
    end do  
    the_data%ElIDs = local_elem_numbers
    do i1 = 1,the_data%NumNd
      local_node_numbers(i1) = NGID_X(i1)
    end do
    the_data%NdIDs = local_node_numbers

    num_global_nodes = MGL
    allocate(node_owner_gl(num_global_nodes))
    allocate(node_owner(m))
!JQI        allocate(the_data%owned_to_present_nodes(the_data%NumOwnedND))
!JQI        allocate(the_data%owned_to_present_halo_nodes(m-the_data%NumOwnedND))
    allocate(tmp_owned_to_present_nodes(the_data%NumOwnedND))
    allocate(tmp_owned_to_present_halo_nodes(the_data%NumOwnedND))

!   read(unit=235100, fmt=*) node_owner
    node_owner = 0
    i11 = 0
    do i1 = 1,num_global_nodes
      if(NLID(i1) /= 0)then
        i11 = i11+1
        node_owner(i11) = MYID
!        write(100+myid,*) node_owner(i11),i11
      end if
    end do  

    if(fvcom_pet_num == 1)then
      node_owner_gl = node_owner
    else
      call acollect(MYID,MSRID,NPROCS,NMAP,NODE_OWNER,NODE_OWNER_GL)
      call ESMF_VMBroadcast(the_data%vm, NODE_OWNER_GL,size(NODE_OWNER_GL),0,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
!JQI     CALL MPI_BCAST(NODE_OWNER_GL,MGL,MPI_INTEGER,0,MPI_FVCOM_GROUP,IERR)
!!!JQI        CALL ESMF_VMBroadcast(the_data%vm, NODE_OWNER_GL,size(NODE_OWNER_GL),0,rc=rc)
    end if

!        do i1 = 1, the_data%NumNd
!            read(unit=23514, fmt=*) local_node_numbers(i1), &
!                the_data%NdCoords((i1-1)*dim1 + 1), &
!                the_data%NdCoords((i1-1)*dim1 + 2), &
!                the_data%bathymetry(i1)
!        end do
    do i1 = 1, the_data%NumNd
       the_data%NdCoords((i1-1)*dim1 + 1) = lon(i1)
       the_data%NdCoords((i1-1)*dim1 + 2) = lat(i1)
       the_data%bathymetry(i1) = h(i1)
    end do
        
!        do i1 = 1, the_data%NumEl
!            read(unit=23514, fmt=*) local_elem_numbers(i1), i_num, &
!                the_data%ElConnect((i1-1)*NumND_per_El+1), &
!                the_data%ElConnect((i1-1)*NumND_per_El+2), &
!                the_data%ElConnect((i1-1)*NumND_per_El+3)
!        end do
    do i1 = 1, the_data%NumEl
      the_data%ElConnect((i1-1)*NumND_per_El+1) = NV(i1,1)
      the_data%ElConnect((i1-1)*NumND_per_El+2) = NV(i1,2)
      the_data%ElConnect((i1-1)*NumND_per_El+3) = NV(i1,3)
!      write(500+myid,*) i1,nv(i1,1),nv(i1,2),nv(i1,3)
    end do

    do i1 = 1, the_data%NumEl
       the_data%ElCoords((i1-1)*dim1 + 1) = lonc(i1)
       the_data%ElCoords((i1-1)*dim1 + 2) = latc(i1)
    end do

    do i1= 1, the_data%NumNd
      the_data%NdOwners(i1) = node_owner_gl(the_data%NdIDs(i1)) - 1
!      write(200+myid,*) i1,the_data%NdOwners(i1)
    end do

    j1 = 0
    j2 = 0
    do i1 = 1, the_data%NumOwnedNd
      if (the_data%NdOwners(i1) == localPet) then
        j1 = j1 + 1
!JQI                the_data%owned_to_present_nodes(j1) = i1
        tmp_owned_to_present_nodes(j1) = i1
!JQI		write(300+myid,*) j1,the_data%owned_to_present_nodes(j1)
!        write(300+myid,*) j1,tmp_owned_to_present_nodes(j1)
      else if (the_data%NdOwners(i1) /= localPet) then
        j2 = j2 + 1
!JQI                the_data%owned_to_present_halo_nodes(j2) = i1
        tmp_owned_to_present_halo_nodes(j2) = i1
!        write(400+myid,*) j2,i1,ngid_x(i1)
      end if
    end do
    the_data%ElTypes = ESMF_MESHELEMTYPE_TRI

    the_data%NumOwnedNd_NOHALO = J1
    the_data%NumOwnedNd_HALO = J2
    allocate(the_data%owned_to_present_nodes(j1))
    allocate(the_data%owned_to_present_halo_nodes(j2))
    the_data%owned_to_present_nodes(1:j1) = tmp_owned_to_present_nodes(1:j1)
    if(j2 > 0)the_data%owned_to_present_halo_nodes(1:j2) = tmp_owned_to_present_halo_nodes(1:j2)
        
!	if(j2 /= m-the_data%NumOwnedND)then
!	  print*,"Not match ...... ",j2,m-the_data%NumOwnedND
!	  stop
!	end if  
!    print*,"Match ...... ",j2,the_data%NumOwnedND,myid

!        close(23514)
!        close(23518)
!        close(235100)
  end subroutine extract_parallel_data_from_mesh

  subroutine read_config()
  ! This subroutine is not used with NEMS system. Because the 
  ! time interval information is passed via nems.configure file 
  ! with time slot definitation.

    implicit none   
    character(ESMF_MAXPATHLEN) :: fname ! config file name
    type(ESMF_Config)          :: cf     ! the Config itself
    integer                    :: rc

    rc = ESMF_SUCCESS
    
    !Initiate reading resource file
    cf = ESMF_ConfigCreate(rc=rc)  ! Create the empty Config
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    fname = "config.rc" ! Name the Resource File
    call ESMF_ConfigLoadFile(cf, fname, rc=rc) ! Load the Resource File
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! read time coupling interval info
   
!JQI    call ESMF_ConfigGetAttribute(cf, adc_cpl_int, label="cpl_int:",default=300, rc=rc)
!JQI    call ESMF_ConfigGetAttribute(cf, adc_cpl_num, label="cpl_num:",default=0  , rc=rc)
!JQI    call ESMF_ConfigGetAttribute(cf, adc_cpl_den, label="cpl_den:",default=1  , rc=rc)
!JQI    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!JQI      line=__LINE__, &
!JQI      file=__FILE__)) &
!JQI      return  ! bail out
        
    call ESMF_ConfigDestroy(cf, rc=rc) ! Destroy the Config
        
  end subroutine read_config

  subroutine eliminate_ghosts(the_data, localPet, dbug)
    implicit none

    type(meshdata), intent(inout) :: the_data
    integer, intent(in)           :: localPet
    logical, intent(in)           :: dbug

    ! local variables
    integer :: i, j, k, indx, indx2
    integer :: n1, i1, j1, pe1, n2, i2, j2, pe2, n3, i3, j3, pe3
    integer :: minNdIDs, maxNdIDs
    integer :: NumEl_nog, NumNd_nog
    integer :: NumND_per_El = 3
    integer :: dim1 = 2
    integer, allocatable  :: ElIDs(:)
    integer, allocatable  :: ElConnect(:)
    real(ESMF_KIND_R8), allocatable :: ElCoords(:)
    integer, allocatable  :: ElMask(:)
    integer, allocatable  :: NdIDs(:)
    real(ESMF_KIND_R8), allocatable :: NdCoords(:)
    integer, allocatable  :: NdOwners(:)
    real(ESMF_KIND_R8), allocatable :: bathymetry(:)
    integer, allocatable  :: NdMask(:)
    character(len=1024)   :: msgString
    character(len=*), parameter :: subname=trim(modName)//':(eliminate_ghosts) '

    ! message for entering call
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! create mask for elements
    if (.not. allocated(ElMask)) allocate(ElMask(the_data%NumEl))
    ElMask = 0

    ! 1 to N are resident elements and N+1 to NT are ghost elements
    do i = 1, N
       ElMask(i) = the_data%ElIDs(i)
    end do

    ! debug
    if (dbug) then
       do i = 1, the_data%NumEl
          ! get node owners
          pe1 = the_data%NdOwners(the_data%ElConnect((i-1)*NumND_per_El+1))
          pe2 = the_data%NdOwners(the_data%ElConnect((i-1)*NumND_per_El+2))
          pe3 = the_data%NdOwners(the_data%ElConnect((i-1)*NumND_per_El+3))

          write(msgString,'(8I8)') the_data%ElIDs(i), ElMask(i), pe1, pe2, pe3, &
            the_data%ElConnect((i-1)*NumND_per_El+1), &
            the_data%ElConnect((i-1)*NumND_per_El+2), &
            the_data%ElConnect((i-1)*NumND_per_El+3)
          call ESMF_LogWrite(subname//' EMASK: '//trim(msgString), ESMF_LOGMSG_INFO)
       end do
    end if

    ElMask = abs(ElMask)

    ! number of elements without ghosts
    NumEl_nog = count(mask=ElMask > 0)
    write(msgString,'(A,I8,A,I8)') 'NumEl = ', size(the_data%ElIDs), ' NumEl_no_ghost = ', NumEl_nog
    if (dbug) call ESMF_LogWrite(subname//' '//trim(msgString), ESMF_LOGMSG_INFO)

    ! allocate temporary arrays for element
    if (.not. allocated(ElIDs)) allocate(ElIDs(NumEl_nog))
    if (.not. allocated(ElConnect)) allocate(ElConnect(NumND_per_El*NumEl_nog)) 
    if (.not. allocated(ElCoords)) allocate(ElCoords(dim1*NumEl_nog))

    ! fill temporary element arrays with non-ghost data
    j = 1
    do i = 1, the_data%NumEl   
       if (ElMask(i) > 0) then
          ElIDs(j) = the_data%ElIDs(i)
          ElConnect((j-1)*NumND_per_El+1) = the_data%ElConnect((i-1)*NumND_per_El+1)
          ElConnect((j-1)*NumND_per_El+2) = the_data%ElConnect((i-1)*NumND_per_El+2)
          ElConnect((j-1)*NumND_per_El+3) = the_data%ElConnect((i-1)*NumND_per_El+3)
          ElCoords((j-1)*dim1+1) = the_data%ElCoords((i-1)*dim1+1)
          ElCoords((j-1)*dim1+2) = the_data%ElCoords((i-1)*dim1+2)
          j = j+1
       end if
    end do
    deallocate(ElMask)

    ! make local node ids in ElConnect monotonic again since we removed
    ! ghosts elements
    minNdIDs = minval(ElConnect, dim=1)
    maxNdIDs = maxval(ElConnect, dim=1)
    write(msgString,'(2I8)') minNdIDs, maxNdIDs
    if (dbug) call ESMF_LogWrite(subname//' minNdIDs, maxNdIDs: '//trim(msgString), ESMF_LOGMSG_INFO)

    ! create mask for unique list of nodes in the element connection
    if (.not. allocated(NdMask)) allocate(NdMask(the_data%NumNd))
    NdMask = 0

    ! 1 to N are resident nodes and N+1 to NT are ghost elements
    do i = minNdIDs, maxNdIDs
       NdMask(i) = i
    end do    

    ! debug, print nodeids and associated masks
    if (dbug) then
       do i = 1, the_data%NumNd
          write(msgString,'(4I8)') the_data%NdIDs(i), NLID_X(the_data%NdIDs(i)), NdMask(i), the_data%NdOwners(i)
          call ESMF_LogWrite(subname//' NMASK: '//trim(msgString), ESMF_LOGMSG_INFO)
       end do
    end if

    ! number of nodes without ghosts
    NumNd_nog = count(mask=NdMask .gt. 0)
    write(msgString,'(A,I8,A,I8)') 'NumNd = ', size(the_data%NdIDs), ' NumNd_no_ghost = ', NumNd_nog
    if (dbug) call ESMF_LogWrite(subname//' '//trim(msgString), ESMF_LOGMSG_INFO)
    write(msgString,'(A,2I8)') 'NdMask min/max = ', minval(NdMask, dim=1, mask=NdMask .gt. 0), &
       maxval(NdMask, dim=1, mask=NdMask .gt. 0)
    if (dbug) call ESMF_LogWrite(subname//' '//trim(msgString), ESMF_LOGMSG_INFO)

    ! allocate temporary arrays for node 
    if (.not. allocated(NdIDs)) allocate(NdIDs(NumNd_nog))
    if (.not. allocated(NdCoords)) allocate(NdCoords(dim1*NumNd_nog))
    if (.not. allocated(NdOwners)) allocate(NdOwners(NumNd_nog))
    if (.not. allocated(bathymetry)) allocate(bathymetry(NumNd_nog))

    ! fill temporary node arrays with non-ghost data
    ! node ids in the node list are global ids
    j = 1
    do i = 1, the_data%NumNd
       if (NdMask(i) .ne. 0) then 
          NdIDs(j) = the_data%NdIDs(i)
          NdCoords((j-1)*dim1+1) = the_data%NdCoords((i-1)*dim1+1)
          NdCoords((j-1)*dim1+2) = the_data%NdCoords((i-1)*dim1+2)
          NdOwners(j) = the_data%NdOwners(i)
          bathymetry(j) = the_data%bathymetry(i)
          j = j+1
       end if
    end do
    deallocate(NdMask)

    ! update data with non-ghost one
    the_data%NumEl = NumEl_nog

    if (allocated(the_data%ElIDs)) then
       deallocate(the_data%ElIDs)
       allocate(the_data%ElIDs(NumEl_nog))
       the_data%ElIDs = ElIDs
       deallocate(ElIDs)
    end if

    if (allocated(the_data%ElConnect)) then
       deallocate(the_data%ElConnect)
       allocate(the_data%ElConnect(NumND_per_El*NumEl_nog))
       the_data%ElConnect = ElConnect
       deallocate(ElConnect)
    end if

    if (allocated(the_data%ElCoords)) then
       deallocate(the_data%ElCoords)
       allocate(the_data%ElCoords(dim1*NumEl_nog))
       the_data%ElCoords = ElCoords
       deallocate(ElCoords)
    end if

    if (allocated(the_data%ElTypes)) then
       deallocate(the_data%ElTypes)
       allocate(the_data%ElTypes(NumEl_nog))
       the_data%ElTypes = ESMF_MESHELEMTYPE_TRI
    end if

    the_data%NumNd = NumNd_nog

    if (allocated(the_data%NdIDs)) then
       deallocate(the_data%NdIDs)
       allocate(the_data%NdIDs(NumNd_nog))
       the_data%NdIDs = NdIDs
       deallocate(NdIDs)
    end if

    if (allocated(the_data%NdCoords)) then
       deallocate(the_data%NdCoords)
       allocate(the_data%NdCoords(dim1*NumNd_nog))
       the_data%NdCoords = NdCoords
       deallocate(NdCoords)
    end if

    if (allocated(the_data%NdOwners)) then
       deallocate(the_data%NdOwners)
       allocate(the_data%NdOwners(NumNd_nog))
       the_data%NdOwners = NdOwners
       deallocate(NdOwners)
    end if

    if (allocated(the_data%bathymetry)) then
       deallocate(the_data%bathymetry)
       allocate(the_data%bathymetry(NumNd_nog))
       the_data%bathymetry = bathymetry
       deallocate(bathymetry)
    end if

    if (dbug) then
       ! debug, output local, global ids and coordinates of nodes accociated
       ! with element
       do i = 1, the_data%NumEl
          write(msgString,'(4I8,2F8.3)') the_data%ElIDs(i), &
            the_data%ElConnect((i-1)*NumND_per_El+1), the_data%ElConnect((i-1)*NumND_per_El+2), &
            the_data%ElConnect((i-1)*NumND_per_El+3), &
            the_data%ElCoords((i-1)*dim1+1), the_data%ElCoords((i-1)*dim1+2)
          call ESMF_LogWrite(subname//' ELEM: '//trim(msgString), ESMF_LOGMSG_INFO)
       end do

       ! debug, output nodeid, owners and coordinates
       do i = 1, the_data%NumNd
          write(msgString,'(2I8,3F10.3)') the_data%NdIDs(i), the_data%NdOwners(i), &
            the_data%NdCoords((i-1)*dim1+1), the_data%NdCoords((i-1)*dim1+2), &
            the_data%bathymetry(i)
          call ESMF_LogWrite(subname//' NODE: '//trim(msgString), ESMF_LOGMSG_INFO)
       end do
    end if

    ! message for exiting call
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine eliminate_ghosts

end module
