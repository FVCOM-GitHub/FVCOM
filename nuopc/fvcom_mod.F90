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
  use ALL_VARS,   only : NV,VX,VY,H,LAT,LON
  use MOD_PAR,    only : NGID_X,EGID_X,NLID,NMAP,ACOLLECT
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
    !! is stored in location <tt> 2*(n-1)+j</tt> of this array.
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

  logical :: NUOPC4WAV

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
    integer                      :: rc

    out_esmf_mesh=ESMF_MeshCreate(parametricDim=dim1, spatialDim=spacedim, &
       nodeIDs=the_data%NdIDs, nodeCoords=the_data%NdCoords, &
       nodeOwners=the_data%NdOwners, elementIDs=the_data%ElIDs, &
       elementTypes=the_data%ElTypes, elementConn=the_data%ElConnect, &
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
	
    print*, the_data%NumEl,NT,the_data%NumNd,MT,the_data%NumOwnedND,M,'NT,MT,M'
    allocate(the_data%NdIDs(the_data%NumNd))
    allocate(local_node_numbers(the_data%NumNd))
    allocate(the_data%ElIDs(the_data%NumEl))
    allocate(local_elem_numbers(the_data%NumEl))
    allocate(the_data%NdCoords(dim1*the_data%NumNd))
    allocate(the_data%bathymetry(the_data%NumNd))
    allocate(the_data%ElConnect(NumND_per_El*the_data%NumEl))
    allocate(the_data%NdOwners(the_data%NumNd))
    allocate(the_data%ElTypes(the_data%NumEl))

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
    the_data%ElIDs = abs(local_elem_numbers)
    do i1 = 1,the_data%NumNd
      local_node_numbers(i1) = NGID_X(i1)
    end do
    the_data%NdIDs = abs(local_node_numbers)

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
!      the_data%NdCoords((i1-1)*dim1 + 1) = lon(i1)
!      the_data%NdCoords((i1-1)*dim1 + 2) = lat(i1)
      the_data%NdCoords((i1-1)*dim1 + 1) = lat(i1)
      the_data%NdCoords((i1-1)*dim1 + 2) = lon(i1)
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
    print*,"Match ...... ",j2,the_data%NumOwnedND,myid

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

end module
