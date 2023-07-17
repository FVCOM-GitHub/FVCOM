!========================================================================
!PostFVCOM --- ncnest_reorder
!  reorder the element-variables in the nesting file, including:
!      UA --- depth-averaged u-velocity
!      VA --- depth-averaged v-velocity
!      U  --- u-velocity
!      V  --- v-velocity
!
! To compile:
!   ifort postFVCOM_ncnest_reorder.F90                           \
!      -L${nc_path}/lib -lnetcdff -lnetcdf  -I${nc_path}/include \
!      -o postFVCOM_ncnest_reorder
!
! To run:
!   ---Write a new output:
!      ./postFVCOM_ncnest_reorder INPUT OUTPUT
!   ---Overwrite the existing input file:
!      ./postFVCOM_ncnest_reorder INPUT
!
!
! Siqi Li, SMAST
! 2023-06-23
!========================================================================
PROGRAM postFVCOM_ncnest_reorder
  !
  USE NETCDF
  !
  IMPLICIT NONE
  !

  ! Read the command line
  character(len=200)   :: fin, fout
  character(len=500)   :: cmd
  integer              :: numarg
  LOGICAL              :: OVERWRITE
  ! Data
  INTEGER              :: nele, kbm1, nt
  INTEGER, ALLOCATABLE :: IDC0(:,:)
  REAL, ALLOCATABLE    :: UA0(:,:), VA0(:,:), U0(:,:,:), V0(:,:,:)
  INTEGER, ALLOCATABLE :: IDC(:,:), IDC1(:), k(:)
  REAL, ALLOCATABLE    :: UA(:,:), VA(:,:), U(:,:,:), V(:,:,:)
  ! NetCDF variables
  INTEGER              :: ncid, dimid, varid
  ! Index
  INTEGER              :: i, j, it

  ! Read the command line
  numarg = iargc()
  IF (numarg == 1) THEN
    CALL GETARG(1, fin)
    fout = fin
    OVERWRITE = .TRUE.
  ELSEIF (numarg == 2) THEN
    CALL GETARG(1, fin)
    CALL GETARG(2, fout)
    OVERWRITE = .FALSE.
    cmd = 'cp '//TRIM(fin)//' '//TRIM(fout)
  ELSE
    WRITE(*,*) 'Usage:'
    WRITE(*,*) './postFVCOM_ncnest_reorder INPUT OUTPUT'
    WRITE(*,*) '      OR'
    WRITE(*,*) './postFVCOM_ncnest_reorder INPUT'
    STOP
  END IF
  
  !----------------------- Read the input file -----------------------
  ! Open the input file
  CALL check( nf90_open(fin, nf90_nowrite, ncid) )
  ! Read dimension lengths
  call check( nf90_inq_dimid(ncid, 'nele', dimid) )
  call check( nf90_inquire_dimension(ncid, dimid, len=nele) )
  call check( nf90_inq_dimid(ncid, 'siglay', dimid) )
  call check( nf90_inquire_dimension(ncid, dimid, len=kbm1) )
  call check( nf90_inq_dimid(ncid, 'time', dimid) )
  call check( nf90_inquire_dimension(ncid, dimid, len=nt) )
  ! Allocate variables
  ALLOCATE(IDC0(nele, nt))
  ALLOCATE(UA0(nele, nt))
  ALLOCATE(VA0(nele, nt))
  ALLOCATE(U0(nele, kbm1, nt))
  ALLOCATE(V0(nele, kbm1, nt))
  ! Read the data
  ! IDC
  call check( nf90_inq_varid(ncid, 'idc', varid) )
  call check( nf90_get_var(ncid, varid, IDC0) )
  ! UA
  call check( nf90_inq_varid(ncid, 'ua', varid) )
  call check( nf90_get_var(ncid, varid, UA0) )
  ! VA
  call check( nf90_inq_varid(ncid, 'va', varid) )
  call check( nf90_get_var(ncid, varid, VA0) )
  ! U
  call check( nf90_inq_varid(ncid, 'u', varid) )
  call check( nf90_get_var(ncid, varid, U0) )
  ! V
  call check( nf90_inq_varid(ncid, 'v', varid) )
  call check( nf90_get_var(ncid, varid, V0) )
  ! Close the file
  call check( nf90_close(ncid) )


  !----------------------- Copy the input to output-------------------
  IF (.not. OVERWRITE) THEN
    WRITE(*,*) 'Copy the input file to a new one:'
    WRITE(*,*) '  --INPUT : '//TRIM(fin)
    WRITE(*,*) '  --OUTPUT: '//TRIM(fout)
    call SYSTEM(cmd)
  ELSE
    WRITE(*,*) 'OVERWRITE the reordered data to input file:'
    WRITE(*,*) '  --INPUT : '//TRIM(fin)
    WRITE(*,*) '  --OUTPUT: '//TRIM(fout)
  END IF

  !----------------------- Reorder the elements ----------------------
  ! Allocate variables
  ALLOCATE(IDC1(nele))
  ALLOCATE(IDC(nele, nt))
  ALLOCATE(k(nele))
  ALLOCATE(UA(nele, nt))
  ALLOCATE(VA(nele, nt))
  ALLOCATE(U(nele, kbm1, nt))
  ALLOCATE(V(nele, kbm1, nt))
  IDC1 = IDC0(:, 1)
  IDC = IDC0
  UA  = UA0
  VA  = VA0
  U   = U0
  V   = V0
  DO it = 2, nt
    IF (SUM( (IDC0(:,it)-IDC1)**2 ) == 0) THEN
      WRITE(*,'(A,I7.7,A)') '-- Time ', it, ': Unchanged'
    ELSE
      k = 0
      DO j = 1, nele
        DO i = 1, nele
          IF (IDC0(i,it) == IDC1(j)) THEN
            k(j) = i
            EXIT
          END IF
        END DO
      END DO
      IDC(:,it) = IDC1
      UA(:,it)  = UA0(k,it)
      VA(:,it)  = VA0(k,it)
      U(:,:,it) = U(k,:,it)
      V(:,:,it) = V(k,:,it)
      WRITE(*,'(A,I7.7,A)') '-- Time ', it, ': Reordered'
    END IF
  END DO


  !----------------------- Write the output file ---------------------- 
  ! Open the output file
  call check( nf90_open(fout, nf90_write, ncid) )
  ! Write the data
  ! IDC
  call check( nf90_inq_varid(ncid, 'idc', varid) )
  call check( nf90_put_var(ncid, varid, IDC) )
  ! UA
  call check( nf90_inq_varid(ncid, 'ua', varid) )
  call check( nf90_put_var(ncid, varid, UA) )
  ! VA
  call check( nf90_inq_varid(ncid, 'va', varid) )
  call check( nf90_put_var(ncid, varid, VA) )
  ! U
  call check( nf90_inq_varid(ncid, 'u', varid) )
  call check( nf90_put_var(ncid, varid, U) )
  ! V
  call check( nf90_inq_varid(ncid, 'v', varid) )
  call check( nf90_put_var(ncid, varid, V) )
  ! Close the output file
  call check( nf90_close(ncid) )        


  ! DEALLOCATE variables
  DEALLOCATE(IDC0)
  DEALLOCATE(UA0)
  DEALLOCATE(VA0)
  DEALLOCATE(U0)
  DEALLOCATE(V0)
  DEALLOCATE(IDC1)
  DEALLOCATE(IDC)
  DEALLOCATE(k)
  DEALLOCATE(UA)
  DEALLOCATE(VA)
  DEALLOCATE(U)
  DEALLOCATE(V)

CONTAINS
  SUBROUTINE check(status)
    INTEGER, INTENT(in) :: status

    IF (status /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(status))
      STOP "STOPPED"
    END IF
  END SUBROUTINE check

END PROGRAM postFVCOM_ncnest_reorder
