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

MODULE MOD_MLD_RHO
  !
  USE MOD_PREC
  USE ALL_VARS, ONLY : GAMMA_MIN, MLD_DEFAULT, DEEPWATER_DEPTH, DEEPWATER_GAMMA ! Siqi Li, @20210809
  !
  IMPLICIT NONE
  !
contains
  !
!---> Siqi Li, @20210809
!  SUBROUTINE MLD_RHO_CAL(RHO_AVG, MLD_RHO)
  SUBROUTINE MLD_RHO_CAL(M, KBM1, RHO_AVG, MLD_RHO)
    USE ALL_VARS, ONLY : ZZ, D
!<--- 
    IMPLICIT NONE
    !
    INTEGER    :: I, K
    INTEGER, PARAMETER    :: STD_NZ=201
!---> Siqi Li, @20210809
    INTEGER, INTENT(IN)                     :: M, KBM1
    REAL(SP), DIMENSION(M,KBM1), INTENT(IN) :: RHO_AVG !DENSITY
    REAL(SP), DIMENSION(M), INTENT(OUT)     :: MLD_RHO !MIXED LAYER DEPTH  
    REAL(SP), DIMENSION(M,KBM1)        :: DEP_SIG  !DEPTH OF EACH SIGMA LAYER
    REAL(SP), DIMENSION(M)     :: GAMMA
    REAL(SP), DIMENSION(M,STD_NZ)   :: STD_DEPTH !STANDARD DEPTH OF EACH NODE
    REAL(SP), DIMENSION(M,STD_NZ)   :: STD_RHO 

!    REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: DEP_SIG  !DEPTH OF EACH SIGMA LAYER
!    REAL(SP), ALLOCATABLE, DIMENSION(:)     :: GAMMA
!    REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: STD_DEPTH !STANDARD DEPTH OF EACH NODE
!    REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: STD_RHO 
!    REAL(SP), ALLOCATABLE, DIMENSION(:,:), INTENT(IN)   :: RHO_AVG !DENSITY
!    REAL(SP), ALLOCATABLE, DIMENSION(:), INTENT(OUT)  :: MLD_RHO !MIXED LAYER DEPTH  
!<---
    !
    ! Calculate the real depth of each sigma layer
!    ALLOCATE(DEP_SIG(M,KBM1)) ! Siqi Li, @20210809
    DO I = 1, M
      DO K = 1, KBM1
        DEP_SIG(I,K)=-ZZ(I,K)*D(I)
      END DO
    END DO
    !
!    ALLOCATE(GAMMA(M))   ! Siqi Li, @20210809
!    ALLOCATE(RHO_AVG(M,KBM1)) 
    ! Calculate the gamma
    CALL CALC_MAX_K(M, KBM1, DEP_SIG, RHO_AVG, GAMMA)
    !
!    ALLOCATE(STD_DEPTH(M,STD_NZ)) ! Siqi Li, @20210809
!    ALLOCATE(STD_RHO(M,STD_NZ))  ! Siqi Li, @20210809
    ! Interolate the rho onto standard depth.
    CALL VERTICAL_INTERP_RHO(M, KBM1, DEP_SIG, RHO_AVG, STD_NZ, STD_DEPTH, STD_RHO)
    !
!    ALLOCATE(MLD_RHO(M))  ! Siqi Li, @20210809
    ! Calculate the mixed layer depth.
    CALL CALC_MLD_DENSITY(M, STD_NZ, STD_DEPTH, STD_RHO, GAMMA, MLD_RHO)
    !
  END SUBROUTINE MLD_RHO_CAL
  !
  subroutine calc_max_k(node, nz, depth, rho, gamma)
    !
    integer, intent(in)      :: node, nz
    ! ------- new: changes 2021 April K. Lettmann ---------
    !real, dimension(node,nz), intent(in) :: depth, rho  ! original line
    !real, dimension(node), intent(out)   :: gamma       ! original line
    REAL(SP), dimension(node,nz), intent(in) :: depth, rho 
    REAL(SP), dimension(node), intent(out)   :: gamma     
    ! -----------------------------------------------------
    !
    real(sp), dimension(nz-1)  :: k
    integer                :: i, j
    !
    do i=1, node
      do j=1, nz-1
        k(j)=(rho(i,j+1)-rho(i,j))/(depth(i,j+1)-depth(i,j))
      end do
      gamma(i)=maxval(k)
    ! Lu Wang modified for deeper region (depth is deeper than 50m) @May22,2019
!---> Siqi Li, 20210809
!      if (depth(i,nz) >100.0) then  !lwang at Jul. 26
!         gamma(i)=0.03e-3
      if (depth(i,nz) >DEEPWATER_DEPTH) then  !lwang at Jul. 26
         gamma(i)=DEEPWATER_GAMMA
!<---
      end if
    ! Lu Wang 
    end do
       
    !
  end subroutine calc_max_k
  !
  subroutine v_interp(n1,h1,var1,n2,h2,var2)
    USE ALL_VARS, only: KBM1
    !
    integer, intent(in)       :: n1, n2
    ! ------- new: changes 2021 April K. Lettmann ---------
    !real, dimension(n1), intent(in) :: h1(n1), var1(n1)  ! original line
    !real, dimension(n2), intent(in) :: h2(n2)            ! original line
    !real, dimension(n2), intent(out):: var2(n2)          ! original line
    REAL(SP), dimension(n1), intent(in) :: h1(n1), var1(n1)  
    REAL(SP), dimension(n2), intent(in) :: h2(n2)            
    REAL(SP), dimension(n2), intent(out):: var2(n2)            
    ! -----------------------------------------------------
    !
    integer :: i, j
    real(sp), parameter          :: missing_value=-99.
    !
    do j=1,n2
      if (h2(j)<h1(1)) then      ! If shallower than the first sigma layer
        var2(j)=var1(1)
      elseif (h2(j)>h1(KBM1)) then ! If deeper than the last sigma layer
        var2(j)=missing_value
      else
        do i=2,n1
          if (h2(j)>=h1(i-1) .and. h2(j)<=h1(i)) then
            var2(j)=var1(i-1)*(h1(i)-h2(j))/(h1(i)-h1(i-1))+var1(i)*(h2(j)-h1(i-1))/(h1(i)-h1(i-1))
          end if
        end do
      end if
    end do       
    !
  end subroutine v_interp
  !
  subroutine vertical_interp_rho(node, nz, depth, rho, std_nz, std_depth, std_rho)
    !
    integer, intent(in)      :: node, nz
    ! ------- new: changes 2021 April K. Lettmann ---------
    !real, dimension(node,nz), intent(in)  :: depth, rho ! original line
    REAL(SP), dimension(node,nz), intent(in)  :: depth, rho
    ! -----------------------------------------------------
    integer, intent(in)      :: std_nz
    ! ------- new: changes 2021 April K. Lettmann ---------
    !real, dimension(node,std_nz), intent(out) :: std_depth, std_rho ! original line
    REAL(SP), dimension(node,std_nz), intent(out) :: std_depth, std_rho
    ! -----------------------------------------------------
    !
    integer                  :: i,j
    real(sp), parameter          :: missing_value=-99.
    !
    do i=1,node
      if (depth(i,nz)<nz) then

        std_depth(i,1)=0.
        std_depth(i,2:nz+1)=depth(i,1:nz)
        std_depth(i,nz+2:std_nz)=missing_value

        std_rho(i,1)=rho(i,1)
        std_rho(i,2:nz+1)=rho(i,1:nz)
        std_rho(i,nz+2:std_nz)=missing_value
      else 
        do j=1,std_nz
          std_depth(i,j)=dble(j-1)
        end do
        call v_interp(nz,depth(i,:),rho(i,:),std_nz,std_depth(i,:),std_rho(i,:))
      end if
    end do
    !
  end subroutine vertical_interp_rho 
  !
  subroutine calc_mld_density(node, nz, depth, rho, gamma, mld)
    !
    integer, intent(in)             :: node, nz
    ! ------- new: changes 2021 April K. Lettmann ---------
    !real, dimension(node), intent(in)     :: gamma       ! original line
    !real, dimension(node, nz), intent(in) :: depth, rho  ! original line
    !real, dimension(node), intent(out)    :: mld         ! original line
    REAL(SP), dimension(node), intent(in)     :: gamma       
    REAL(SP), dimension(node, nz), intent(in) :: depth, rho  
    REAL(SP), dimension(node), intent(out)    :: mld         
    ! -----------------------------------------------------
    !
!    real, parameter   :: gamma_min=0.04e-3 !0.005 Lu may22,2019
    integer           :: iz
    real(sp)              :: h, vdif
    integer           :: i, j
    !
    do i=1, node
!    do i=119778,119778
      !
      ! Get the last layer of water for this node point.
      iz=nz
      do j=1,nz
        if (rho(i,j)<-90.) then
          iz=j-1
          exit
        end if
      end do
      !
      ! Judge the gamma
      ! Lu added the depth limit may22,2019
!      if (gamma(i)<gamma_min .and. depth(i,iz)<100.0) then ! lwang Jul. 26
      if (gamma(i)<gamma_min .and. depth(i,iz)<DEEPWATER_DEPTH) then ! Siqi Li, @20210809
        mld(i)=depth(i,iz)
      else 
        ! a.
        h=0. 
        do j=2,iz
          h=h+(rho(i,j-1)+rho(i,j))/2*(depth(i,j)-depth(i,j-1))
        end do
        ! b.
        vdif=h-rho(i,1)*(depth(i,iz)-depth(i,1))
! Lu modified in order @ Jun 28, 2019
        if (vdif<0.) then
!!          write(*,*) '===================================================='
!!          write(*,*) 'Warning: vdif is negative in MLD calculation!'
!!          write(*,*) 'NODE ID: ',i
!!          write(*,*) gamma(i),gamma_min
!!          write(*,'(F10.2,F20.8)') (depth(i,j),rho(i,j),j=1,iz)
!          mld(i)=5.0 
          mld(i)=min(MLD_DEFAULT, depth(i,iz))   ! Siqi Li, @20210809
        else
          mld(i)=depth(i,iz)-sqrt(2.*vdif/gamma(i)) 
        end if
        if (mld(i)<0.) then
!           mld(i)=5.0
           mld(i)=min(MLD_DEFAULT, depth(i,iz)) ! Siqi Li, @20210809
        end if
! Lu
      end if
    end do
    !
  end subroutine calc_mld_density
  !
END MODULE MOD_MLD_RHO
