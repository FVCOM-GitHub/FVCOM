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

!==================================================== John C. Warner ===
!                                                                      !
!  This routine computes floc transformation.                          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Verney, R., Lafite, R., Claude Brun-Cottan, J., & Le Hir, P. (2011).!
!    Behaviour of a floc population during a tidal cycle: laboratory   !
!    experiments and numerical modelling. Continental Shelf Research,  !
!    31(10), S64-S83.                                                  !
!=======================================================================

Module Mod_FlocMod
#if defined (SEDIMENT) && (CSTMS_SED)

# if defined (WAVE_CURRENT_INTERACTION)
USE MOD_WAVE_CURRENT_INTERACTION
# endif

Use All_vars
Use Mod_Par
Use Mod_Prec 
Use Mod_Types


Use Mod_CSTMS_vars

! 
!  switched to activate the flocmod model
!
  implicit none
  
  public

!  l_ADS           : boolean set to .true. if differential settling aggregation 
!
!  l_ASH           : boolean set to .true. if shear aggregation
!
!  l_COLLFRAG      : boolean set to .true. if collision-induced fragmentation enable
!  f_dp0           : Primary particle size (m), typically 4e-6 m
!
!  f_nf            : Floc fractal dimension, typically ranging from 1.6 to 2.6
!
!  f_dmax          : (unused) Maximum diameter (m)
!
!  f_nb_frag       : number of fragments by shear erosion. If binary/ternary : 2.
!
!  f_alpha         :  flocculation efficiency, ranging from 0 .to 1.
!
!  f_beta          : shear fragmentation rate
!
!  f_ater          : for ternary breakup, use 0.5, for binary : 0.
!                   (a boolean could be better)
!
!  f_ero_frac      : fraction of the shear fragmentation term 
!                    transfered to shear erosion. 
!                    Ranging from 0. (no erosion) to 1. (all erosion)
!
!  f_ero_nbfrag    : Number of fragments induced by shear erosion.   
!
!  f_ero_iv        : fragment size class (could be changed to a particle 
!                    size or a particle distribution (INTEGER)
!
!  f_collfragparam : Fragmentation rate for collision-induced breakup
!
! f_clim           : min concentration below which flocculation
!                    processes are not calculated
!
! l_testcase - if .TRUE. sets G(t) to values from Verney et al., 2011 lab experiment
!
! GLS_P    : Stability exponent (non-dimensional).
!
! GLS_M    : Turbulent kinetic energy exponent (non-dimensional).
!
! GLS_N    : Turbulent length scale exponent (non-dimensional).
!
! GLS_CMU0 : Stability coefficient.


  logical :: l_ASH
  logical :: l_ADS
  logical :: l_COLLFRAG
  logical :: l_testcase
  INTEGER  :: f_ero_iv
  real(sp) :: f_dp0,f_alpha,f_beta,f_nb_frag
  real(sp) :: f_dmax,f_ater,f_clim
  real(sp) :: f_ero_frac,f_ero_nbfrag
  real(sp) :: f_nf
  real(sp) :: f_frag
  real(sp) :: f_fter
  real(sp) :: f_collfragparam
  real(sp) :: gls_p
  real(sp) :: gls_m
  real(sp) :: gls_n
  real(sp) :: gls_cmu0
  real(sp), parameter :: rhoref = 1030.0_sp

  
  ! Dynamics selecting switches
  logical :: FLOC_TURB_DISS
  logical :: FLOC_BBL_DISS


  contains


  SUBROUTINE allocate_sedflocs

!  Imported variable declarations.
!
!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
!
  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: allocate_sedflocs"

      allocate ( SEDFLOCS % f_diam(NCS) )
      allocate ( SEDFLOCS % f_vol(NCS) )
      allocate ( SEDFLOCS % f_rho(NCS) )
      allocate ( SEDFLOCS % f_cv(NCS) )
      allocate ( SEDFLOCS % f_l3(NCS) )
      allocate ( SEDFLOCS % f_mass(0:NCS+1) )
      allocate ( SEDFLOCS % f_coll_prob_sh(NCS,NCS) )
      allocate ( SEDFLOCS % f_coll_prob_ds(NCS,NCS) )
      allocate ( SEDFLOCS % f_l1_sh(NCS,NCS) )
      allocate ( SEDFLOCS % f_l1_ds(NCS,NCS) )
      allocate ( SEDFLOCS % f_g3(NCS,NCS) )
      allocate ( SEDFLOCS % f_l4(NCS,NCS) )
      allocate ( SEDFLOCS % f_g1_sh(NCS,NCS,NCS) )
      allocate ( SEDFLOCS % f_g1_ds(NCS,NCS,NCS) )
      allocate ( SEDFLOCS % f_g4(NCS,NCS,NCS) )
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: allocate_sedflocs"

  RETURN

  END SUBROUTINE allocate_sedflocs


  SUBROUTINE initialize_sedflocs_param (f_mass,f_diam,f_g1_sh,          &
     &                              f_g1_ds,f_g3,f_l1_sh,f_l1_ds,       &
     &                              f_coll_prob_sh,f_coll_prob_ds,      &
     &                              f_l3)


  Use Scalar
  Use Control, only : ireport,iint,msr,par,pi
  Use Lims,    only : m,mt,nt,kbm1,numqbc,kb,nprocs,myid

  implicit none 
!
!  Imported variable declarations.
!
      real(sp), intent(inout) :: f_mass(0:NCS+1)
      real(sp), intent(inout) :: f_diam(NCS)
      real(sp), intent(inout) :: f_g1_sh(NCS,NCS,NCS)
      real(sp), intent(inout) :: f_g1_ds(NCS,NCS,NCS)
      real(sp), intent(inout) :: f_g3(NCS,NCS)
      real(sp), intent(inout) :: f_l1_sh(NCS,NCS)
      real(sp), intent(inout) :: f_l1_ds(NCS,NCS)
      real(sp), intent(inout) :: f_coll_prob_sh(NCS,NCS)
      real(sp), intent(inout) :: f_coll_prob_ds(NCS,NCS)
      real(sp), intent(inout) :: f_l3(NCS)
!
!  Local variable declarations.
!
   logical  :: f_test
   real(sp) :: f_weight,mult,dfragmax
   integer  :: iv1,iv2,iv3,iv,itrc
   real(sp) :: f_vol(NCS),f_rho(NCS)
   real(sp), parameter :: mu = 0.001_sp
   real(sp) :: eps 
   REAL(DP), PARAMETER :: g = 9.81_SP
   eps = epsilon(1.0)
!
! ALA the rest of the initialization
!
!      l_ADS=.false.
!      l_ASH=.true.
!      l_COLLFRAG=.false.
!      f_dp0=0.000004_sp
!      f_nf=1.9_sp
!      f_dmax=0.001500_sp
!      f_nb_frag=2.0_sp
!      f_alpha=0.35_sp
!      f_beta=0.15_sp
!      f_ater=0.0_sp
!      f_ero_frac=0.0_sp
!      f_ero_nbfrag=2.0_sp
!      f_ero_iv=1
!      f_collfragparam=0.01_sp

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: initialize_sedflocs_param"


!!--------------------------------------------------
!! floc characteristics
   DO itrc=1,NCS
      f_diam(itrc)=sed(itrc)%Sd50
      f_vol(itrc)=pi/6.0_sp*(f_diam(itrc))**3.0_sp
      f_rho(itrc)=rhoref+(2650.0_sp-rhoref)*                          &
    &     (f_dp0/f_diam(itrc))**(3.0_sp-f_nf)
      f_mass(itrc)=f_vol(itrc)*(f_rho(itrc)-rhoref)
   ENDDO
   f_mass(NCS+1)=f_mass(NCS)*2.0_sp+1.0_sp  
   IF (f_diam(1).eq.f_dp0)  THEN
       f_mass(1)=f_vol(1)*sed(1)%Srho
   ENDIF

   IF(MSR)THEN
!TODO - This wont parallelize 
      WRITE(*,*) ' '
      WRITE(*,*) '*** FLOCMOD INIT *** '
      write(*,*) 'NCS, NNS:',NCS,NNS   
      WRITE(*,*) 'class diameter (um)  volume (m3)  density (kg/m3)  mass (kg)'
      DO itrc=1,NCS
         WRITE(*,*) itrc,f_diam(itrc)*1e6,f_vol(itrc),f_rho(itrc),f_mass(itrc)
      ENDDO
      write(*,*) 'f_mass(0) and f_mass(NCS+1): ',f_mass(0),f_mass(NCS+1)
      WRITE(*,*) ' '
      WRITE(*,*) ' *** PARAMETERS ***'
      WRITE(*,*) 'Primary particle size (f_dp0)                                : ',f_dp0
      WRITE(*,*) 'Fractal dimension (f_nf)                                     : ',f_nf
      WRITE(*,*) 'Flocculation efficiency (f_alpha)                            : ',f_alpha
      WRITE(*,*) 'Floc break up parameter (f_beta)                             : ',f_beta
      WRITE(*,*) 'Nb of fragments (f_nb_frag)                                  : ',f_nb_frag
      WRITE(*,*) 'Ternary fragmentation (f_ater)                               : ',f_ater
      WRITE(*,*) 'Floc erosion (% of mass) (f_ero_frac)                        : ',f_ero_frac
      WRITE(*,*) 'Nb of fragments by erosion (f_ero_nbfrag)                    : ',f_ero_nbfrag
      WRITE(*,*) 'fragment class (f_ero_iv)                                    : ',f_ero_iv
!      WRITE(*,*) 'negative mass tolerated before redistribution (f_mneg_param) : ',f_mneg_param
      WRITE(*,*) 'Boolean for differential settling aggregation (L_ADS)        : ',l_ADS
      WRITE(*,*) 'Boolean for shear aggregation (L_ASH)                        : ',l_ASH
      WRITE(*,*) 'Boolean for collision fragmenation (L_COLLFRAG)              : ',l_COLLFRAG
      WRITE(*,*) 'Collision fragmentation parameter (f_collfragparam)          : ',f_collfragparam
      WRITE(*,*) ' '
      WRITE(*,*) 'Value of eps                                                 : ',eps
      WRITE(*,*) ' '
   END IF
! kernels computation former SUBROUTINE flocmod_kernels


!!--------------------------------------------------------------------------
!! * Executable part
      f_test=.true.
      dfragmax=0.00003_sp
! compute collision probability former SUBROUTINE flocmod_agregation_statistics

      DO iv1=1,NCS
        DO iv2=1,NCS

          f_coll_prob_sh(iv1,iv2)=1.0_sp/6.0_sp*(f_diam(iv1)+           &
     &            f_diam(iv2))**3.0_sp

          f_coll_prob_ds(iv1,iv2)=0.25_sp*pi*(f_diam(iv1)+              &
     &            f_diam(iv2))**2.0_sp*g/mu*abs((f_rho(iv1)-         &
     &            rhoref)*f_diam(iv1)**2.0_sp-(f_rho(iv2)-rhoref)*      &
     &            f_diam(iv2)**2.0_sp)

        ENDDO
      ENDDO

  !********************************************************************************
  ! agregation : GAIN : f_g1_sh and f_g1_ds
  !********************************************************************************
      DO iv1=1,NCS
       DO iv2=1,NCS
        DO iv3=iv2,NCS
           IF((f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1-1)    .and.     &
     &            ((f_mass(iv2)+f_mass(iv3)) .le. f_mass(iv1))) THEN

              f_weight=(f_mass(iv2)+f_mass(iv3)-f_mass(iv1-1))/         &
     &            (f_mass(iv1)-f_mass(iv1-1))

           ELSEIF ((f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1) .and.     &
     &            ((f_mass(iv2)+f_mass(iv3)) .lt. f_mass(iv1+1))) THEN

              IF (iv1 .eq. NCS) THEN
                 f_weight=1.0_sp
              ELSE
                 f_weight=1.0_sp-(f_mass(iv2)+f_mass(iv3)-              &
     &                 f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))
              ENDIF

           ELSE
              f_weight=0.0_sp
           ENDIF

           f_g1_sh(iv2,iv3,iv1)=f_weight*f_alpha*                       &
     &            f_coll_prob_sh(iv2,iv3)*(f_mass(iv2)+                 &
     &            f_mass(iv3))/f_mass(iv1)
           f_g1_ds(iv2,iv3,iv1)=f_weight*f_alpha*                       &
     &            f_coll_prob_ds(iv2,iv3)*(f_mass(iv2)+                 &
     &            f_mass(iv3))/f_mass(iv1)

        ENDDO
       ENDDO
      ENDDO

  !********************************************************************************
  ! Shear fragmentation : GAIN : f_g3
  !********************************************************************************

      DO iv1=1,NCS
       DO iv2=iv1,NCS

        IF (f_diam(iv2).gt.dfragmax) THEN
           ! binary fragmentation

           IF (f_mass(iv2)/f_nb_frag .gt. f_mass(iv1-1) &
                .and. f_mass(iv2)/f_nb_frag .le. f_mass(iv1)) THEN

              IF (iv1 .eq. 1) THEN 
                 f_weight=1.0_sp
              ELSE
                 f_weight=(f_mass(iv2)/f_nb_frag-f_mass(iv1-1))/        &
     &                (f_mass(iv1)-f_mass(iv1-1))
              ENDIF

           ELSEIF (f_mass(iv2)/f_nb_frag .gt. f_mass(iv1)               &
     &            .and. f_mass(iv2)/f_nb_frag .lt. f_mass(iv1+1)) THEN

              f_weight=1.0_sp-(f_mass(iv2)/f_nb_frag-f_mass(iv1))/      &
     &              (f_mass(iv1+1)-f_mass(iv1))

           ELSE

              f_weight=0.0_sp

           ENDIF

        ELSE
           f_weight=0.0_sp
        ENDIF

        f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0_sp-f_ero_frac)*(1.0_sp-        &
     &            f_ater)*f_weight*f_beta*f_diam(iv2)*((f_diam(iv2)-    &
     &            f_dp0)/f_dp0)**(3.0_sp-f_nf)*f_mass(iv2)/             &
     &            f_mass(iv1)

        ! ternary fragmentation
        IF (f_diam(iv2).gt.dfragmax) THEN
           IF (f_mass(iv2)/(2.0_sp*f_nb_frag) .gt. f_mass(iv1-1) .and.  &
     &          f_mass(iv2)/(2.0_sp*f_nb_frag) .le. f_mass(iv1)) THEN

              IF (iv1 .eq. 1) THEN 
                 f_weight=1.0_sp
              ELSE
                 f_weight=(f_mass(iv2)/(2.0_sp*f_nb_frag)-              &
     &                f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
              ENDIF

           ELSEIF (f_mass(iv2)/(2.0_sp*f_nb_frag) .gt. f_mass(iv1)      &
     &           .and. f_mass(iv2)/(2.0_sp*f_nb_frag) .lt.              &
     &            f_mass(iv1+1)) THEN

              f_weight=1.0_sp-(f_mass(iv2)/(2.0_sp*f_nb_frag)-          &
     &                 f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))

           ELSE
              f_weight=0.0_sp

           ENDIF
           ! update for ternary fragments
           f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0_sp-f_ero_frac)*(f_ater)*    &
     &            f_weight*f_beta*f_diam(iv2)*((f_diam(iv2)-f_dp0)/     &
     &            f_dp0)**(3.0_sp-f_nf)*f_mass(iv2)/f_mass(iv1)   

           ! Floc erosion

           IF ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt.         &
     &            f_mass(f_ero_iv)) THEN

              IF (((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt.     &
     &            f_mass(iv1-1)) .and. (f_mass(iv2)-f_mass(f_ero_iv)*   &
     &            f_ero_nbfrag) .le. f_mass(iv1)) THEN

                 IF (iv1 .eq. 1) THEN
                    f_weight=1.0_sp
                 ELSE
                    f_weight=(f_mass(iv2)-f_mass(f_ero_iv)*             &
     &                      f_ero_nbfrag-f_mass(iv1-1))/(f_mass(iv1)-   &
     &                      f_mass(iv1-1))
                 ENDIF

              ELSEIF ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt.  &
     &            f_mass(iv1) .and. (f_mass(iv2)-f_mass(f_ero_iv)*      &
     &            f_ero_nbfrag) .lt. f_mass(iv1+1)) THEN

                 f_weight=1.0_sp-(f_mass(iv2)-f_mass(f_ero_iv)*         &
     &                    f_ero_nbfrag-f_mass(iv1))/(f_mass(iv1+1)-     &
     &                    f_mass(iv1))

              ELSE
                 f_weight=0.0_sp
              ENDIF

              ! update for eroded floc masses 

              f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_weight*f_beta*   &
     &            f_diam(iv2)*(max(eps,(f_diam(iv2)-f_dp0))/f_dp0)**    &
     &            (3.0_sp-f_nf)*(f_mass(iv2)-f_mass(f_ero_iv)*          &
     &            f_ero_nbfrag)/f_mass(iv1)

              IF (iv1 .eq. f_ero_iv) THEN

                 f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_beta*         &
     &              f_diam(iv2)*(max(eps,(f_diam(iv2)-f_dp0))/f_dp0)**  &
     &              (3.0_sp-f_nf)*f_ero_nbfrag*f_mass(f_ero_iv)/        &
     &              f_mass(iv1)
              ENDIF
           ENDIF
        ENDIF ! condition on dfragmax
       ENDDO
      ENDDO

  !********************************************************************************
  !  Shear agregation : LOSS : f_l1
  !********************************************************************************

      DO iv1=1,NCS
       DO iv2=1,NCS

        IF(iv2 .eq. iv1) THEN
           mult=2.0_sp
        ELSE
           mult=1.0_sp
        ENDIF

        f_l1_sh(iv2,iv1)=mult*f_alpha*f_coll_prob_sh(iv2,iv1) 
        f_l1_ds(iv2,iv1)=mult*f_alpha*f_coll_prob_ds(iv2,iv1) 

       ENDDO
      ENDDO

  !********************************************************************************
  !  Shear fragmentation : LOSS : f_l2
  !********************************************************************************

      DO iv1=1,NCS
       f_l3(iv1)=0.0_sp
       IF (f_diam(iv1).gt.dfragmax) THEN
        ! shear fragmentation
           f_l3(iv1)=f_l3(iv1)+(1.0_sp-f_ero_frac)*f_beta*f_diam(iv1)*  &
     &            ((f_diam(iv1)-f_dp0)/f_dp0)**(3.0_sp-f_nf)

        ! shear erosion
        IF ((f_mass(iv1)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt.            &
     &            f_mass(f_ero_iv)) THEN
           f_l3(iv1)=f_l3(iv1)+f_ero_frac*f_beta*f_diam(iv1)*           &
     &            ((f_diam(iv1)-f_dp0)/f_dp0)**(3.0_sp-f_nf)
        ENDIF
       ENDIF
      ENDDO

   IF(MSR)THEN
      WRITE(*,*) ' '
      write(*,*) 'Sum of kernal coefficients:'
      write(*,*) 'f_coll_prob_sh',sum(f_coll_prob_sh)
      write(*,*) 'f_coll_prob_ds',sum(f_coll_prob_ds)
      write(*,*) 'f_g1_sh',sum(f_g1_sh)
      write(*,*) 'f_g1_ds',sum(f_g1_ds)
      write(*,*) 'f_l1_sh',sum(f_l1_sh)
      write(*,*) 'f_l1_ds',sum(f_l1_ds)
      write(*,*) 'f_g3',sum(f_g3)
      write(*,*) 'f_l3',sum(f_l3)
      WRITE(*,*) ' '
      WRITE(*,*) '*** END FLOCMOD INIT *** '   
   END IF

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: initialize_sedflocs_param"

!  END FORMER SUBROUTINE flocmod_kernels

  RETURN

  END SUBROUTINE initialize_sedflocs_param




!=================================================================== ===
!                                                                      !
!  This routine computes floc transformation.                          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Verney, R., Lafite, R., Claude Brun-Cottan, J., & Le Hir, P. (2011).!
!    Behaviour of a floc population during a tidal cycle: laboratory   !
!    experiments and numerical modelling. Continental Shelf Research,  !
!    31(10), S64-S83.                                                  !
!=======================================================================
!
  Subroutine calc_sed_flocmod

  implicit none

  integer :: i, indx, ised, j, k, ks

!
! Variable declarations for floc model
!
  integer  :: iv1
  real(sp), dimension(0:kb) :: Hz_inv
  real(sp) :: Gval,diss,mneg,dttemp,f_dt
  real(sp) :: dt1,f_csum,epsilon8
  real(sp) :: cvtotmud,tke_av, gls_av, exp1, exp2, exp3, ustr2
  real(sp), dimension(0:kb,ncs) :: susmud
  real(sp), dimension(0:kb,0:mt) :: f_davg
  real(sp), dimension(0:kb,0:mt) :: f_d50
  real(sp), dimension(0:kb,0:mt) :: f_d90
  real(sp), dimension(0:kb,0:mt) :: f_d10
  real(sp), dimension(1:NCS)  :: cv_tmp,NNin,NNout
! f_mneg_param : negative mass tolerated to avoid small sub time step (g/l)
  real(sp), parameter :: f_mneg_param=0.000_sp
  real(sp), parameter :: vonKar = 0.41_sp

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: calc_sed_flocmod "


      epsilon8=epsilon(1.0)
!     epsilon8=1.e-8
!
!--------------------------------------------------------------------------
! * Executable part
!
  do i=1,m 
!
!  Extract mud variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!
!  Compute inverse thicknesses of vertical layer to avoid repeated divisions. 
     DO k=1,kbm1
        Hz_inv(k)=1.0_sp/(dt(i)*dz(i,k))
     END DO
     DO ised=1,NCS
        DO k=1,kbm1
           susmud(k,ised)=sed(ised)%conc(i,k)*Hz_inv(k)
         ENDDO
     ENDDO
 
! min concentration below which flocculation processes are not calculated
!      f_clim=0.001_sp 
     exp1 = 3.0_sp+gls_p/gls_n
     exp2 = 1.5_sp+gls_m/gls_n
     exp3 = -1.0_sp/gls_n

     DO k=kbm1,1,-1
 
        f_dt=DTsed
        dttemp=0.0_sp

        ! concentration of all mud classes in one grid cell
        cvtotmud=0.0_sp
        DO ised=1,NCS
           cv_tmp(ised)=susmud(k,ised) 
           cvtotmud=cvtotmud+cv_tmp(ised)

           NNin(ised)=cv_tmp(ised)/sedflocs%f_mass(ised)
        ENDDO

        DO iv1=1,NCS
           IF (NNin(iv1).lt.0.0_sp) THEN
             WRITE(*,*) '***************************************'
             WRITE(*,*) 'CAUTION, negative mass at node i,k :',  &
     &                          i,k
             WRITE(*,*) '***************************************'        
           ENDIF
        ENDDO

        IF (cvtotmud .gt. f_clim) THEN

          if(FLOC_TURB_DISS .and. (.not.FLOC_BBL_DISS))then
!
!ALA dissipation from turbulence clossure
!

!  tke : Turbulent energy squared (m2/s2)                        !
!  gls : Turbulent energy squared times turbulent length scale   !

#  if defined (GOTM)   
            ! use gotm's dissipation rate with caution.
            diss = teps(i,k+1)
#  else
    ! The original turbulence closure scheme is applied and has been tested.
	    IF (k.eq.kbm1) THEN
               tke_av = q2(i,k)
               gls_av = q2l(i,k)

	    ELSEIF (k.eq.1) THEN
               tke_av = q2(i,k)
               gls_av = q2l(i,k)
	    ELSE 
               tke_av = q2(i,k+1)+q2(i,k)
               gls_av = q2l(i,k+1)+q2l(i,k)
            ENDIF
            !tke_av = q2(i,k)
            !gls_av = q2l(i,k)
            diss = gls_cmu0**exp1*tke_av**exp2*gls_av**exp3
#  endif
          elseif(FLOC_BBL_DISS .and. (.not.FLOC_TURB_DISS))then
!
! ALA dissipation from wavecurrent bottom stress
! NOT READY FOR PRIME TIME
! NEEDS VERTICAL DISTRIBUTION
!
              WRITE(*,*) 'CAUTION :'
              WRITE(*,*) 'CHOOSE A DISSIPATION MODEL FOR FLOCS'
              WRITE(*,*) 'BBL MODEL IS NOT AVAILABLE'
              WRITE(*,*) 'SIMULATION STOPPED'
              STOP
          else 
            diss = epsilon8
            IF (l_testcase) THEN
!             if (j.eq.1.and.i.eq.1.and.k.eq.1) then
!               WRITE(*,*) 'VERNEY ET AL TESTCASE FOR FLOCS'
!             endif
            ELSE    
              WRITE(*,*) 'CAUTION :'
              WRITE(*,*) 'CHOOSE A DISSIPATION MODEL FOR FLOCS'
              WRITE(*,*) 'SIMULATION STOPPED'
              STOP
            ENDIF                 
          endif	
      
          CALL flocmod_comp_g(k,i,Gval,diss) 

          DO WHILE (dttemp .le. DTsed)
             CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt)
             CALL flocmod_mass_control(NNout,mneg)
             IF (mneg .gt. f_mneg_param) THEN
                DO WHILE (mneg .gt. f_mneg_param)
                  f_dt=MIN(f_dt/2.0_sp,DTsed-dttemp)
                  IF (f_dt.lt.epsilon8) THEN
	             CALL flocmod_mass_redistribute(NNin)
                     dttemp=DTsed
                     exit
                  ENDIF
                  CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt)
                  CALL flocmod_mass_control(NNout,mneg)   
                ENDDO
             ELSE

                IF (f_dt.lt.DTsed) THEN
                   DO while (mneg .lt.f_mneg_param) 
                      IF (dttemp+f_dt .eq. DTsed) THEN
                         CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt)
                         exit
                      ELSE
                         dt1=f_dt
                         f_dt=MIN(2.0_sp*f_dt,DTsed-dttemp)
                         CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt)
                         CALL flocmod_mass_control(NNout,mneg)
                         IF (mneg .gt. f_mneg_param) THEN 
                           f_dt=dt1
                           CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt)
                           exit
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF
             dttemp=dttemp+f_dt

             NNin(:)=NNout(:) ! update new Floc size distribution
             ! redistribute negative masses IF any on positive classes, 
             ! depends on f_mneg_param
	     CALL flocmod_mass_redistribute(NNin)

             IF (dttemp.eq.DTsed) exit

          ENDDO ! loop on full dt

        ENDIF ! only if cvtotmud > f_clim

        999 continue
        ! update mass concentration for all mud classes
        DO ised=1,NCS
           susmud(k,ised)=NNin(ised)*sedflocs%f_mass(ised)
        ENDDO

     ENDDO
!
!-----------------------------------------------------------------------
!  Update global tracer variables.
!-----------------------------------------------------------------------
!
     DO ised=1,NCS
        indx = idsed(ised)
        DO k=1,kbm1
             sed(ised)%conc(i,k)=susmud(k,ised)*(dt(i)*dz(i,k))
        ENDDO
     ENDDO
      
      
  End do
!      WRITE(*,*) 'END flocmod_main'

  if(dbg_set(dbg_sbr)) write(ipt,*) "End:  calc_sed_flocmod"

  End Subroutine calc_sed_flocmod




!!===========================================================================
      SUBROUTINE flocmod_collfrag(Gval,f_g4,f_l4)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_collfrag  ***
  !&E
  !&E ** Purpose : computation of collision fragmentation term, 
  !&E              based on McAnally and Mehta, 2001
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_comp_fsd
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
  !
      implicit none 
 
      real(sp),intent(in) :: Gval

  !! * Local declarations
      integer      :: iv1,iv2,iv3
      real(sp) :: f_fp,f_fy,f_cfcst,gcolfragmin,gcolfragmax
      real(sp) :: gcolfragiv1,gcolfragiv2,f_weight,mult
      real(sp) :: cff1, cff2
      real(sp),DIMENSION(NCS,NCS,NCS) :: f_g4
      real(sp),DIMENSION(NCS,NCS) :: f_l4
      real(sp),DIMENSION(NCS) :: fdiam,fmass
 
      if(dbg_set(dbg_sbr)) write(ipt,*) "Start:  flocmod_collfrag"


! shorten names
      fdiam=SEDFLOCS%f_diam
      fmass=SEDFLOCS%f_mass
!  !!--------------------------------------------------------------------------
  !! * Executable part

      f_fp=0.1_sp
      f_fy=1e-10
      f_cfcst=3.0_sp/16.0_sp
      f_g4(1:NCS,1:NCS,1:NCS)=0.0_sp
      f_l4(1:NCS,1:NCS)=0.0_sp
      cff1=2.0_sp/(3.0_sp-f_nf)
      cff2=1.0_sp/rhoref

      DO iv1=1,NCS
       DO iv2=1,NCS
        DO iv3=iv2,NCS
           ! fragmentation after collision probability based on 
	   ! Gval for particles iv2 and iv3
           ! gcolfrag=(collision induced shear) / (floc strength)

           gcolfragmin=2.0_sp*(Gval*(fdiam(iv2)+fdiam(iv3)))            &
     &          **2.0_sp*fmass(iv2)*fmass(iv3)/(pi*f_fy*f_fp*           &
     &          fdiam(iv3)**2.0_sp*(fmass(iv2)+fmass(iv3))              &
     &          *((SEDFLOCS%f_rho(iv3)-rhoref)*cff2)**cff1)

           gcolfragmax=2.0_sp*(Gval*(fdiam(iv2)+fdiam(iv3)))            &
     &          **2.0_sp*fmass(iv2)*fmass(iv3)/(pi*f_fy*f_fp*           &
     &          fdiam(iv2)**2.0_sp*(fmass(iv2)+fmass(iv3))              &
     &          *((SEDFLOCS%f_rho(iv2)-rhoref)*cff2)**cff1)


           ! first case : iv3 not eroded, iv2 eroded forming 2 particles
	   ! : iv3+f_cfcst*iv2 / iv2-f_cfcst*iv2
           IF (gcolfragmin.lt.1.0_sp .and. gcolfragmax.ge.1.0_sp) THEN  

              IF (((fmass(iv3)+f_cfcst*fmass(iv2)).gt.fmass(iv1-1))     &
     &             .and. ((fmass(iv3)+f_cfcst*fmass(iv2)).le.           &
     &             fmass(iv1))) THEN

                 f_weight=((fmass(iv3)+f_cfcst*fmass(iv2)-              &
     &                     fmass(iv1-1))/(fmass(iv1)-fmass(iv1-1)))

              ELSEIF (fmass(iv3)+f_cfcst*fmass(iv2).gt.fmass(iv1)       &
     &              .and. fmass(iv3)+f_cfcst*fmass(iv2).lt.             &
     &              fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_sp
                 ELSE

                    f_weight=1.0_sp-((fmass(iv3)+f_cfcst*fmass(iv2)-    &
     &                      fmass(iv1))/(fmass(iv1+1)-fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_sp
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &            (SEDFLOCS%f_coll_prob_sh(iv2,iv3))*(fmass(iv3)+   &
     &            f_cfcst*fmass(iv2))/fmass(iv1)

              IF (fmass(iv2)-f_cfcst*fmass(iv2).gt.fmass(iv1-1)         &
     &            .and. fmass(iv2)-f_cfcst*fmass(iv2).le.               &
     &            fmass(iv1)) THEN

                 f_weight=((fmass(iv2)-f_cfcst*fmass(iv2)-              &
     &                     fmass(iv1-1))/(fmass(iv1)-fmass(iv1-1)))

              ELSEIF (fmass(iv2)-f_cfcst*fmass(iv2).gt.fmass(iv1)       &
     &            .and. fmass(iv2)-f_cfcst*fmass(iv2).lt.               &
     &            fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_sp
                 ELSE

                    f_weight=1.0_sp-((fmass(iv2)-f_cfcst*fmass(iv2)-    &
     &                       fmass(iv1))/(fmass(iv1+1)-fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_sp
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &            (SEDFLOCS%f_coll_prob_sh(iv2,iv3))*(fmass(iv2)-   &
     &            f_cfcst*fmass(iv2))/fmass(iv1)


              ! second case : iv3 eroded and iv2 eroded forming 3 particles : 
	      !iv3-f_cfcst*iv3 / iv2-f_cfcst*iv2 / f_cfcst*iv3+f_cfcst*iv2
           ELSEIF (gcolfragmin.ge.1.0_sp .and. gcolfragmax.ge.          &
     &            1.0_sp) THEN  

              IF (f_cfcst*fmass(iv2)+f_cfcst*fmass(iv3).gt.             &
     &            fmass(iv1-1) .and. f_cfcst*fmass(iv2)+f_cfcst*        &
     &            fmass(iv3).le.fmass(iv1)) THEN

                 f_weight=((f_cfcst*fmass(iv2)+f_cfcst*fmass(iv3)-      &
     &                   fmass(iv1-1))/(fmass(iv1)-fmass(iv1-1)))

              ELSEIF (f_cfcst*fmass(iv2)+f_cfcst*fmass(iv3).gt.         &
     &            fmass(iv1) .and. f_cfcst*fmass(iv2)+f_cfcst*          & 
     &            fmass(iv3).lt.fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_sp
                 ELSE

                    f_weight=1.0_sp-((f_cfcst*fmass(iv2)+f_cfcst*       &
     &                 fmass(iv3)-fmass(iv1))/(fmass(iv1+1)-            &
     &                 fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_sp
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &            (SEDFLOCS%f_coll_prob_sh(iv2,iv3))*(f_cfcst*      &
     &            fmass(iv2)+f_cfcst*fmass(iv3))/fmass(iv1)

              IF ((1.0_sp-f_cfcst)*fmass(iv2).gt.fmass(iv1-1)           &
     &            .and. (1.0_sp-f_cfcst)*fmass(iv2).le.                 &
     &            fmass(iv1)) THEN

                 f_weight=((1.0_sp-f_cfcst)*fmass(iv2)-                 &
     &                    fmass(iv1-1))/(fmass(iv1)-fmass(iv1-1))

              ELSEIF ((1.0_sp-f_cfcst)*fmass(iv2).gt.fmass(iv1)         &
     &            .and. (1.0_sp-f_cfcst)*fmass(iv2).lt.                 &
     &            fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_sp
                 ELSE

                    f_weight=1.0_sp-(((1.0_sp-f_cfcst)*fmass(iv2)-      &
     &                      fmass(iv1))/(fmass(iv1+1)-fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_sp
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &                  (SEDFLOCS%f_coll_prob_sh(iv2,iv3))*         &
     &                  ((1.0_sp-f_cfcst)*fmass(iv2))/fmass(iv1) 


              IF ((1.0_sp-f_cfcst)*fmass(iv3).gt.fmass(iv1-1) .and.     &
     &            (1.0_sp-f_cfcst)*fmass(iv3).le.fmass(iv1)) THEN

                 f_weight=((1.0_sp-f_cfcst)*fmass(iv3)-fmass(iv1-1))    &
     &                   /(fmass(iv1)-fmass(iv1-1))

              ELSEIF ((1.0_sp-f_cfcst)*fmass(iv3).gt.fmass(iv1)         &
     &            .and. (1.0_sp-f_cfcst)*fmass(iv3).lt.                 &
     &            fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_sp
                 ELSE

                    f_weight=1.0_sp-(((1.0_sp-f_cfcst)*fmass(iv3)-      &
     &                       fmass(iv1))/(fmass(iv1+1)-fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_sp
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &                  (SEDFLOCS%f_coll_prob_sh(iv2,iv3))*         &
     &                  ((1.0_sp-f_cfcst)*fmass(iv3))/fmass(iv1)


           ENDIF ! end collision test case
        ENDDO
       ENDDO
      ENDDO

      DO iv1=1,NCS
       DO iv2=1,NCS

        gcolfragiv1=2.0_sp*(Gval*(fdiam(iv1)+fdiam(iv2)))**2.0_sp*      &
     &            fmass(iv1)*fmass(iv2)/(pi*f_fy*f_fp*fdiam(iv1)        &
     &            **2.0_sp*(fmass(iv1)+fmass(iv2))*                     &
     &            ((SEDFLOCS%f_rho(iv1)-                            &
     &            rhoref)*cff2)**cff1)

        gcolfragiv2=2.0_sp*(Gval*(fdiam(iv1)+fdiam(iv2)))**2.0_sp*      &
     &            fmass(iv1)*fmass(iv2)/(pi*f_fy*f_fp*fdiam(iv2)        &
     &            **2.0_sp*(fmass(iv1)+fmass(iv2))*                     &
     &            ((SEDFLOCS%f_rho(iv2)-                            &
     &            rhoref)*cff2)**cff1)    

        mult=1.0_sp
        IF (iv1.eq.iv2) mult=2.0_sp

        IF (iv1.eq.MAX(iv1,iv2) .and. gcolfragiv1.ge.1.0_sp) THEN
           f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*                            &
     &            (SEDFLOCS%f_coll_prob_sh(iv1,iv2))
        ELSEIF (iv1.eq.MIN(iv1,iv2) .and. gcolfragiv2.ge.1.0_sp) THEN
           f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*                            &
     &            (SEDFLOCS%f_coll_prob_sh(iv1,iv2))
        ENDIF

       ENDDO
      ENDDO

      f_g4(1:NCS,1:NCS,1:NCS)=                                          &
     &            f_g4(1:NCS,1:NCS,1:NCS)*f_collfragparam
      f_l4(1:NCS,1:NCS)=f_l4(1:NCS,1:NCS)*f_collfragparam


      if(dbg_set(dbg_sbr)) write(ipt,*) "End:  flocmod_collfrag"

      RETURN
      END SUBROUTINE flocmod_collfrag

!!===========================================================================
      SUBROUTINE flocmod_comp_fsd(NNin,NNout,Gval,f_dt)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_comp_fsd  ***
  !&E
  !&E ** Purpose : computation of floc size distribution
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !
      implicit none 

  !! * Arguments
      real(sp),intent(in) :: Gval
      real(sp),dimension(1:NCS),intent(in)  :: NNin
      real(sp),dimension(1:NCS),intent(out) :: NNout
      real(sp),intent(in) :: f_dt

  !! * Local declarations
      integer      :: iv1,iv2,iv3
      real(sp) :: tmp_g1,tmp_g3,tmp_l1,tmp_l3,tmp_l4,tmp_g4
      real(sp),dimension(1:NCS,1:NCS,1:NCS)     :: f_g1_tmp,f_g4
      real(sp),dimension(1:NCS,1:NCS)           :: f_l1_tmp,f_l4
      
  !!--------------------------------------------------------------------------
  !! * Executable part


      if(dbg_set(dbg_sbr)) write(ipt,*) "Start: flocmod_comp_fsd "

      tmp_g1=0.0_sp
      tmp_g3=0.0_sp
      tmp_g4=0.0_sp
      tmp_l1=0.0_sp
      tmp_l3=0.0_sp
      tmp_l4=0.0_sp    
      f_g1_tmp(1:NCS,1:NCS,1:NCS)=0.0_sp
      f_l1_tmp(1:NCS,1:NCS)=0.0_sp

      IF (l_COLLFRAG) CALL flocmod_collfrag(Gval,f_g4,f_l4)

      DO iv1=1,NCS
       DO iv2=1,NCS
        DO iv3=1,NCS
           IF (l_ASH) THEN
              f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+              &
     &            SEDFLOCS%f_g1_sh(iv2,iv3,iv1)*Gval   
           ENDIF
           IF (l_ADS) THEN
              f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+              &
     &            SEDFLOCS%f_g1_ds(iv2,iv3,iv1)   
           ENDIF

           tmp_g1=tmp_g1+(NNin(iv3)*(f_g1_tmp(iv2,iv3,iv1))*NNin(iv2))

           IF (l_COLLFRAG) THEN
              tmp_g4=tmp_g4+(NNin(iv3)*                                 &
     &            (SEDFLOCS%f_g4(iv2,iv3,iv1)*Gval)*NNin(iv2))
           ENDIF
        ENDDO

        tmp_g3=tmp_g3+SEDFLOCS%f_g3(iv2,iv1)*NNin(iv2)*             &
     &            Gval**1.5_sp

        IF (l_ASH) THEN
           f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+                         &
     &            SEDFLOCS%f_l1_sh(iv2,iv1)*Gval
        ENDIF
        IF (l_ADS) THEN
           f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+                         &
     &            SEDFLOCS%f_l1_ds(iv2,iv1)*Gval
        ENDIF

        tmp_l1=tmp_l1+(f_l1_tmp(iv2,iv1))*NNin(iv2)

        IF (l_COLLFRAG) THEN
           tmp_l4=tmp_l4+(SEDFLOCS%f_l4(iv2,iv1)*Gval)*NNin(iv2)
        ENDIF
       ENDDO

       tmp_l1=tmp_l1*NNin(iv1)
       tmp_l4=tmp_l4*NNin(iv1)
       tmp_l3=SEDFLOCS%f_l3(iv1)*Gval**1.5_sp*NNin(iv1)

       NNout(iv1)=NNin(iv1)+ f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+        &
     &            tmp_l3+tmp_l4))

       tmp_g1=0.0_sp
       tmp_g3=0.0_sp
       tmp_g4=0.0_sp
       tmp_l1=0.0_sp
       tmp_l3=0.0_sp
       tmp_l4=0.0_sp    
      ENDDO

      if(dbg_set(dbg_sbr)) write(ipt,*) "End: flocmod_comp_fsd "
      RETURN
      END SUBROUTINE flocmod_comp_fsd



!!===========================================================================
      SUBROUTINE flocmod_mass_control(NN,mneg)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_mass_control  ***
  !&E
  !&E ** Purpose : Compute mass in every class after flocculation and 
  !&E              returns negative mass if any
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------

!
      implicit none 
 
  !! * Local declarations
      integer      :: iv1
      real(sp),intent(out)     :: mneg
      real(sp),dimension(1:NCS),intent(in)     :: NN
      !real(sp),DIMENSION(0:NCS+1)      :: f_mass

  !!--------------------------------------------------------------------------
  !! * Executable part
      if(dbg_set(dbg_sbr)) write(ipt,*) "Start: flocmod_mass_control "
      mneg=0.0_sp

      DO iv1=1,NCS
       IF (NN(iv1).lt.0.0_sp) THEN
         mneg=mneg-NN(iv1)*SEDFLOCS%f_mass(iv1)
       ENDIF
      ENDDO

      if(dbg_set(dbg_sbr)) write(ipt,*) "End: flocmod_mass_control "
      RETURN
      END SUBROUTINE flocmod_mass_control

!!===========================================================================
      SUBROUTINE flocmod_mass_redistribute(NN)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_mass_redistribute  ***
  !&E
  !&E ** Purpose : based on a tolerated negative mass parameter, negative masses  
  !&E              are redistributed linearly towards remaining postive masses 
  !&E              and negative masses are set to 0
  !&E                   
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !
      implicit none 


  !! * Local declarations
      integer      :: iv
      real(sp)     :: npos
      real(sp)     :: mneg
      real(sp),dimension(1:NCS),intent(inout)     :: NN
      real(sp),dimension(1:NCS)                   :: NNtmp
      !real(sp),DIMENSION(0:NCS+1)      :: f_mass
  !!--------------------------------------------------------------------------
  !! * Executable part

      if(dbg_set(dbg_sbr)) write(ipt,*) "Start:  flocmod_mass_redistribute"

      mneg=0.0_sp
      npos=0.0_sp
      NNtmp(:)=NN(:)

      DO iv=1,NCS
       IF (NN(iv).lt.0.0_sp) THEN
         mneg=mneg-NN(iv)*SEDFLOCS%f_mass(iv)
         NNtmp(iv)=0.0_sp
       ELSE
         npos=npos+1.0_sp
       ENDIF
      ENDDO

      IF (mneg.gt.0.0_sp) THEN 
       IF (npos.eq.0.0_sp) THEN
         WRITE(*,*) 'CAUTION : all floc sizes have negative mass!'
         WRITE(*,*) 'SIMULATION STOPPED'
         STOP    
       ELSE
         DO iv=1,NCS
           IF (NN(iv).gt.0.0_sp) THEN
              NN(iv)=NN(iv)-mneg/sum(NNtmp)*NN(iv)/                     &
     &            SEDFLOCS%f_mass(iv)
           ELSE
              NN(iv)=0.0_sp
           ENDIF

         ENDDO

       ENDIF
      ENDIF

      if(dbg_set(dbg_sbr)) write(ipt,*) "End:  flocmod_mass_redistribute"
      RETURN
      END SUBROUTINE flocmod_mass_redistribute

!!===========================================================================    
  SUBROUTINE flocmod_comp_g(k,i,Gval,diss)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_comp_g  ***
  !&E
  !&E ** Purpose : compute shear rate to estimate shear aggregation and erosion  
  !&E 
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !
  implicit none 

  !! * Local declarations
  integer,  intent(in)      :: k,i
  real(sp),intent(out)     :: Gval
  real(sp)     :: htn,ustar,z,diss,nueau
! l_testcase - if .TRUE. sets G(t) to values from lab experiment
!      logical, parameter :: l_testcase = .TRUE.
  !!--------------------------------------------------------------------------
  !! * Executable part

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: flocmod_comp_g "

  !    nueau=1.06e-6_sp
  !    ustar=sqrt(tenfon(i)/rhoref)
  !    htn=h0(i)+ssh(i)
  !    z=(1.0_sp+sig(k))*htn

!
! ALA from CRS
!
  nueau = 1.5E-6_sp
  IF (l_testcase) THEN
    ! reproducing flocculation experiment Verney et al., 2011
!    Gval=0.0_sp
!    IF (time .lt. 7201.0_sp) THEN
!      Gval=1.0_sp
!    ELSEIF (time .lt. 8401.0_sp) THEN
!      Gval=2.0_sp
!    ELSEIF (time .lt. 9601.0_sp) THEN
!      Gval=3.0_sp  
!    ELSEIF (time .lt. 10801.0_sp) THEN
!      Gval=4.0_sp
!    ELSEIF (time .lt. 12601.0_sp) THEN
!      Gval=12.0_sp
!    ELSEIF (time .lt. 13801.0_sp) THEN
!      Gval=4.0_sp
!    ELSEIF (time .lt. 15001.0_sp) THEN
!      Gval=3.0_sp
!    ELSEIF (time .lt. 16201.0_sp) THEN
!      Gval=2.0_sp
!    ELSEIF (time .lt. 21601.0_sp) THEN
!      Gval=1.0_sp
!    ELSEIF (time .lt. 25201.0_sp) THEN
!      Gval=0.0_sp
!    ELSEIF (time .lt. 30601.0_sp) THEN
!      Gval=1.0_sp
!    ELSEIF (time .lt. 31801.0_sp) THEN
!      Gval=2.0_sp                     
!    ELSEIF (time .lt. 33001.0_sp) THEN
!      Gval=3.0_sp       
!    ELSEIF (time .lt. 34201.0_sp) THEN
!      Gval=4.0_sp
!    ELSEIF (time .lt. 36001.0_sp) THEN
!      Gval=12.0_sp
!    ELSEIF (time .lt. 37201.0_sp) THEN
!      Gval=4.0_sp
!    ELSEIF (time .lt. 38401.0_sp) THEN
!      Gval=3.0_sp
!    ELSEIF (time .lt. 39601.0_sp) THEN
!      Gval=2.0_sp
!    ELSEIF (time .lt. 45001.0_sp) THEN
!      Gval=1.0_sp
!    ELSEIF (time .lt. 48601.0_sp) THEN
!      Gval=0.0_sp
!    ELSEIF (time .lt. 54001.0_sp) THEN
!      Gval=1.0_sp
!    ELSEIF (time .lt. 55201.0_sp) THEN
!      Gval=2.0_sp                     
!    ELSEIF (time .lt. 56401.0_sp) THEN
!      Gval=3.0_sp       
!    ELSEIF (time .lt. 57601.0_sp) THEN
!      Gval=4.0_sp
!    ELSE 
!      Gval=12.0_sp
!    ENDIF
  ELSE
    Gval=sqrt(diss/nueau)
  ENDIF
! NO KLUDGE
!      Gval = 12.0_sp

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: flocmod_comp_g "
  RETURN

  END SUBROUTINE flocmod_comp_g
      
      
# endif

End Module Mod_FlocMod

