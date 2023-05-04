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
!
! Note:  You must link with gotm 4.x libraries/modules  
!
!/===========================================================================/

MODULE MOD_GOTM 
#  if defined (GOTM)
   USE MOD_TYPES

# if defined (VEGETATION)
   USE MOD_VEGETATION
# endif

   IMPLICIT NONE
   SAVE

   CONTAINS !------------------------------------------------------------------!
            ! INIT_GOTM           :   INITIALIZE TMODEL USING GOTM LIBRARIES   !
            ! ADVANCE_GOTM        :   ADVANCE TMODEL USING GOTM LIBRARIES      !
            !------------------------------------------------------------------!

!==============================================================================|
   SUBROUTINE INIT_GOTM
!==============================================================================|
!  A WRAPPER ROUTINE TO INITIALIZE TMODEL USING GOTM LIBRARIES                 |
!==============================================================================|
   use lims,         only: kbm1
   use control,      only: input_dir,casename,msr,ipt
   use turbulence,   only: init_turbulence
   use mtridiagonal, only: init_tridiagonal
   USE MOD_UTILS
   IMPLICIT NONE
   CHARACTER(LEN=80) :: FNAME
   LOGICAL           :: FEXIST
   INTEGER, PARAMETER :: igotm = 59
!------------------------------------------------------------------------------!

   FNAME = TRIM(INPUT_DIR)//"/"//trim(casename)//'_gotmturb.inp'
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'GOTM PARAMETER FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF
   IF(MSR)THEN
     WRITE(IPT,*)'============== INITIATING GOTM ==============================='
     WRITE(IPT,*)
   END IF
   CALL init_turbulence(igotm,trim(fname),kbm1)
   CALL init_tridiagonal(kbm1)
   IF(MSR)THEN
     WRITE(IPT,*)'============== GOTM INITIATION COMPLETE======================='
     WRITE(IPT,*)
   END IF

   RETURN
   END SUBROUTINE INIT_GOTM
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!==============================================================================|
   SUBROUTINE ADVANCE_GOTM
!==============================================================================|
!  A WRAPPER ROUTINE TO CALL THE GOTM LIBRARIES                                |
!==============================================================================|

!--Variables from FVCOM Modules------------------------------
   use lims, only : m,kbm1,kb
   use all_vars, only : rho1,dz,dzz,d,t1,u,v,km,kh,tke,teps,wusurf,wvsurf,l,  &
                        z,zz,nbve,ntve,grav_n,cc_z0b
   use control, only: dti,umol,iint,iend

!--Variables from GOTM Modules-------------------------------
   use turbulence, only: do_turbulence,cde
   use turbulence, only: tke1d => tke, eps1d => eps, L1d => L
   use turbulence, only: num,nuh 

!--Local Variables (Temporary)-------------------------------
   IMPLICIT NONE
   INTEGER                   :: I,J,K
   REAL(SP)                  :: rr,ztemp

!--Variables for Interfacing with GOTM-----------------------
   DOUBLE PRECISION          :: depth,time_step
   DOUBLE PRECISION          :: u_taus,u_taub,z0s_gotm,z0b_gotm,umold
   DOUBLE PRECISION          :: h(0:KBM1)
   REAL(SP)                  :: NN(M,1:KB),SS(M,1:KB)
   DOUBLE PRECISION          :: NN1d(0:KBM1),SS1d(0:KBM1)
!--Parameters------------------------------------------------
   REAL(SP), PARAMETER       :: KAPPA  = .4          !!Von Karman's Constant
   REAL(SP), PARAMETER       :: CHARNOK_VAL = 1400.  !!Charnok Constant 
   REAL(SP), PARAMETER       :: z0s_min = .00        !!Minimum Surface Roughness 
   INTEGER , PARAMETER       :: MaxItz0b = 10

!------------------------------------------------------------
   REAL(SP) :: UU_NODE_K,VV_NODE_K,UU_NODE_KM1,VV_NODE_KM1
   REAL(SP) :: UU_NODE_KBM1,VV_NODE_KBM1
   REAL(SP) :: WUSURF_NODE,WVSURF_NODE

!==============================================================================|


!------------------Calculate Buoyancy Frequency Squared (NN)-------------------!

   DO K=2,KBM1
     DO I=1,M
       NN(I,K) = -GRAV_N(I) * (RHO1(I,K-1)-RHO1(I,K))/(DZZ(I,K-1)*D(I))
     END DO
   END DO
   !!Set BC's as GOTM Does
   NN(:,1)  = NN(:,2)   
   NN(:,KB) = NN(:,KBM1)

!------------------Calculate Shear Frequency Squared (SS)----------------------!

   DO K=2,KBM1
     DO I=1,M
       UU_NODE_K   = 0.0_SP
       VV_NODE_K   = 0.0_SP
       UU_NODE_KM1 = 0.0_SP
       VV_NODE_KM1 = 0.0_SP
       DO J=1,NTVE(I)
         UU_NODE_K   = UU_NODE_K + U(NBVE(I,J),K)
         VV_NODE_K   = VV_NODE_K + V(NBVE(I,J),K)
         UU_NODE_KM1 = UU_NODE_KM1 + U(NBVE(I,J),K-1)
         VV_NODE_KM1 = VV_NODE_KM1 + V(NBVE(I,J),K-1)
       END DO
       UU_NODE_K   = UU_NODE_K/FLOAT(NTVE(I))
       VV_NODE_K   = VV_NODE_K/FLOAT(NTVE(I))
       UU_NODE_KM1 = UU_NODE_KM1/FLOAT(NTVE(I))
       VV_NODE_KM1 = VV_NODE_KM1/FLOAT(NTVE(I))

!       SS(I,K) = ((U(I,K)-U(I,K-1))**2+(V(I,K)-V(I,K-1))**2)/(DZZ(I,K-1)*D(I))**2 
       SS(I,K) = ((UU_NODE_K-UU_NODE_KM1)**2+(VV_NODE_K-VV_NODE_KM1)**2)/(DZZ(I,K-1)*D(I))**2 
     END DO
   END DO

   !!Set BC's as GOTM Does
   SS(:,1)  = SS(:,2)
   SS(:,KB) = SS(:,KBM1)

!------------------Main Loop Over Elements-------------------------------------!

   DO I=1,M
      !Surface Friction Velocity [m/s]
     WUSURF_NODE = 0.0_SP
     WVSURF_NODE = 0.0_SP
     DO J=1,NTVE(I)
       WUSURF_NODE = WUSURF_NODE + WUSURF(NBVE(I,J))
       WVSURF_NODE = WVSURF_NODE + WVSURF(NBVE(I,J))
     END DO
     WUSURF_NODE = WUSURF_NODE/FLOAT(NTVE(I))
     WVSURF_NODE = WVSURF_NODE/FLOAT(NTVE(I))

!      u_taus=(wusurf(i)**2+wvsurf(i)**2)**(1./4.)
      u_taus=(wusurf_node**2+wvsurf_node**2)**(1./4.)
  
      !Set Surface Roughness Height Using Charnok Formula[m]
!      z0s_gotm = charnok_val*u_taus**2/grav
!      if(z0s_gotm < z0s_min) z0s_gotm = z0s_min
      z0s_gotm = z0s_min
      
      !Set Bottom Roughness Height [m] and Friction Velocity [m/s] using GOTM Method 
      ztemp=(zz(I,kbm1)-z(I,kb))*d(i)
      umold = umol

      UU_NODE_KBM1 = 0.0_SP
      VV_NODE_KBM1 = 0.0_SP
      DO J=1,NTVE(I)
        UU_NODE_KBM1 = UU_NODE_KBM1 + U(NBVE(I,J),KBM1)
        VV_NODE_KBM1 = VV_NODE_KBM1 + V(NBVE(I,J),KBM1)
      END DO
      UU_NODE_KBM1 = UU_NODE_KBM1/FLOAT(NTVE(I))
      VV_NODE_KBM1 = VV_NODE_KBM1/FLOAT(NTVE(I))

      u_taub = 0.0_SP
      DO J=1,MaxItz0b
        z0b_gotm = 0.1*umol/max(umold,u_taub)+0.03*(cc_z0b(i)/.03)
        rr = kappa/(log((z0b_gotm+ztemp)/z0b_gotm)) 
!        u_taub = rr*sqrt( u(i,kbm1)*u(i,kbm1) + v(i,kbm1)*v(i,kbm1) )
        u_taub = rr*sqrt(uu_node_kbm1*uu_node_kbm1+vv_node_kbm1*vv_node_kbm1)
      END DO
      
      !Set up Depth [m] 
      depth = d(i)

      !Set up Time Step [s] 
      time_step = dti

      !Set up Layer Thicknesses [m]
      h(1:kbm1) = dz(I,kbm1:1:-1)*depth 
      
      !Set Up 1-D Arrays                                      
      NUM(0:KBM1)   = KM(I,KB:1:-1)    !Vertical Kinematic Viscosity  [m^2/s]
      NUH(0:KBM1)   = KH(I,KB:1:-1)    !Vertical Kinematic Viscosity  [m^2/s]
      SS1d(0:KBM1)  = SS(I,KB:1:-1)    !Shear Frequency Squared       [1/s^2]        
      NN1d(0:KBM1)  = NN(I,KB:1:-1)    !Buoyancy Frequency Squared    [1/s^2]
      tke1d(0:KBM1) = tke(I,KB:1:-1)   !Turbulent Kinetic Energy      [m^2/s^2]
      eps1d(0:KBM1) = teps(I,KB:1:-1)  !Turbulence Dissipation Rate   [m^2/s^3] 
      l1d(0:kbm1)   = l(i,kb:1:-1)

# if defined (VEGETATION)
      !
      !  Add the effect of vegetation on tke and gls 
      ! 
      if(VEGETATION_MODEL.and.VEG_TURB)then
        tke1d(0:KBM1) = tke1d(0:KBM1) +  VEG % tke_veg(i,KB:1:-1)
        eps1d(0:KBM1) = eps1d(0:KBM1) +  VEG % gls_veg(i,KB:1:-1)
      end if
      !
      !-----------------------------------------------------------------------
      !  Add the effect of vegetation on horizontal viscosity.
      !-----------------------------------------------------------------------
      !
      if(VEGETATION_MODEL.and.VEG_HMIXING)then
          NUM(0:KBM1) = NUM(0:KBM1) + VEG % visc3d_r_veg(i,KB:1:-1)
      end if

# endif

      !Update Turbulence Model
      call do_turbulence(kbm1,time_step,depth,u_taus,u_taub,z0s_gotm,z0b_gotm,h,nn1D,ss1D) 

      !Update 3D Fields of TKE,EPS,KM,KH
       tke(I,1:KB) = tke1d(KBM1:0:-1)
      teps(I,1:KB) = eps1d(KBM1:0:-1)
      teps(I,1)    = 0.0 !required, as gotm can produce "Inf" which will crash NetCDF on output 
        km(I,1:KB) = num(KBM1:0:-1)
        kh(I,1:KB) = nuh(KBM1:0:-1)
         l(i,1:kb) = l1d(kbm1:0:-1)
   END DO


   RETURN
   END SUBROUTINE ADVANCE_GOTM
!==============================================================================|
#  endif
END MODULE MOD_GOTM
