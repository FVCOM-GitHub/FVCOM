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

!=======================================================================
! FVCOM Scalar Module  
!
!    contains methods:
!        Adv_Scal            => Advect a Scalar Quantity 
!        Vdif_Scal           => Vertical Diffusion of Scalar Quantity
!        Bcond_Scal_OBC      => Open Boundary Condition for Scalar
!        Bcond_Scal_PTsource => Point Sources of Scalar
!=======================================================================
Module Scalar

  logical, parameter :: debug = .true. 

  contains
!==============================================================================|
! Calculate Horizontal Advection and Diffusion For Scalar (f)                  |
!==============================================================================|
!  Subroutine Adv_Scal(f,fn,d_fdis,fdis,d_fflux,fflux_obc,deltat,source)
!J. Ge for tracer advection
  Subroutine Adv_Scal(f,f0,f2,fn,d_fdis,fdis,d_fflux,fflux_obc,deltat,source)
!J. Ge for tracer advection

!------------------------------------------------------------------------------|

  use all_vars
  use lims, only: m,mt,n,nt,kbm1,kb
  use bcs
  use mod_obcs
# if defined (MULTIPROCESSOR)
  use mod_par
# endif
# if defined (WET_DRY)
  use mod_wd
# endif

# if defined (THIN_DAM)
  use mod_dam, only : kdam,N_DAM_MATCH,IS_DAM
# endif

  implicit none
  real(sp), intent(in ), dimension(0:mt,kb)      :: f
!J. Ge for tracer advection
  real(sp), intent(in ), dimension(0:mt,kb)      :: f0,f2
!J. Ge for tracer advection 
  real(sp), intent(out), dimension(0:mt,kb)      :: fn
  integer , intent(in )                          :: d_fdis
  real(sp), intent(in ), dimension(d_fdis)       :: fdis
  integer , intent(in )                          :: d_fflux
  real(sp), intent(out), dimension(d_fflux,kbm1) :: fflux_obc 
  real(sp), intent(in )                          :: deltat
  logical , intent(in )                          :: source

  !----------------local--------------------------------------
  real(sp), dimension(0:mt,kb)   :: xflux,xflux_adv
  real(sp), dimension(m)         :: pupx,pupy,pvpx,pvpy  
  real(sp), dimension(m)         :: pfpx,pfpy,pfpxd,pfpyd,viscoff
  real(sp), dimension(3*nt,kb)      :: dtij 
  real(sp), dimension(3*nt,kbm1) :: uvn
  real(sp), dimension(kb)        :: vflux
  real(sp) :: utmp,vtmp,sitai,ffd,ff1,x11,y11,x22,y22,x33,y33
  real(sp) :: tmp1,tmp2,xi,yi
  real(sp) :: dxa,dya,dxb,dyb,fij1,fij2,un
  real(sp) :: txx,tyy,fxx,fyy,viscof,exflux,temp,fpoint
  real(sp) :: fact,fm1,fmean
  integer  :: i,i1,i2,ia,ib,j,j1,j2,k,jtmp,jj
!J. Ge for tracer advectio
  real(sp),dimension(kb) :: ftmp
!J. Ge for tracer advection

# if defined (SPHERICAL)
  real(sp) :: ty,txpi,typi
# else
  real(sp) :: da,db,ds
# endif

# if defined (THIN_DAM)
  INTEGER  :: NX
  real(sp) :: tmpflx
  real(sp),dimension(kb) :: wvel
# endif


!------------------------------------------------------------------------------!

!-------------------------------------------------------
!Calculate Mean Values
!-------------------------------------------------------

  fmean = sum(f(1:m,1:kbm1))/float(m*kbm1)

!-------------------------------------------------------
!Initialize Multipliers to Control Horizontal Diff
!-------------------------------------------------------

  fact = 0.0_sp
  fm1  = 1.0_sp
  if(HORIZONTAL_MIXING_TYPE == 'closure') then
    fact = 1.0_sp
    fm1  = 0.0_sp
  end if
     
!-------------------------------------------------------
!Initialize Fluxes
!-------------------------------------------------------
  xflux     = 0.0_sp
  xflux_adv = 0.0_sp

!-------------------------------------------------------
!Calculate Normal Velocity on Control Volume Edges
!-------------------------------------------------------
!!# if !defined (WET_DRY)
  do i=1,ncv
    i1=ntrg(i)
    !dtij(i)=dt1(i1)
    do k=1,kbm1
      dtij(i,k)=dt1(i1)*dz1(i1,k) 
      uvn(i,k) = v(i1,k)*dltxe(i) - u(i1,k)*dltye(i)

#  if defined(PLBC)
      uvn(i,k) =  - u(i1,k)*dltye(i)
#  endif

    end do
  end do
!!# else
!!  do i=1,ncv
!!    i1=ntrg(i)
!!    dtij(i)=dt1(i1)
!!    do k=1,kbm1
!!      uvn(i,k) = vs(i1,k)*dltxe(i) - us(i1,k)*dltye(i)
!!    end do
!!  end do
!!# endif

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   do k=1,kbm1
      pfpx  = 0.0_sp 
      pfpy  = 0.0_sp 
      pfpxd = 0.0_sp 
      pfpyd = 0.0_sp
     do i=1,m
       do j=1,ntsn(i)-1
         i1=nbsn(i,j)
         i2=nbsn(i,j+1)

#    if defined (WET_DRY)
!         IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
!          FFD=0.5_SP*(f(I,K)+f(I2,K))
!          FF1=0.5_SP*(f(I,K)+f(I2,K))
!	 ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
!          FFD=0.5_SP*(f(I1,K)+f(I,K))
!          FF1=0.5_SP*(f(I1,K)+f(I,K))
!	 ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
!          FFD=0.5_SP*(f(I,K)+f(I,K))
!          FF1=0.5_SP*(f(I,K)+f(I,K))
!	 ELSE
!          FFD=0.5_SP*(f(I1,K)+f(I2,K))
!          FF1=0.5_SP*(f(I1,K)+f(I2,K))
!	 END IF
 
!J. Ge for tracer advection
         IF(BACKWARD_ADVECTION .NEQV. .TRUE.)THEN
           IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
            FFD=0.5_SP*(f(I,K)+f(I2,K))
            FF1=0.5_SP*(f(I,K)+f(I2,K))
	   ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
            FFD=0.5_SP*(f(I1,K)+f(I,K))
            FF1=0.5_SP*(f(I1,K)+f(I,K))
	   ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
            FFD=0.5_SP*(f(I,K)+f(I,K))
            FF1=0.5_SP*(f(I,K)+f(I,K))
	   ELSE
            FFD=0.5_SP*(f(I1,K)+f(I2,K))
            FF1=0.5_SP*(f(I1,K)+f(I2,K))
	   END IF
         ELSE
           IF(BACKWARD_STEP==1)THEN
             IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
               FFD=0.5_SP*((f0(I,K)+f(I,K))*0.5+(f0(I2,K)+f(I2,K))*0.5)
               FF1=0.5_SP*((f0(I,K)+f(I,K))*0.5+(f0(I2,K)+f(I2,K))*0.5)
             ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((f0(I1,K)+f(I1,K))*0.5+(f0(I,K)+f(I,K))*0.5)
               FF1=0.5_SP*((f0(I1,K)+f(I1,K))*0.5+(f0(I,K)+f(I,K))*0.5)
	     ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((f0(I,K)+f(I,K))*0.5+(f0(I,K)+f(I,K))*0.5)
               FF1=0.5_SP*((f0(I,K)+f(I,K))*0.5+(f0(I,K)+f(I,K))*0.5)
	     ELSE
               FFD=0.5_SP*((f0(I1,K)+f(I1,K))*0.5+(f0(I2,K)+f(I2,K))*0.5)
               FF1=0.5_SP*((f0(I1,K)+f(I1,K))*0.5+(f0(I2,K)+f(I2,K))*0.5)
             END IF
           ELSEIF(BACKWARD_STEP==2)THEN
             IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
               FFD=0.5_SP*((f2(I,K)+f0(I,K)+f(I,K))/3.0_SP+(f2(I2,K)+f0(I2,K)+f(I2,K))/3.0_SP)
               FF1=0.5_SP*((f2(I,K)+f0(I,K)+f(I,K))/3.0_SP+(f2(I2,K)+f0(I2,K)+f(I2,K))/3.0_SP)
             ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((f2(I1,K)+f0(I1,K)+f(I1,K))/3.0_SP+(f2(I,K)+f0(I,K)+f(I,K))/3.0_SP)
               FF1=0.5_SP*((f2(I1,K)+f0(I1,K)+f(I1,K))/3.0_SP+(f2(I,K)+f0(I,K)+f(I,K))/3.0_SP)
	     ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((f2(I,K)+f0(I,K)+f(I,K))/3.0_SP+(f2(I,K)+f0(I,K)+f(I,K))/3.0_SP)
               FF1=0.5_SP*((f2(I,K)+f0(I,K)+f(I,K))/3.0_SP+(f2(I,K)+f0(I,K)+f(I,K))/3.0_SP)
	     ELSE
               FFD=0.5_SP*((f2(I1,K)+f0(I1,K)+f(I1,K))/3.0_SP+(f2(I2,K)+f0(I2,K)+f(I2,K))/3.0_SP)
               FF1=0.5_SP*((f2(I1,K)+f0(I1,K)+f(I1,K))/3.0_SP+(f2(I2,K)+f0(I2,K)+f(I2,K))/3.0_SP)
	     END IF  
           ENDIF

         END IF
!J. Ge for tracer advection
#    else	 
!         ffd=0.5_sp*(f(i1,k)+f(i2,k)) !-fmean1(i1,k)-fmean1(i2,k))
!         ff1=0.5_sp*(f(i1,k)+f(i2,k))
!J. Ge for tracer advection	 
         IF(BACKWARD_ADVECTION .NEQV. .TRUE.)THEN
           FFD=0.5_SP*(f(I1,K)+f(I2,K)) !-SMEAN1(I1,K)-SMEAN1(I2,K))
           FF1=0.5_SP*(f(I1,K)+f(I2,K))
         ELSE
           IF(BACKWARD_STEP==1)THEN
             FFD=0.5_SP*((f0(I1,K)+f(I1,K))*0.5+(f0(I2,K)+f(I2,K))*0.5)
             FF1=0.5_SP*((f0(I1,K)+f(I1,K))*0.5+(f0(I2,K)+f(I2,K))*0.5)
           ELSEIF(BACKWARD_STEP==2)THEN
             FFD=0.5_SP*((f2(I1,K)+f0(I1,K)+f(I1,K))*0.5+(f2(I2,K)+f0(I2,K)+f(I2,K))*0.5)
             FF1=0.5_SP*((f2(I1,K)+f0(I1,K)+f(I1,K))*0.5+(f2(I2,K)+f0(I2,K)+f(I2,K))*0.5)
           ENDIF
         ENDIF
!J. Ge for tracer advection
#    endif	 
	 
#        if defined (SPHERICAL)
         ty=0.5_sp*(vy(i1)+vy(i2))
         txpi=(vx(i2)-vx(i1))*tpi*cos(deg2rad*ty)
         typi=(vy(i1)-vy(i2))*tpi
         pfpx(i)=pfpx(i)+ff1*typi
         pfpy(i)=pfpy(i)+ff1*txpi
         pfpxd(i)=pfpxd(i)+ffd*typi
         pfpyd(i)=pfpyd(i)+ffd*txpi
#        else
         pfpx(i) = pfpx(i) +ff1*(vy(i1)-vy(i2))
         pfpy(i) = pfpy(i) +ff1*(vx(i2)-vx(i1))
         pfpxd(i)= pfpxd(i)+ffd*(vy(i1)-vy(i2))
         pfpyd(i)= pfpyd(i)+ffd*(vx(i2)-vx(i1))
#        endif
       end do

! gather all neighboring control volumes connecting at dam node 
# if defined (THIN_DAM)
       IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
         DO NX=1,N_DAM_MATCH(I,1)
           DO J=1,NTSN(N_DAM_MATCH(I,NX+1))-1
             I1=NBSN(N_DAM_MATCH(I,NX+1),J)
             I2=NBSN(N_DAM_MATCH(I,NX+1),J+1)
!             FFD=0.5_SP*(f(I1,K)+f(I2,K)) !-SMEAN1(I1,K)-SMEAN1(I2,K))
!             FF1=0.5_SP*(f(I1,K)+f(I2,K))
!J. Ge for tracer advection
             IF(BACKWARD_ADVECTION .NEQV. .TRUE.)THEN
               FFD=0.5_SP*(f(I1,K)+f(I2,K)) !-SMEAN1(I1,K)-SMEAN1(I2,K))
               FF1=0.5_SP*(f(I1,K)+f(I2,K))
             ELSE
               IF(BACKWARD_STEP==1)THEN
                FFD=0.5_SP*((f0(I1,K)+f(I1,K))*0.5+(f0(I2,K)+f(I2,K))*0.5)
                FF1=0.5_SP*((f0(I1,K)+f(I1,K))*0.5+(f0(I2,K)+f(I2,K))*0.5)
               ELSEIF(BACKWARD_STEP==2)THEN
                FFD=0.5_SP*((f2(I1,K)+f0(I1,K)+f(I1,K))/3.0_SP+(f2(I2,K)+f0(I2,K)+f(I2,K))/3.0_SP)
                FF1=0.5_SP*((f2(I1,K)+f0(I1,K)+f(I1,K))/3.0_SP+(f2(I2,K)+f0(I2,K)+f(I2,K))/3.0_SP)
               ENDIF
             ENDIF
!J. Ge for tracer advection
#        if defined (SPHERICAL)
             ty=0.5_sp*(vy(i1)+vy(i2))
             txpi=(vx(i2)-vx(i1))*tpi*cos(deg2rad*ty)
             typi=(vy(i1)-vy(i2))*tpi
             pfpx(i)=pfpx(i)+ff1*typi
             pfpy(i)=pfpy(i)+ff1*txpi
             pfpxd(i)=pfpxd(i)+ffd*typi
             pfpyd(i)=pfpyd(i)+ffd*txpi
#        else
             pfpx(i) = pfpx(i) +ff1*(vy(i1)-vy(i2))
             pfpy(i) = pfpy(i) +ff1*(vx(i2)-vx(i1))
             pfpxd(i)= pfpxd(i)+ffd*(vy(i1)-vy(i2))
             pfpyd(i)= pfpyd(i)+ffd*(vx(i2)-vx(i1))
#        endif
           END DO
         END DO
       END IF
# endif

# if !defined (THIN_DAM)
       pfpx(i)  =pfpx(i )/art2(i)
       pfpy(i)  =pfpy(i )/art2(i)
       pfpxd(i) =pfpxd(i)/art2(i)
       pfpyd(i) =pfpyd(i)/art2(i)
# else
       IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
         PFPX(I)=PFPX(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
         PFPY(I)=PFPY(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
         PFPXD(I)=PFPXD(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
         PFPYD(I)=PFPYD(I)/(ART2(I)+SUM(ART2(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))
       ELSE
         PFPX(I)=PFPX(I)/ART2(I)
         PFPY(I)=PFPY(I)/ART2(I)
         PFPXD(I)=PFPXD(I)/ART2(I)
         PFPYD(I)=PFPYD(I)/ART2(I)
       END IF
# endif

     end do
          
     if(k == kbm1)then
       do i=1,m
         pfpxb(i) = pfpx(i)
         pfpyb(i) = pfpy(i)
       end do
     end if

     do i=1,m
       pupx(i)=0.0_sp
       pupy(i)=0.0_sp
       pvpx(i)=0.0_sp
       pvpy(i)=0.0_sp
       j=1
       i1=nbve(i,j)
       jtmp=nbvt(i,j)
       j1=jtmp+1-(jtmp+1)/4*3
       j2=jtmp+2-(jtmp+2)/4*3
       x11=0.5_sp*(vx(i)+vx(nv(i1,j1)))
       y11=0.5_sp*(vy(i)+vy(nv(i1,j1)))
       x22=xc(i1)
       y22=yc(i1)
       x33=0.5_sp*(vx(i)+vx(nv(i1,j2)))
       y33=0.5_sp*(vy(i)+vy(nv(i1,j2)))

#      if defined (SPHERICAL)
       ty  =0.5_sp*(y11+y33)
       txpi=(x33-x11)*tpi*cos(deg2rad*ty)
       typi=(y11-y33)*tpi
       pupx(i)=pupx(i)+u(i1,k)*typi 
       pupy(i)=pupy(i)+u(i1,k)*txpi
       pvpx(i)=pvpx(i)+v(i1,k)*typi
       pvpy(i)=pvpy(i)+v(i1,k)*txpi
#      else
       pupx(i)=pupx(i)+u(i1,k)*(y11-y33)
       pupy(i)=pupy(i)+u(i1,k)*(x33-x11)
       pvpx(i)=pvpx(i)+v(i1,k)*(y11-y33)
       pvpy(i)=pvpy(i)+v(i1,k)*(x33-x11)
#      endif

       if(isonb(i) /= 0) then
#        if defined (SPHERICAL)
         ty=0.5_sp*(vy(i)+y11)
         txpi=(x11-vx(i))*tpi*cos(deg2rad*ty)
         typi=(vy(i)-y11)*tpi
         pupx(i)=pupx(i)+u(i1,k)*typi
         pupy(i)=pupy(i)+u(i1,k)*txpi
         pvpx(i)=pvpx(i)+v(i1,k)*typi
         pvpy(i)=pvpy(i)+v(i1,k)*txpi
#        else
         pupx(i)=pupx(i)+u(i1,k)*(vy(i)-y11)
         pupy(i)=pupy(i)+u(i1,k)*(x11-vx(i))
         pvpx(i)=pvpx(i)+v(i1,k)*(vy(i)-y11)
         pvpy(i)=pvpy(i)+v(i1,k)*(x11-vx(i))
#        endif
       end if

       do j=2,ntve(i)-1
         i1=nbve(i,j)
         jtmp=nbvt(i,j)
         j1=jtmp+1-(jtmp+1)/4*3
         j2=jtmp+2-(jtmp+2)/4*3
         x11=0.5_sp*(vx(i)+vx(nv(i1,j1)))
         y11=0.5_sp*(vy(i)+vy(nv(i1,j1)))
         x22=xc(i1)
         y22=yc(i1)
         x33=0.5_sp*(vx(i)+vx(nv(i1,j2)))
         y33=0.5_sp*(vy(i)+vy(nv(i1,j2)))

#        if defined (SPHERICAL)
         ty=0.5_sp*(y11+y33)
         txpi=(x33-x11)*tpi*COS(deg2rad*TY)
         typi=(y11-y33)*tpi
         pupx(i)=pupx(i)+u(i1,k)*typi
         pupy(i)=pupy(i)+u(i1,k)*txpi
         pvpx(i)=pvpx(i)+v(i1,k)*typi
         pvpy(i)=pvpy(i)+v(i1,k)*txpi
#        else
         pupx(i)=pupx(i)+u(i1,k)*(y11-y33)
         pupy(i)=pupy(i)+u(i1,k)*(x33-x11)
         pvpx(i)=pvpx(i)+v(i1,k)*(y11-y33)
         pvpy(i)=pvpy(i)+v(i1,k)*(x33-x11)
#        endif
       end do
       j=ntve(i)
       i1=nbve(i,j)
       jtmp=nbvt(i,j)
       j1=jtmp+1-(jtmp+1)/4*3
       j2=jtmp+2-(jtmp+2)/4*3
       x11=0.5_sp*(vx(i)+vx(nv(i1,j1)))
       y11=0.5_sp*(vy(i)+vy(nv(i1,j1)))
       x22=xc(i1)
       y22=yc(i1)
       x33=0.5_sp*(vx(i)+vx(nv(i1,j2)))
       y33=0.5_sp*(vy(i)+vy(nv(i1,j2)))

#      if defined (SPHERICAL)
       ty=0.5*(Y11+Y33)
       txpi=(x33-x11)*tpi*cos(deg2rad*TY)
       typi=(y11-y33)*tpi
       pupx(i)=pupx(i)+u(i1,k)*typi
       pupy(i)=pupy(i)+u(i1,k)*txpi
       pvpx(i)=pvpx(i)+v(i1,k)*typi
       pvpy(i)=pvpy(i)+v(i1,k)*txpi
#      else
       pupx(i)=pupx(i)+u(i1,k)*(y11-y33)
       pupy(i)=pupy(i)+u(i1,k)*(x33-x11)
       pvpx(i)=pvpx(i)+v(i1,k)*(y11-y33)
       pvpy(i)=pvpy(i)+v(i1,k)*(x33-x11)
#      endif

       if(isonb(i) /= 0) then
#      if defined (SPHERICAL)
         ty=0.5*(Y11+VY(I))
         txpi=(VX(I)-X11)*tpi*COS(deg2rad*ty)
         typi=(Y11-VY(I))*tpi
         pupx(i)=pupx(i)+u(i1,k)*typi
         pupy(i)=pupy(i)+u(i1,k)*txpi
         pvpx(i)=pvpx(i)+v(i1,k)*typi
         pvpy(i)=pvpy(i)+v(i1,k)*txpi
#        else
         pupx(i)=pupx(i)+u(i1,k)*(y11-vy(i))
         pupy(i)=pupy(i)+u(i1,k)*(vx(i)-x11)
         pvpx(i)=pvpx(i)+v(i1,k)*(y11-vy(i))
         pvpy(i)=pvpy(i)+v(i1,k)*(vx(i)-x11)
#        endif
       end if
       pupx(i)=pupx(i)/art1(i)
       pupy(i)=pupy(i)/art1(i)
       pvpx(i)=pvpx(i)/art1(i)
       pvpy(i)=pvpy(i)/art1(i)
       tmp1=pupx(i)**2+pvpy(i)**2
       tmp2=0.5_sp*(pupy(i)+pvpx(i))**2
       viscoff(i)=sqrt(tmp1+tmp2)*art1(i)
     end do
!     if(k == kbm1) then
!       ah_bottom(1:m) = horcon*(fact*viscoff(1:m) + fm1)
!     end if


     do i=1,ncv_i
       ia=niec(i,1)
       ib=niec(i,2)
       xi=0.5_sp*(xije(i,1)+xije(i,2))
       yi=0.5_sp*(yije(i,1)+yije(i,2))
#      if defined (SPHERICAL)
       ty=0.5_sp*(yi+vy(ia))
       dxa=(xi-vx(ia))*tpi*cos(deg2rad*ty)
       dya=(yi-vy(ia))*tpi
       ty=0.5*(YI+VY(IB))
       DXB=(XI-VX(IB))*tpi*COS(deg2rad*ty)
       DYB=(YI-VY(IB))*tpi
#      else
       dxa=xi-vx(ia)
       dya=yi-vy(ia)
       dxb=xi-vx(ib)
       dyb=yi-vy(ib)
#      endif

!       fij1=f(ia,k)+dxa*pfpx(ia)+dya*pfpy(ia)
!       fij2=f(ib,k)+dxb*pfpx(ib)+dyb*pfpy(ib)
!J. Ge for tracer advection
       IF(BACKWARD_ADVECTION .NEQV. .TRUE.)THEN
         fij1=f(ia,k)+dxa*pfpx(ia)+dya*pfpy(ia)
         fij2=f(ib,k)+dxb*pfpx(ib)+dyb*pfpy(ib)
       ELSE
         IF(BACKWARD_STEP==1)THEN
           fij1=(f0(ia,k)+f(ia,k))*0.5+dxa*pfpx(ia)+dya*pfpy(ia)
           fij2=(f0(ib,k)+f(ib,k))*0.5+dxb*pfpx(ib)+dyb*pfpy(ib)
         ELSEIF(BACKWARD_STEP==2)THEN
           fij1=(f2(ia,k)+f0(ia,k)+f(ia,k))/3.0_SP+dxa*pfpx(ia)+dya*pfpy(ia)
           fij2=(f2(ib,k)+f0(ib,k)+f(ib,k))/3.0_SP+dxb*pfpx(ib)+dyb*pfpy(ib)
         ENDIF
       ENDIF
!J. Ge for tracer advection
       un=uvn(i,k)

!       viscof=horcon*(fact*(viscoff(ia)+viscoff(ib))*0.5_sp + fm1)
        VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))

       txx=0.5_sp*(pfpxd(ia)+pfpxd(ib))*viscof
       tyy=0.5_sp*(pfpyd(ia)+pfpyd(ib))*viscof

       fxx=-dtij(i,k)*txx*dltye(i)
       fyy= dtij(i,k)*tyy*dltxe(i)

# if defined (PLBC)
       fyy=0.0_SP
# endif

       exflux=-un*dtij(i,k)* &
          ((1.0_sp+sign(1.0_sp,un))*fij2+(1.0_sp-sign(1.0_sp,un))*fij1)*0.5_sp+fxx+fyy


       ! --------- new: 14-02-2012 ,K.Lettmann ------
       ! Limit the total flux according to the amount
       ! of mass within the respective neigbour control volumnes

# if !defined (SPHERICAL)
       Da = dz(ia,k)*dt(ia)
       Db = dz(ib,k)*dt(ib)
       Ds = 0.5_sp*( dz(ia,k) + dz(ib,k) )  ! Approximation of water depth for segment in sigma units

       exflux = exflux / Ds ! trafo to old eflux unit (FVCOM 2.7)
       call limit_hor_flux_Scal(f(ia,k),f(ib,k),ia,ib,Da,Db,ntrg(i),deltat,exflux)
       exflux = exflux * Ds ! trafo back to new eflux unit (>FVCOM 3.1.4)
# endif
       ! --------------------------------------------


       xflux(ia,k)=xflux(ia,k)+exflux
       xflux(ib,k)=xflux(ib,k)-exflux

       xflux_adv(ia,k)=xflux_adv(ia,k)+(exflux-fxx-fyy)
       xflux_adv(ib,k)=xflux_adv(ib,k)-(exflux-fxx-fyy)

#      if defined (THIN_DAM)
       IF(K<=KDAM(IA).AND.IS_DAM(IA)==1)THEN
         IF(N_DAM_MATCH(IA,1)==1)THEN
           XFLUX(N_DAM_MATCH(IA,2),K) = XFLUX(N_DAM_MATCH(IA,2),K) + EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IA,2),K) = XFLUX_ADV(N_DAM_MATCH(IA,2),K) +(EXFLUX-FXX-FYY)
         END IF
         IF(N_DAM_MATCH(IA,1)==2)THEN
           XFLUX(N_DAM_MATCH(IA,2),K) = XFLUX(N_DAM_MATCH(IA,2),K) + EXFLUX
           XFLUX(N_DAM_MATCH(IA,3),K) = XFLUX(N_DAM_MATCH(IA,3),K) + EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IA,2),K) = XFLUX_ADV(N_DAM_MATCH(IA,2),K) +(EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IA,3),K) = XFLUX_ADV(N_DAM_MATCH(IA,3),K) +(EXFLUX-FXX-FYY)
         END IF
         IF(N_DAM_MATCH(IA,1)==3)THEN
           XFLUX(N_DAM_MATCH(IA,2),K) = XFLUX(N_DAM_MATCH(IA,2),K) + EXFLUX
           XFLUX(N_DAM_MATCH(IA,3),K) = XFLUX(N_DAM_MATCH(IA,3),K) + EXFLUX
           XFLUX(N_DAM_MATCH(IA,4),K) = XFLUX(N_DAM_MATCH(IA,4),K) + EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IA,2),K) = XFLUX_ADV(N_DAM_MATCH(IA,2),K) +(EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IA,3),K) = XFLUX_ADV(N_DAM_MATCH(IA,3),K) +(EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IA,4),K) = XFLUX_ADV(N_DAM_MATCH(IA,4),K) +(EXFLUX-FXX-FYY)
         END IF
       END IF
       IF(K<=KDAM(IB).AND.IS_DAM(IB)==1)THEN
         IF(N_DAM_MATCH(IB,1)==1)THEN
           XFLUX(N_DAM_MATCH(IB,2),K) = XFLUX(N_DAM_MATCH(IB,2),K) - EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IB,2),K) = XFLUX_ADV(N_DAM_MATCH(IB,2),K) - (EXFLUX-FXX-FYY)
         END IF
         IF(N_DAM_MATCH(IB,1)==2)THEN
           XFLUX(N_DAM_MATCH(IB,2),K) = XFLUX(N_DAM_MATCH(IB,2),K) - EXFLUX
           XFLUX(N_DAM_MATCH(IB,3),K) = XFLUX(N_DAM_MATCH(IB,3),K) - EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IB,2),K) = XFLUX_ADV(N_DAM_MATCH(IB,2),K) - (EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IB,3),K) = XFLUX_ADV(N_DAM_MATCH(IB,3),K) - (EXFLUX-FXX-FYY)
         END IF
         IF(N_DAM_MATCH(IB,1)==3)THEN
           XFLUX(N_DAM_MATCH(IB,2),K) = XFLUX(N_DAM_MATCH(IB,2),K) - EXFLUX
           XFLUX(N_DAM_MATCH(IB,3),K) = XFLUX(N_DAM_MATCH(IB,3),K) - EXFLUX
           XFLUX(N_DAM_MATCH(IB,4),K) = XFLUX(N_DAM_MATCH(IB,4),K) - EXFLUX
           XFLUX_ADV(N_DAM_MATCH(IB,2),K) = XFLUX_ADV(N_DAM_MATCH(IB,2),K) - (EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IB,3),K) = XFLUX_ADV(N_DAM_MATCH(IB,3),K) - (EXFLUX-FXX-FYY)
           XFLUX_ADV(N_DAM_MATCH(IB,4),K) = XFLUX_ADV(N_DAM_MATCH(IB,4),K) - (EXFLUX-FXX-FYY)
         END IF
       END IF
#      endif
     end do
  end do !!sigma loop

!---------------------------------------------------------------------------------
! Accumulate Fluxes at Boundary Nodes
!---------------------------------------------------------------------------------
 
# if defined (MULTIPROCESSOR)
  if(par)call node_match(0,nbn,bn_mlt,bn_loc,bnc,mt,kb,myid,nprocs,xflux,xflux_adv)
# endif

!---------------------------------------------------------------------------------
! Store Advective Fluxes at the Boundary
!---------------------------------------------------------------------------------
  do k=1,kbm1
     if(iobcn > 0) then
       do i=1,iobcn
         i1=i_obc_n(i)
         fflux_obc(i,k)=xflux_adv(i1,k)
       end do
     end if
  end do

!---------------------------------------------------------------------------------
! Calculate Vertical Advection Terms 
!---------------------------------------------------------------------------------

   do i=1,m 
#    if defined (WET_DRY)
     if(iswetn(i)*iswetnt(i) == 1) then
#    endif
#    if defined (THIN_DAM)
     if(IS_DAM(I)==1)then
       wvel(1:kb)=0.0_sp
!YC    call calc_vflux(kbm1,f(i,1:kbm1),wvel(1:kb),vflux)
!J. Ge for tracer advection
       !w_tmp = wvel(1:kb)
       !dz_tmp= dz(i,1:kbm1)
       IF(BACKWARD_ADVECTION .NEQV. .TRUE.)THEN
         ftmp = f(i,1:kbm1)
         !call calc_vflux(kbm1,ftmp,w_tmp,dz_tmp,vflux)
         call calc_vflux(kbm1,f(i,1:kbm1),wvel(1:kb),dz(i,1:kbm1),vflux)
       ELSE
         IF(BACKWARD_STEP==1)THEN
           ftmp=(f0(i,1:kbm1)+f(i,1:kbm1))*0.5
         ELSEIF(BACKWARD_STEP==2)THEN
           ftmp=(f2(i,1:kbm1)+f0(i,1:kbm1)+f(i,1:kbm1))/3.0_SP
         ENDIF
         !call calc_vflux(kbm1,ftmp,w_tmp,dz_tmp,vflux)
         call calc_vflux(kbm1,ftmp,wvel(1:kb),dz(i,1:kbm1),vflux)
       ENDIF
!J. Ge for tracer advection
     else
!YC    call calc_vflux(kbm1,f(i,1:kbm1),wts(i,1:kb),vflux)
!J. Ge for tracer advection
       !w_tmp = wts(i,1:kb)
       !dz_tmp= dz(i,1:kbm1)
       IF(BACKWARD_ADVECTION .NEQV. .TRUE.)THEN
         ftmp = f(i,1:kbm1)
         !call calc_vflux(kbm1,ftmp,w_tmp,dz_tmp,vflux)
         call calc_vflux(kbm1,f(i,1:kbm1),wts(i,1:kb),dz(i,1:kbm1),vflux)
       ELSE
         IF(BACKWARD_STEP==1)THEN
           ftmp=(f0(i,1:kbm1)+f(i,1:kbm1))*0.5
         ELSEIF(BACKWARD_STEP==2)THEN
           ftmp=(f2(i,1:kbm1)+f0(i,1:kbm1)+f(i,1:kbm1))/3.0_SP
         ENDIF
         !call calc_vflux(kbm1,ftmp,w_tmp,dz_tmp,vflux)
         call calc_vflux(kbm1,ftmp,wts(i,1:kb),dz(i,1:kbm1),vflux)
       ENDIF
!J. Ge for tracer advection
     end if
#    else
!J. Ge for tracer advection
       !w_tmp = wts(i,1:kb)
       !dz_tmp= dz(i,1:kbm1)
       IF(BACKWARD_ADVECTION .NEQV. .TRUE.)THEN
         !ftmp = f(i,1:kbm1)
         !call calc_vflux(kbm1,ftmp,w_tmp,dz_tmp,vflux)
         call calc_vflux(kbm1,f(i,1:kbm1),wts(i,1:kb),dz(i,1:kbm1),vflux)
       ELSE
         IF(BACKWARD_STEP==1)THEN
           ftmp=(f0(i,1:kbm1)+f(i,1:kbm1))*0.5
         ELSEIF(BACKWARD_STEP==2)THEN
           ftmp=(f2(i,1:kbm1)+f0(i,1:kbm1)+f(i,1:kbm1))/3.0_SP
         ENDIF
         !call calc_vflux(kbm1,ftmp,w_tmp,dz_tmp,vflux)
         call calc_vflux(kbm1,ftmp,wts(i,1:kb),dz(i,1:kbm1),vflux)
       ENDIF
!J. Ge for tracer advection
!!!      vflux = 0.0_SP
#    endif

     do k=1,kbm1
       if(isonb(i) == 2) then
         xflux(i,k)= vflux(k)*art1(i)   !/dz(i,k)
!JQI         xflux(i,k)= (vflux(k)-vflux(k+1))*art1(i)/dz(i,k)
       else
!JQI         xflux(i,k)=xflux(i,k)+ (vflux(k)-vflux(k+1))*art1(i)/dz(i,k)
         xflux(i,k)=xflux(i,k)+ vflux(k)*art1(i)    !/dz(i,k)
       end if
#    if defined (THIN_DAM)
       IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
         tmpflx = (vflux(k)-vflux(k+1))*art1(i)/dz(i,k)
         IF(N_DAM_MATCH(I,1)==1)THEN
            XFLUX(N_DAM_MATCH(I,2),K) = XFLUX(N_DAM_MATCH(I,2),K)+tmpflx
         END IF
         IF(N_DAM_MATCH(I,1)==2)THEN
            XFLUX(N_DAM_MATCH(I,2),K) = XFLUX(N_DAM_MATCH(I,2),K)+tmpflx
            XFLUX(N_DAM_MATCH(I,3),K) = XFLUX(N_DAM_MATCH(I,3),K)+tmpflx
         END IF
         IF(N_DAM_MATCH(I,1)==3)THEN
            XFLUX(N_DAM_MATCH(I,2),K) = XFLUX(N_DAM_MATCH(I,2),K)+tmpflx
            XFLUX(N_DAM_MATCH(I,3),K) = XFLUX(N_DAM_MATCH(I,3),K)+tmpflx
            XFLUX(N_DAM_MATCH(I,4),K) = XFLUX(N_DAM_MATCH(I,4),K)+tmpflx
         END IF                
       END IF
#    endif
     end do
#    if defined (WET_DRY)
     end if
#    endif
   end do

!-------------------------------------------------------
!Point Source                                      
!-------------------------------------------------------
  if(source)then  !!user specified

  if(RIVER_TS_SETTING == 'calculated') then
    if(RIVER_INFLOW_LOCATION == 'node') then
        do j=1,numqbc
          jj=inodeq(j)
          fpoint=fdis(j)
          do k=1,kbm1
            xflux(jj,k)=xflux(jj,k) - qdis(j)*vqdist(j,k)*fpoint !/dz(jj,k)
          end do
        end do
    else if(RIVER_INFLOW_LOCATION == 'edge') then
      write(*,*)'scalar advection not setup for "edge" point source'
      stop
    end if
  end if

  else

  if(RIVER_TS_SETTING == 'calculated')then   
    if(RIVER_INFLOW_LOCATION == 'node') then
        do j=1,numqbc
          jj=inodeq(j)
          do k=1,kbm1
            fpoint = f(jj,k)
!J. Ge for tracer advection
            IF(BACKWARD_ADVECTION .NEQV. .TRUE.)THEN
              fpoint = f(jj,k)
            ELSE
              IF(BACKWARD_STEP==1)THEN
                fpoint = (f0(jj,k)+f(jj,k))*0.5
              ELSEIF(BACKWARD_STEP==2)THEN
                fpoint = (f2(jj,k)+f0(jj,k)+f(jj,k))/3.0_SP
              ENDIF
            ENDIF
!J. Ge for tracer advection
            xflux(jj,k)=xflux(jj,k) - qdis(j)*vqdist(j,k)*fpoint! /dz(jj,k)
          end do
        end do
    else if(RIVER_INFLOW_LOCATION == 'edge') then
      write(*,*)'scalar advection not setup for "edge" point source'
      stop
    end if
  end if

  endif
!------------------------------------------------------------------------
!Update Scalar Quantity
!------------------------------------------------------------------------

  do i=1,m
#   if defined (WET_DRY)
    if(iswetn(i)*iswetnt(i) == 1 )then
#   endif
    do k=1,kbm1
#   if !defined (THIN_DAM)
      !fn(i,k)=(f(i,k)-xflux(i,k)/art1(i)*(deltat/dt(i)))*(dt(i)/dtfa(i))
      fn(i,k)=(f(i,k)-xflux(i,k)/art1(i)*(deltat/(dt(i)*dz(i,k))))*(dt(i)/dtfa(i)) 
#   else
      IF(IS_DAM(I)==1.AND.K<=KDAM(I))THEN
        fn(i,k)=(f(i,k)-xflux(i,k)/(ART1(I)&
        &+SUM(ART1(N_DAM_MATCH(I,2:1+N_DAM_MATCH(I,1)))))*(deltat/dt(i)))*(dt(i)/dtfa(i))
      ELSE
        !fn(i,k)=(f(i,k)-xflux(i,k)/art1(i)*(deltat/dt(i)))*(dt(i)/dtfa(i))
        fn(i,k)=(f(i,k)-xflux(i,k)/art1(i)*(deltat/(dt(i)*dz(i,k))))*(dt(i)/dtfa(i)) 
      END IF
#   endif
    end do
#   if defined (WET_DRY)
    else
    do k=1,kbm1
      fn(i,k)=f(i,k)
    end do
    end if
#   endif
  end do

  return
  End Subroutine Adv_Scal
!==============================================================================|

!==============================================================================|
! Vertical Diffusion of Scalar                                                 |
!==============================================================================|
  Subroutine Vdif_Scal(f,deltat)

  use mTridiagonal_scal
  use all_vars 
# if defined (THIN_DAM)
  use mod_dam,only : NODE_DAM1_N,NODE_DAM2_N,NODE_DAM3_N, &
                   &I_NODE_DAM1_N,I_NODE_DAM2_N,I_NODE_DAM3_N, &
                   &kdam
# endif

  Implicit None 
  Real(sp), intent(inout) :: f(0:mt,kb)
  Real(sp), intent(in   ) :: deltat
  !--local--------------------
  integer  :: i,k,ll
  real(sp) :: dsqrd,dfdz,visb
  real(sp) :: fsol(0:kb)

# if defined (THIN_DAM)
  real(sp) :: ftmp,stmp
# endif

  call init_tridiagonal_scal(kb)

  Do i=1,m
     dsqrd = d(i)*d(i)

    !----------------------------------------------------------------
    !  Set up Diagonals of Matrix (lower=au,diag=bu,upper=cu)
    !----------------------------------------------------------------
    

    !Surface
    au(1) = 0.0
    cu(1)=      - deltat*(kh(i,2)+umol)/(dzz(i,1)*dz(i,1)*dsqrd)
    bu(1)=  1.0 - cu(1) 

    !Interior
    do k=2,kbm1-1
      au(k) =     - deltat*(kh(i,k  )+umol)/(dzz(i,k-1)*dz(i,k)*dsqrd)
      cu(k) =     - deltat*(kh(i,k+1)+umol)/(dzz(i,k  )*dz(i,k)*dsqrd)
      bu(k) = 1.0 - cu(k) - au(k) 
    end do

    !Bottom
     au(kbm1) =     - deltat*(kh(i,kbm1)+umol)/(dzz(i,kbm1-1)*dz(i,kbm1)*dsqrd)
     cu(kbm1) = 0.0
     bu(kbm1) = 1.0 - au(kbm1) 

    !----------------------------------------------------------------
    ! Set up RHS forcing vector and boundary conditions 
    !----------------------------------------------------------------
    do k=1,kbm1
      du(k) = f(i,k)
    end do

    !Free Surface: No flux

    !Bottom: No flux
      

    !----------------------------------------------------------------
    ! Solve 
    !----------------------------------------------------------------

     call tridiagonal_scal(kb,1,kbm1,fsol)
    
     !Transfer
     f(i,1:kbm1) = fsol(1:kbm1)

  End Do

#  if defined (THIN_DAM)
   DO K=1,KBM1
     DO I=1,NODE_DAM1_N
       IF(K<=KDAM(I_NODE_DAM1_N(I,1)).AND.K<=KDAM(I_NODE_DAM1_N(I,2)) )THEN
          FTMP=F(I_NODE_DAM1_N(I,1),K)*ART1(I_NODE_DAM1_N(I,1)) &
            & +F(I_NODE_DAM1_N(I,2),K)*ART1(I_NODE_DAM1_N(I,2))
          STMP=ART1(I_NODE_DAM1_N(I,1))+ART1(I_NODE_DAM1_N(I,2))
          F(I_NODE_DAM1_N(I,1),K)=FTMP/STMP
          F(I_NODE_DAM1_N(I,2),K)=FTMP/STMP
       END IF
     END DO

     DO I=1,NODE_DAM2_N
       IF(K<=KDAM(I_NODE_DAM2_N(I,1)).AND.K<=KDAM(I_NODE_DAM2_N(I,2)) &
          & .AND.K<=KDAM(I_NODE_DAM2_N(I,2)) )THEN
          FTMP= F(I_NODE_DAM2_N(I,1),K)*ART1(I_NODE_DAM2_N(I,1)) &
           &   +F(I_NODE_DAM2_N(I,2),K)*ART1(I_NODE_DAM2_N(I,2)) &
           &   +F(I_NODE_DAM2_N(I,3),K)*ART1(I_NODE_DAM2_N(I,3)) 
          STMP=ART1(I_NODE_DAM2_N(I,1))+ART1(I_NODE_DAM2_N(I,2)) &
           &   +ART1(I_NODE_DAM2_N(I,3))
          F(I_NODE_DAM2_N(I,1),K)=FTMP/STMP
          F(I_NODE_DAM2_N(I,2),K)=FTMP/STMP
          F(I_NODE_DAM2_N(I,3),K)=FTMP/STMP
       END IF
     END DO

     DO I=1,NODE_DAM3_N
       IF(K<=KDAM(I_NODE_DAM3_N(I,1)).AND.K<=KDAM(I_NODE_DAM3_N(I,2)) &
   & .AND.K<=KDAM(I_NODE_DAM3_N(I,3)).AND.K<=KDAM(I_NODE_DAM3_N(I,4)) )THEN
          FTMP =F(I_NODE_DAM3_N(I,1),K)*ART1(I_NODE_DAM3_N(I,1)) &
           &   +F(I_NODE_DAM3_N(I,2),K)*ART1(I_NODE_DAM3_N(I,2)) &
           &   +F(I_NODE_DAM3_N(I,3),K)*ART1(I_NODE_DAM3_N(I,3)) &
           &   +F(I_NODE_DAM3_N(I,4),K)*ART1(I_NODE_DAM3_N(I,4))  
          STMP =ART1(I_NODE_DAM3_N(I,1)) + ART1(I_NODE_DAM3_N(I,2)) &
           &  + ART1(I_NODE_DAM3_N(I,3)) + ART1(I_NODE_DAM3_N(I,4))
          F(I_NODE_DAM3_N(I,1),K)=FTMP/STMP
          F(I_NODE_DAM3_N(I,2),K)=FTMP/STMP
          F(I_NODE_DAM3_N(I,3),K)=FTMP/STMP
          F(I_NODE_DAM3_N(I,4),K)=FTMP/STMP
       END IF
     END DO
   END DO
#  endif


  End Subroutine Vdif_Scal


!==============================================================================|
! Set Point Source Conditions for Scalar Function                              |
!==============================================================================|

  Subroutine Bcond_Scal_PTsource(f,fn,fdis)

!------------------------------------------------------------------------------|
  use all_vars
  use bcs
  use mod_obcs
  implicit none
  real(sp), intent(in ), dimension(0:mt,kb)      :: f 
  real(sp), intent(out), dimension(0:mt,kb)      :: fn
  real(sp), intent(in ), dimension(numqbc )      :: fdis
!--local-------------------------------------------
  integer  :: i,j,k,j1,j11,j22
!------------------------------------------------------------------------------|


!--------------------------------------------
! Set Source Terms
!--------------------------------------------
  if(RIVER_TS_SETTING == 'specified') then
    if(numqbc > 0) then
      if(RIVER_INFLOW_LOCATION == 'node') then
        do i=1,numqbc
          j11=inodeq(i)
          do k=1,kbm1
            fn(j11,k)=fdis(i)
          end do
        end do
      else if(RIVER_INFLOW_LOCATION == 'edge') then
        do i=1,numqbc
          j11=n_icellq(i,1)
          j22=n_icellq(i,2)
          do k=1,kbm1
            fn(j11,k)=fdis(i)
            fn(j22,k)=fdis(i)
          end do
        end do
      end if
    end if
  end if

  return
  End Subroutine Bcond_Scal_PTSource 
!==============================================================================|
!==============================================================================|

!==============================================================================|
!   Set Boundary Conditions for Scalar Function on Open Boundary               |
!==============================================================================|

  Subroutine Bcond_Scal_OBC(f,fn,fflux_obc,f_obc,deltat,alpha_nudge)

!------------------------------------------------------------------------------|
  use all_vars
  use bcs
  use mod_obcs
  implicit none
  real(sp), intent(in   ), dimension(0:mt,kb)      :: f 
  real(sp), intent(inout), dimension(0:mt,kb)      :: fn
  real(sp), intent(in   ), dimension(iobcn+1,kbm1) :: fflux_obc 
  real(sp), intent(in   ), dimension(iobcn       ) :: f_obc 
  real(sp), intent(in   )                          :: deltat
  real(sp), intent(in   )                          :: alpha_nudge 
!--local-------------------------------------------
  real(sp) :: f2d,f2d_next,f2d_obc,xflux2d,tmp
  integer  :: i,j,k,j1,j11,j22
!------------------------------------------------------------------------------|
       
!--------------------------------------------
! Set Scalar Value on Open Boundary
!--------------------------------------------
  if(iobcn > 0) then
    do i=1,iobcn
      j=i_obc_n(i)
      j1=next_obc(i)
      f2d=0.0_sp
      f2d_next=0.0_sp
      xflux2d=0.0_sp
      do k=1,kbm1
        f2d=f2d+f(j,k)*dz(j,k)
        f2d_next=f2d_next+fn(j1,k)*dz(j1,k)
        xflux2d=xflux2d+fflux_obc(i,k)*dz(j,k)
      end do
  
      if(uard_obcn(i) > 0.0_sp) then
        tmp=xflux2d+f2d*uard_obcn(i)
        f2d_obc=(f2d*dt(j)-tmp*deltat/art1(j))/d(j)
        do k=1,kbm1
          fn(j,k)=fn(j1,k) !f2d_obc+(fn(j1,k)-f2d_next)
        end do
      else
        do k=1,kbm1
          fn(j,k) = f(j,k)-alpha_nudge*(f(j,k)-f_obc(i))
        end do
      end if
    end do
  endif

  return
  End Subroutine Bcond_Scal_OBC 
!==============================================================================|
!==============================================================================|

  Subroutine fct_sed(f,fn)
  !==============================================================================|
  USE ALL_VARS
  USE MOD_UTILS
  USE BCS
  USE MOD_OBCS
  IMPLICIT NONE
  real(sp), intent(inout), dimension(0:mt,kb)      :: fn
  real(sp), intent(in), dimension(0:mt,kb)      :: f
  REAL(SP):: SMAX,SMIN
  INTEGER :: I,J,K,K1
  !==============================================================================|
  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: fct_sed"

  nodes: DO I=1,M

     ! SKIP OPEN BOUNDARY NODES
     IF(IOBCN > 0)THEN
        DO J=1,IOBCN
           IF(I == I_OBC_N(J)) CYCLE nodes
        END DO
     END IF

     ! SKIP RIVER INFLOW POINTS
     IF(NUMQBC > 0)THEN
        DO J=1,NUMQBC
           IF(RIVER_INFLOW_LOCATION == 'node')THEN
              IF(I == INODEQ(J)) CYCLE nodes
           END IF
           IF(RIVER_INFLOW_LOCATION == 'edge')THEN
              IF(I == N_ICELLQ(J,1) .OR. I == N_ICELLQ(J,2)) CYCLE nodes
           END IF
        END DO
     END IF

     ! SKIP GROUND WATER INFLOW POINTS
     IF(BFWDIS(I) .GT. 0.0_SP .and. GROUNDWATER_SALT_ON) CYCLE nodes

     K1 = 1
     IF(PRECIPITATION_ON) K1 = 2
!     DO K=1,KBM1
     DO K=K1,KBM1
        SMAX = MAXVAL(f(NBSN(I,1:NTSN(I)),K))
        SMIN = MINVAL(f(NBSN(I,1:NTSN(I)),K))

        IF(K == 1)THEN
           SMAX = MAX(SMAX,(f(I,K)*DZ(I,K+1)+f(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
           SMIN = MIN(SMIN,(f(I,K)*DZ(I,K+1)+f(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
        ELSE IF(K == KBM1)THEN
           SMAX = MAX(SMAX,(f(I,K)*DZ(I,K-1)+f(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
           SMIN = MIN(SMIN,(f(I,K)*DZ(I,K-1)+f(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
        ELSE
           SMAX = MAX(SMAX,(f(I,K)*DZ(I,K-1)+f(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (f(I,K)*DZ(I,K+1)+f(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
           SMIN = MIN(SMIN,(f(I,K)*DZ(I,K-1)+f(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (f(I,K)*DZ(I,K+1)+f(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
        END IF

        IF(SMIN-fn(I,K) > 0.0_SP)fn(I,K) = SMIN
        IF(fn(I,K)-SMAX > 0.0_SP)fn(I,K) = SMAX

     END DO
  END DO nodes

  WHERE(fn < 0.0_SP)fn=0.0_SP
  
  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: fct_sed"
  End Subroutine fct_sed

!==========================================================================
! Calculate Fluxes for Vertical Advection Equation                            
! n: number of cells
! c: scalar variable (1:n)
! w: velocity field at cell interfaces (1:n+1)
! note:  we dont use face normals to construct inflow/outflow
!        thus we add dissipation term instead of subtracting because 
!        positive velocity is up while computational coordinates increase
!        down towards bottom.  
!==========================================================================
  Subroutine Calc_VFlux(n,c,w,zk,flux) 
  use mod_prec
  use mod_utils, only : LIMLED2
  implicit none
  integer , intent(in ) :: n
  real(sp), intent(in ) ::  c(n)
  real(sp), intent(in ) ::  w(n+1) 
  real(sp), intent(in ) ::  zk(n)
  real(sp), intent(out) ::  flux(n+1)
  real(sp) :: conv(n+1),diss(n+1)
  real(sp) :: cin(0:n+1)
  real(sp) :: sl_h(0:n+1)
  real(sp) :: dis4,sl_u,sl_f
  integer  :: i

  !transfer to working array
  cin(1:n) = c(1:n)
  sl_h(1:n) = zk(1:n)

  !surface bcs (no flux)
  cin(0)  =  -cin(1) 
!  cin(-1) =  -cin(2)
  sl_h(0) = zk(1)
  
  !bottom bcs (no flux)
  cin(n+1) = -cin(n) 
!  cin(n+2) = -cin(n-1)
  sl_h(n+1) = zk(n)

  !flux computation
  do i=2,n
    dis4    = .5*abs(w(i))
    conv(i) = w(i)*(cin(i)+cin(i-1))*0.5_SP
    sl_u = 2.0_SP*(cin(i)-cin(i+1))/(sl_h(i)+sl_h(i+1))
    sl_f = 2.0_SP*(cin(i-2)-cin(i-1))/(sl_h(i-2)+sl_h(i-1))
    diss(i) = dis4*(cin(i-1)-cin(i)-0.5_SP*LIMLED2(sl_u,sl_f,2.0_SP)*(sl_h(i-1)+sl_h(i))) 
!!    flux(i) = conv(i)+diss(i)
  end do
  conv(1) = 0.0_SP
  diss(1) = 0.0_SP
  conv(n+1) = 0.0_SP
  diss(n+1) = 0.0_SP
  
  do i=1,n
    flux(i) = conv(i)-conv(i+1)+diss(i+1)-diss(i)
  end do

  End Subroutine Calc_VFlux
  
!==========================================================================
! Calculate LED Limiter L(u,v)  
!==========================================================================
  Function Lim(a,b)
  use mod_prec
  real(sp) lim,a,b
  real(sp) q,R
  real(sp) eps
  eps = epsilon(eps)
  
  q = 0. !1st order
  q = 1. !minmod
  q = 2. !van leer

  R = abs(   (a-b)/(abs(a)+abs(b)+eps) )**q
  lim = .5*(1-R)*(a+b)

  End Function Lim


! =============================================================================
! K.lettmann, Feb 2012, lettmann@icbm.de
!
! Methods for flux limiting for problems witht the diffusion step
! to omit negative concentrations
! ==============================================================================
!
!==============================================================================|
! Horizontal flux Limitation For Scalar (f)                  |
! 
! K.lettmann, Feb 2012, lettmann@icbm.de
!==============================================================================|
  Subroutine limit_hor_flux_Scal(fka,fkb,ia,ib,Da,Db,ele_num,delt,exflux)
!
! This subroutine limits the total flux about a control volume segment
! in case there is not enough tracer in the source element, 
! the flux is coming from.
!------------------------------------------------------------------------------|

  use all_vars
  use lims, only: m,mt,n,nt,kbm1,kb
  use bcs
  use mod_obcs
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  endif

  implicit none
  real(sp), intent(in )                  :: fka        ! scalar quantity on k-sigma level in cv a
  real(sp), intent(in )                  :: fkb        ! scalar quantity on k-sigma level in cv b
  integer , intent(in )                  :: ele_num    ! The element/triangle we are working
  integer , intent(in )                  :: ia,ib      ! Index of the left and right node to the SCV side S
  real(sp), intent(in )                  :: delt       ! time step for diffusion equation
  real(sp), intent(in )                  :: Da         ! water depth  of sigma layer [m] on side a of interface
  real(sp), intent(in )                  :: Db         ! water depth  of sigma layer [m] on side b of interface
  real(sp), intent(inout)                :: exflux     ! total flux over side S (diffusion + advection) []

  !----------------local--------------------------------------
  integer  :: i_source
  real(sp) :: x_ce,xa,xb,y_ce,ya,yb !coordinates of triangle nodes
  real(sp) :: f_source              !concentrations in source element
  real(sp) :: A_cv                  !area of part of control volume
  real(sp) :: delh                  !hight of sigma-layer
  real(sp) :: D_source              !water depth of sigma layer [m] at source node
  real(sp) :: amount_source         !amount of tracer in source node control volume
  real(sp) :: flux_max              !maximal possible outflux of tracer from source node control volume
# if defined (SPHERICAL)
  real(sp) :: side1,side2,side3     ! length of triangle sides for spherical triangle
# endif
!------------------------------------------------------------------------------!

  ! ----- the coordinates of the corner points of the double control volume ---
  x_ce = xc(ele_num) ;   y_ce = yc(ele_num)  ! center node
  xa = vx(ia) ;   ya = vy(ia)                ! left node
  xb = vx(ib) ;   yb = vy(ib)                ! right node
  !print*,'center node: ',x_ce,y_ce
  ! ---------------------------------------------------------------------------


  ! -------- calculate the area of the tracer control volume -------------------
  ! it is only the small part of the total control volume / not art1
#  if defined (SPHERICAL)
     ! see routine cell_area for an example
     ! how to calculate the area of triangles on a sphere
     CALL ARC(x_ce,y_ce,xa,ya,side1)  
     CALL ARC(xa,ya,xb,yb,side2)
     CALL ARC(xb,yb,xa,ya,side3)

     ! calculate area of triangle from triangle sides
     CALL AREA(side1,side2,side3,A_cv)
     A_cv = 0.5_SP * abs(A_cv)
#  else
   ! area for cartesian triangle on a plane
   A_cv = 0.5* abs( 0.5*( (xb-xa)*(y_ce-ya) - (x_ce-xa)*(yb-ya) ) )
#  endif
  ! ----------------------------------------------------------------------------

  ! ------- find the source element, where the flux comes from -----------------
  if (exflux >= 0.0) then 
       i_source = ia
       f_source = fka
       D_source = Da
  else
       i_source = ib
       f_source = fkb
       D_source = Db
  end if
  !print*,'source node:',i_source
  ! -----------------------------------------------------------------------------

  ! ------- calculate the maximal possible outflow from source element ---------
  ! it cannot flow out more tracer than being in the element
  !print*,'f_source',f_source
  !print*,'D:',Dn(ia),Dn(ib)

  amount_source = f_source * A_cv
  flux_max = amount_source*D_source/delt

  !print*,'flux_max:',flux_max
  ! -----------------------------------------------------------------------------


  ! ------ limit the exflux if necessary ----------------------------------------
  if (exflux >= 0.0) then
    exflux = min(flux_max,exflux)
  else
    exflux = max(-flux_max,exflux)
  end if
  ! -----------------------------------------------------------------------------


  return
  End Subroutine limit_hor_flux_Scal
!==============================================================================|

End Module Scalar
