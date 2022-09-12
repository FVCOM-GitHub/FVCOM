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
module mod_infovar
  use mod_prec

  type infovar
    character(len=20) :: var   
    integer  :: vtype   !=1, node, =2, element
    real(sp) :: val
    real(sp) :: bnd
    real(sp) :: xpos 
    real(sp) :: ypos 
    integer  :: pid
    integer  :: pt
    integer  :: lay
  end type infovar

  contains

  subroutine set_infovar(v,vname,vtype,vval,vbound,xpos,ypos,vpid,vpt,vlay)
    type(infovar), intent(inout) :: v
    character(len=20)    :: vname
    integer,  intent(in) :: vtype
    real(sp), intent(in) :: vval
    real(sp), intent(in) :: vbound
    real(sp), intent(in) :: xpos 
    real(sp), intent(in) :: ypos 
    integer,  intent(in) :: vpid
    integer,  intent(in) :: vpt
    integer,  intent(in) :: vlay 

    v%var   = vname
    v%vtype = vtype
    v%val   = vval
    v%bnd   = vbound
    v%xpos  = xpos
    v%ypos  = ypos 
    v%pt    = vpt
    v%pid   = vpid
    v%lay   = vlay

  end subroutine set_infovar

  subroutine print_infovar(v,iunit) 
    type(infovar), intent(inout) :: v
    integer, intent(in) :: iunit

    write(iunit,*)'variable: ',trim(v%var)
    if(v%vtype==1)then
      write(iunit,*)'type    :  node' 
    else
      write(iunit,*)'type    :  elem' 
    endif
    write(iunit,*)'value   : ',v%val
    write(iunit,*)'bound   : ',v%bnd
    write(iunit,*)'x-pos   : ',v%xpos
    write(iunit,*)'y-pos   : ',v%ypos
    if(v%vtype==1)then
      write(iunit,*)'vertex  :  ',v%pt 
    else
      write(iunit,*)'cell    :  ',v%pt 
    endif
    write(iunit,*)'prod id : ',v%pid
    write(iunit,*)'layer   : ',v%lay

  end subroutine print_infovar


end module mod_infovar
  
module mod_boundschk
  use mod_prec
  use mod_infovar
  implicit none

  type(infovar), allocatable :: vlist(:)
  
  contains

  !==============================================================================|
  subroutine setup_boundschk
    use all_vars
    integer :: nvars

    !report the bounds
    if(boundschk_on)then
      write(ipt,*)'bounds checking  :  on'
      write(ipt,*)'checking interval: ',chk_interval
      write(ipt,*)'veloc_mag_max    : ',veloc_mag_max
      write(ipt,*)'zeta_mag_max     : ',zeta_mag_max
      write(ipt,*)'temp_max         : ',temp_max
      write(ipt,*)'temp_min         : ',temp_min
      write(ipt,*)'salt_max         : ',salt_max
      write(ipt,*)'salt_min         : ',salt_min
#     if defined (WAVE_CURRENT_INTERACTION)
      write(ipt,*)'hs_max           : ',hs_max   
#     endif

      !allocate vlist
      nvars = 7
      allocate(vlist(nvars))

    else
      allocate(vlist(1))
    endif

  end subroutine setup_boundschk


  !==============================================================================|
  !==============================================================================|
  subroutine boundschk
    use all_vars
    use mod_utils, only : pstop
    use mod_wd
    use mod_ncdio, only : archive
    use mod_par,   only : egid,ngid
#   if defined(WAVE_CURRENT_INTERACTION)
    use vars_wave, ONLY : hsc1
#   endif
    implicit none
    type(infovar), allocatable :: vlist_global(:)
    integer  :: icnt,pt,lay,i,j,k,nviolations,printproc,ierr
    integer  :: iiint 
    real(sp) :: val
    integer, allocatable :: violations(:),tmp(:)
    character(len=20) :: vname

    
    iiint=iint
    if(.not.boundschk_on) return
    if(mod(iiint,chk_interval)/= 0) return
!    if(.not.boundschk_on .or. (mod(iint,chk_interval)/= 0)) return
    allocate(tmp(nprocs))
    tmp = 0

    !------------------------------------------------------
    !Check if variables are in user-defined bounds
    !------------------------------------------------------
    icnt = 0 

    ! vert-averaged x-velocity
    val = maxval(abs(ua(1:n)))
    if(val > veloc_mag_max)then
      icnt = icnt + 1
      pt = 1
      vname = 'vert-averaged u'
      ualoop: do i=1,n
        if(abs(ua(i)) > veloc_mag_max)then 
          pt = i
# if defined(MULTIPROCESSOR)
         if(par)then
           pt = egid(pt) 
         endif
# endif
          val = abs(ua(i))
          exit ualoop
        endif
      end do ualoop 
      call set_infovar(vlist(icnt),vname,2,val,veloc_mag_max,xc(pt)+vxmin,yc(pt)+vymin,myid,pt,1)
    endif

    ! vert-averaged y-velocity
    val = maxval(abs(va(1:n)))
    if(val > veloc_mag_max)then      
      icnt = icnt + 1
      pt = 1
      vname = 'vert-averaged v'
      valoop: do i=1,n
        if(abs(va(i)) > veloc_mag_max)then
          pt = i
# if defined(MULTIPROCESSOR)
         if(par)then
           pt = egid(pt) 
         endif
# endif
          val = abs(va(i))
          exit valoop
        endif
      end do valoop
      call set_infovar(vlist(icnt),vname,2,val,veloc_mag_max,xc(pt)+vxmin,yc(pt)+vymin,myid,pt,1)
    endif

    ! x-velocity
    val = maxval(abs(u(1:n,1:kbm1)))
    if(val > veloc_mag_max)then
      icnt = icnt + 1
      pt = 1
      vname = 'u'
      uloop: do k=1,kbm1
      do i=1,n
        if(abs(u(i,k)) > veloc_mag_max)then
          pt = i
# if defined(MULTIPROCESSOR)
         if(par)then
           pt = egid(pt) 
         endif
# endif
          lay = k
          val = abs(u(i,k))
          exit uloop
        endif
      end do 
      end do uloop
      call set_infovar(vlist(icnt),vname,2,val,veloc_mag_max,xc(pt)+vxmin,yc(pt)+vymin,myid,pt,lay)
    endif

    !hydrographic vars
    if(.not. barotropic)then

    ! salinity - max 
    val = maxval(s1(1:m,1:kbm1))
    if(val > salt_max)then
      icnt = icnt + 1
      pt = 1
      vname = 'salinity'
      smaxloop: do k=1,kbm1
      do i=1,m
        if(s1(i,k) > salt_max)then
          pt = i
# if defined(MULTIPROCESSOR)
         if(par)then
           pt = ngid(pt) 
         endif
# endif
          lay = k
          val = s1(i,k)
          exit smaxloop
        endif
      end do
      end do smaxloop
      call set_infovar(vlist(icnt),vname,1,val,salt_max,vx(pt)+vxmin,vy(pt)+vymin,myid,pt,lay)
    endif

    ! salinity - min 
    val = minval(s1(1:m,1:kbm1))
    if(val < salt_min)then
      icnt = icnt + 1
      pt = 1
      vname = 'salinity'
      sminloop: do k=1,kbm1
      do i=1,m
        if(s1(i,k) < salt_min)then
          pt = i
# if defined(MULTIPROCESSOR)
         if(par)then
           pt = ngid(pt) 
         endif
# endif
          lay = k
          val = s1(i,k)
          exit sminloop
        endif
      end do
      end do sminloop
      call set_infovar(vlist(icnt),vname,1,val,salt_min,vx(pt)+vxmin,vy(pt)+vymin,myid,pt,lay)
    endif

    ! temp - max 
    val = maxval(t1(1:m,1:kbm1))
    if(val > temp_max)then
      icnt = icnt + 1
      pt = 1
      vname = 'temperature'
      tmaxloop: do k=1,kbm1
      do i=1,m
        if(t1(i,k) > temp_max)then
          pt = i
# if defined(MULTIPROCESSOR)
         if(par)then
           pt = ngid(pt) 
         endif
# endif
          lay = k
          val = t1(i,k)
          exit tmaxloop
        endif
      end do
      end do tmaxloop
      call set_infovar(vlist(icnt),vname,1,val,temp_max,vx(pt)+vxmin,vy(pt)+vymin,myid,pt,lay)
    endif

    ! temperature - min 
    val = minval(t1(1:m,1:kbm1))
    if(val < temp_min)then
      icnt = icnt + 1
      pt = 1
      vname = 'temperature'
      tminloop: do k=1,kbm1
      do i=1,m
        if(t1(i,k) < temp_min)then
          pt = i
          lay = k
          val = t1(i,k)
          exit tminloop
        endif
      end do
      end do tminloop
      call set_infovar(vlist(icnt),vname,1,val,temp_min,vx(pt)+vxmin,vy(pt)+vymin,myid,pt,lay)
    endif


    endif !.not. barotropic

#   if defined (WAVE_CURRENT_INTERACTION)
    ! significant wave height - max 
    val = maxval(hsc1(1:m))
    if(val > hs_max)then
      icnt = icnt + 1
      pt = 1
      vname = 'signifianct_wave_height'
      hsmaxloop: do i=1,m    
        if(hsc1(i) > hs_max)then
          pt = i
          val = hsc1(i)
          exit hsmaxloop
        endif
      end do hsmaxloop
      call set_infovar(vlist(icnt),vname,1,val,hs_max,vx(pt)+vxmin,vy(pt)+vymin,myid,pt,1)
    endif
#   endif




    !-----------------------------------------------------------
    !Collect number of violations from each proc into vector tmp
    !-----------------------------------------------------------
    nviolations = icnt
    tmp(1) = icnt 

#   if defined (MULTIPROCESSOR)
    IF(PAR)THEN
     call mpi_allgather(nviolations,1,mpi_integer,tmp,1,mpi_integer,mpi_fvcom_group,ierr)
    ENDIF
#   endif

    !check total violations and return if none found
    if(sum(tmp) == 0)then
      deallocate(tmp)
      return 
    endif

    !look at violations, just pick lowest procid with violation to print violation
    printproc = 1
    do i=1,nprocs
      if(tmp(i) > 0)then
        printproc = i
        exit
      endif
    end do

    !write violations to screen
    if(printproc==myid)then
      write(ipt,*)'WARNING: Variable(s) have exceeded user-defined thresholds'
      do i=1,nviolations
        call print_infovar(vlist(i),ipt)
      end do 
      write(ipt,*)'ARCHIVING FRAME AND HALTING'
    endif

    !dump frame to netcdf file
    force_archive = .true.
    call archive

    !shutdown
#   if defined(MULTIPROCESSOR)
    call mpi_barrier(MPI_FVCOM_GROUP,ierr)
#   endif
    call pstop

    deallocate(tmp)    

  end subroutine boundschk

end module mod_boundschk
