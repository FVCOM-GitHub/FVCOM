!=======================================================================
! FVCOM Sediment Initialization Routine (User Defined)
!   Configuration:  default    
!=======================================================================
  Subroutine Init_Sed  
# if defined (SEDIMENT)
  USE MOD_UTILS
  USE MOD_NCTOOLS
  USE MOD_INPUT
  USE ALL_VARS
  USE EQS_OF_STATE
  

# if defined (ORIG_SED)
  Use Mod_Prec 
  Use Mod_Sed
  Use Lims, only: m,kbm1,kb
  implicit none 
  integer :: i,k,ised
  real(sp) :: bed_thickness
  real(sp), allocatable :: tmp(:)

   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCDIM),  POINTER :: DIM
   LOGICAL :: FOUND

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Init_Sed" 

  if(SEDIMENT_PARAMETER_TYPE/=UNIFORM)then
    do i=1,nsed
      call LOAD_SEDIMENT_PARAMETER(sed(i)%t_cd,sed(i)%t_ce,sed(i)%rate)
    end do
  else
    do i=1,nsed
      sed(i)%t_cd = sed(i)%tau_cd
      sed(i)%t_ce = sed(i)%tau_ce
      sed(i)%rate = sed(i)%erate
    end do
  end if

  !set up the bedflag 
  ! bedflag = 0. => deactive erosion/deposition
  ! bedflag = 1. => normal erosion/deposition
  if(BEDFLAG_TYPE/=CNSTNT)then
    allocate(tmp(0:mt)) ; tmp = 0.0
    call LOAD_BEDFLAG(tmp) 
    bottom(0:mt,bflag) = tmp(0:mt) 
    deallocate(tmp) 
  else
    do i=1,m   
      bottom(i,bflag) = 1.0
    end do
  endif;


  !--------------------------------------------------
  !Initialize critical shear stress for erosion
  !--------------------------------------------------
 ! do i=1,nsed
 !   if(sed(i)%stype=='cohesive')then
 !     VAR => FIND_VAR(NC_START,'tce',FOUND)
 !     IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tce'&
 !          & IN THE STARTUP FILE OBJECT")
 !     CALL NC_CONNECT_AVAR(VAR, sed(i)%t_ce)
 !     CALL NC_READ_VAR(VAR)
 !     CALL NC_DISCONNECT(VAR)
 !   end if
 ! end do
  !--------------------------------------------------
  !Initialize critical shear stress for deposition
  !--------------------------------------------------
 ! do i=1,nsed
 !   if(sed(i)%stype=='cohesive')then
 !     VAR => FIND_VAR(NC_START,'tcd',FOUND)
 !     IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tcd'&
 !        & IN THE STARTUP FILE OBJECT")
 !     CALL NC_CONNECT_AVAR(VAR, sed(i)%t_cd)
 !     CALL NC_READ_VAR(VAR)
 !     CALL NC_DISCONNECT(VAR)
 !   end if
 ! end do
  !--------------------------------------------------
  !Initialize surface erosion rate
  !--------------------------------------------------
 ! do i=1,nsed
 !   VAR => FIND_VAR(NC_START,'ero',FOUND)
 !   if(sed(i)%stype=='cohesive')then
 !     IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'ero'&
 !        & IN THE STARTUP FILE OBJECT")
 !     CALL NC_CONNECT_AVAR(VAR, sed(i)%rate)
 !     CALL NC_READ_VAR(VAR)
 !     CALL NC_DISCONNECT(VAR)
 !   end if
 ! end do
  !--------------------------------------------------
  !Initialize Bed Properties
  !--------------------------------------------------

  Do k=1,Nbed
    Do i=1,m
       bed(i,k,iaged) = 0.00
       bed(i,k,ithck) = init_bed_thickness(k) 
       bed(i,k,iporo) = init_bed_porosity
     End Do
   End do

  !--------------------------------------------------
  !Initialize Bed_Frac properties
  !--------------------------------------------------

  !Do k=1,Nbed
  !  Do i=1,m
  !     Do ised=1,Nsed
  !       sed(ised)%frac(i,k)=1.0/float(nsed)   
  !     End Do
  !   End Do
  !End Do
  !--------------------------------------------------
  !Initialize Bed_Frac properties
  !Version 2.0 Skagit Sediment Setup
  !Mostly fine sand (dominant size class on flats)
  !--------------------------------------------------

  Do k=1,Nbed
    Do i=1,m
       Do ised=1,Nsed
         sed(ised)%frac(i,k)=init_bed_fraction(ised)  !1.0/float(nsed)   
       End Do
     End Do
  End Do

  !------------------------------------------------------
  ! store the fraction at bed surface
  !------------------------------------------------------
  !do ised=1,nst
  !   sed(ised)%frac(1:m,1:nbed)=sedbed%bed_frac(1:m,1:nbed,ised)
  !end do

  !--------------------------------------------------
  !Initialize Sediment Concentrations
  !--------------------------------------------------
  Do k=1,Kb
    Do i=1,m
      Do ised=1,Nsed
        sed(ised)%conc(i,k) = 0.0 
        sed(ised)%cnew(i,k) = sed(ised)%conc(i,k) 
      End Do
    End do
  End do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Init_Sed" 

  Return


# elif defined (CSTMS_SED)

  Use Mod_Prec 
  Use Mod_Sed_CSTMS
  Use Lims, only: m,kbm1,kb
  implicit none 
  integer  :: i,k,ised,kbed
  real(sp) :: cff1
!
!-----------------------------------------------------------------------
!  Initialize sediment bed properties.
!-----------------------------------------------------------------------
!
  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Init_Sed" 

  nnew=1
  nstp=2

  Do k=1,Nbed
     bed(1:m,k,iaged) = init_bed_age(k)
     bed(1:m,k,ithck) = init_bed_thickness(k)
     bed(1:m,k,iporo) = init_bed_porosity(k)
     if(COHESIVE_BED.or.SED_BIODIFF.or.MIXED_BED)then
       bed(1:m,k,ibtcr) = init_bed_tau_crit(k)
     endif
     if(SED_BIODIFF)then
       bed(1:m,k,idiff) = init_bed_biodiff(k)
     end if
     do i=1,nst
       sedbed%bed_frac(1:m,k,i)     = init_bed_fraction(i,k)

!      sedbed%bed_mass(1:m,k,i) = init_bed_mass(i,k)
!
!  Calculate mass so it is consistent with density, thickness, and
!  porosity.
!

       sedbed%bed_mass(1:m,k,1:2,i)=init_bed_thickness(k)*sed(i)%Srho&
            &*(1.0-init_bed_porosity(k))*init_bed_fraction(i,k)

     end do
  End do

  !--------------------------------------------------
  !Initialize Sediment Concentrations
  !--------------------------------------------------
  Do ised=1,Nsed
    sed(ised)%conc = 0.0_sp
    sed(ised)%cnew = sed(ised)%conc
  End do
  
!
!-----------------------------------------------------------------------
!  Compute initial total thickness for all sediment bed layers.
!-----------------------------------------------------------------------
!
  if(SED_MORPH)then

    do i=1,m
       sedbed%bed_thick0(i)=0.0_sp
       DO kbed=1,Nbed
         sedbed%bed_thick0(i)=sedbed%bed_thick0(i)+sedbed%bed(i,kbed,ithck)
       END DO
       sedbed%bed_thick(i,nnew)=sedbed%bed_thick0(i)
    end do

  end if

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Init_Sed" 

  Return

# endif


# endif
  End Subroutine Init_Sed
!=======================================================================
