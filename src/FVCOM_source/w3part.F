!/ ------------------------------------------------------------------- /
      MODULE W3PARTMD
#     if defined (WAVE_CURRENT_INTERACTION)
!/
!/                  +---------------------------------------+
!/                  | WAVE PARTITION FOR FVCOM    YANG ZHANG|
!/                  +---------------------------------------+
!  1. Purpose :
!
!     Interface to watershed partitioning routines.
!     Codes in this subroutine is made by Yang Zhang for wind definition and swells 
!
!  2. Variables and types :
!
!      Name                Type  Scope    Description
!     ----------------------------------------------------------------
!      MK, MTH             Int.  Private  Dimensions of stored neighour array.
!      NEIGH               I.A.  Private  Nearest Neighbor array.
!      hs_wind             R.A.  Out      Significant Wave Height of Wind Sea
!      dirdeg_wind         R.A.  Out      Wave Direction of Wind Sea
!      tpeak_wind          R.A.  Out      Peak Period  of Wind Sea
!      tpeak_wind_pos      R.A.  Out      IS(MSC) of Peak Period of Wind Sea
!      hs_swell_all        R.A.  Out      Significant Wave Height of Swell
!      dirdeg_swell_all    R.A.  Out      Wave Direction of Swell
!      tpeak_wind_all      R.A.  Out      Peak Period  of Swell
!      tpeak_wind_pos_all  R.A.  Out      IS(MSC) of Peak Period of Swell
!     ----------------------------------------------------------------
!      Note: IHMAX, HSPMIN, WSMULT, WSCUT and FLCOMB used from W3ODATMD.
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3PART    Subr. Public   Interface to watershed routines.
!      PTSORT    Subr. Public   Sort discretized image.
!      PTNGHB    Subr. Public   Defeine nearest neighbours.
!      PT_FLD    Subr. Public   Incremental flooding algorithm.
!      FIFO_ADD, FIFO_EMPTY, FIFO_FIRST
!                Subr. PT_FLD   Queue management.
!      PART_INFO Subr. Public   Compute wave parameters.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine traceing.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!
!zhangyang      USE W3ODATMD, ONLY: IHMAX, HSPMIN, WSMULT 
!
      PUBLIC
!
      INTEGER, PRIVATE              :: MK = -1, MTH = -1
      INTEGER, ALLOCATABLE, PRIVATE :: NEIGH(:,:)
      INTEGER, PRIVATE              :: IHMAX = 10000       !zhangyang 
      REAL,PRIVATE                  :: HSPMIN=max(0.0001,0.05), WSMULT = max(1.0,1.7)!zhangyang
      REAL                          :: WSCUT
      LOGICAL                       ::FLCOMB
!      REAL     ,ALLOCATABLE         :: hs_wind(:),dirdeg_wind(:),tpeak_wind(:)
!      INTEGER  ,ALLOCATABLE         :: tpeak_wind_pos(:)
!      REAL, ALLOCATABLE             :: hs_swell_all(:,:),dirdeg_swell_all(:,:),tpeak_swell_all(:,:)
!      INTEGER,ALLOCATABLE           :: tpeak_swell_pos_all(:,:)

!      REAL  ,PARAMETER              :: pi_w  = 4.*atan(1.)
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
!zhangyang      SUBROUTINE W3PART ( SPEC, UABS, UDIR, DEPTH, WN, NP, XP, DIMXP )
!      SUBROUTINE W3PART_NEW (SPEC,NK,NTH,NSPEC,DEPTH,DIMXP,WN,NIP)
      SUBROUTINE W3PART (SPEC,NK,NTH,NSPEC,DEPTH,DIMXP,WN,NIP)
!/
!/                  +-----------------------------------+
!/                  | Original from WAVEWATCH III       |   
!/                  |               USACE/NOAA          |
!/                  |          Barbara  Tracy           |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    28-Oct-2006 : Origination.                       ( version 3.10 )
!/    02-Dec-2010 : Adding a mapping PMAP between      ( version 3.14 )
!/                  original and combined partitions
!/                  ( M. Szyszka )
!/
!  1. Purpose :
!
!     Interface to watershed partitioning routines.
!
!  2. Method :
!
!     Watershed Algorithm of Vincent and Soille, 1991, implemented by
!     Barbara Tracy (USACE/ERDC) do NOAA/NCEP.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       SPEC    R.A.   I   2-D spectrum E(f,theta).
!       UABS    Real   I   Wind speed.
!       UDIR    Real   I   Wind direction.
!       DEPTH   Real   I   Water depth.
!       WN      R.A.   I   Wavenumebers do each frequency.
!       NP      Int.   O   Number of partitions.
!                           -1 : Spectrum without minumum energy.
!                            0 : Spectrum with minumum energy.
!                                but no partitions.
!       XP      R.A.   O   Parameters describing partitions.
!                          Entry '0' contains entire spectrum.
!       DIMXP   Int.   I   Second dimension of XP.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - To achieve minimum storage but guaranteed storage of all 
!       partitions DIMXP = ((NK+1)/2) * ((NTH-1)/2)
!
!  8. Structure :
!
!  9. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE SWCOMM2               !yzhang_w3
      USE SWCOMM3
      USE VARS_WAVE
      USE M_GENARR, ONLY : SPCSIG
#     if defined (MULTIPROCESSOR)
      USE MOD_PAR
#     endif
#     if defined (EXPLICIT)
      USE MOD_ACTION_EX
#     else             
      USE MOD_ACTION_IM
#     endif
!      USE MOD_SPARSE_TIMESERIES
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER  ::NSPEC
      INTEGER          :: NP
      INTEGER                        :: DIMXP,INDX!zhangyang
      REAL, INTENT(IN)              :: SPEC(NK,NTH), WN(NK)    
                                      
      REAL             :: XP(6,0:DIMXP)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      REAL     ,PARAMETER     :: AA=4e-5 ,BB=2e-3,gf=1.5
      INTEGER                 :: ITH, IMI(NSPEC), IMD(NSPEC),         &
                                 IMO(NSPEC), IND(NSPEC), NP_MAX,      &
                                 NIP, IT(1), INDEX(DIMXP), NWS,        &
                                 IPW, IPTT, ISP
      INTEGER                 :: PMAP(DIMXP),NK,NTH,IMO_2D(MDC,MSC),IMO_2D_SW(MDC,MSC)
      INTEGER                 :: J,K,I,MAX_IMO
!      INTEGER, SAVE           :: IENT = 0
      REAL                    :: ZP(NSPEC), ZMIN, ZMAX, ZZZ(NSPEC),     &
                                 FACT, WSMAX, HSMAX
      REAL                    :: TP(6,DIMXP)
      REAL                    :: AWX_TMP,AWY_TMP,UDIR,UABS,DEPTH
      Real(SP) ,allocatable   :: parts(:,:,:),hs_parts(:),dirdeg_parts(:),tpeak_parts(:)
      Real     ,allocatable   :: order(:),energy_max(:),parts_swell(:,:,:),fpeak_parts(:)
      Real     ,allocatable   :: dtheta(:),fpeak_threshold(:)
!      Real     ,allocatable   :: hs_swell_all(:),tpeak_swell_all(:),dirdeg_swell_all(:)
!      integer  ,allocatable   :: tpeak_swell_pos_all(:)
      Real(SP)                :: hs_swell,tpeak_swell,dirdeg_swell,part(MDC,MSC)
      integer                 :: tpeak_swell_pos
      REAL(SP)                :: part_wind(MDC,MSC)!,hs_wind,dirdeg_wind,tpeak_wind,tpeak_wind_pos
      REAL                    :: etot_swell,fy_mean,fy2_mean,fx_mean,fx2_mean,gama_f2,etot_tmp
!      REAL                    :: tpeak_energy,tmp_r,part_tmp(MDC,MSC),part_swell(MDC,MSC)
      REAL                    :: tpeak_energy,tmp_r,part_tmp(MDC,MSC)
      REAL(SP)                :: part_swell(MDC,MSC)
      integer                 :: tpeak_pos,tmp_i,count_swells
      integer                 :: id,is,ii,oo,jj,kk,try_time,try_t
      real                    :: energy_min,energe_edge,part1(MDC,MSC),part2(MDC,MSC),theta1,theta2,fp1,fp2
      integer                 :: adjacent,mergee,adjacent_pos(8,2),theta1_pos,theta2_pos,fp1_pos,fp2_pos
      real                    :: fpx1,fpx2,fpy1,fpy2,d_f2,etot_the
      real                    :: parts_swell_tmp1(MDC,MSC),parts_swell_tmp2(MDC,MSC)
!      real                    :: peak_energy1,peak_energy2,hs_parto,dirdeg_parto,tpeak_parto
      real                    :: peak_energy1,peak_energy2
      real(sp)                :: hs_parto,dirdeg_parto,tpeak_parto
!/
!/ ------------------------------------------------------------------- /
! 0.  Initializations
!
!/S      CALL STRACE (IENT, 'W3PART')
!
      FLCOMB = .TRUE.!zhangyang
      WSCUT  = MIN(1.0001,MAX(0.,0.333))!zhangyang
      NP     = 0
      XP     = 0.
      IF(VARWI)THEN                              !yzhang_w3
      AWX_TMP=COMPDA(NIP,16)!JWX2 zhangyang      !yzhang_w3
      AWY_TMP=COMPDA(NIP,17)!JWY2 zhangyang      !yzhang_w3
      ELSE                                       !yzhang_w3
        AWX_TMP = U10 * COS(WDIC)                !yzhang_w3
        AWY_TMP = U10 * SIN(WDIC)                !yzhang_w3
      END IF                                     !yzhang_w3

      UABS=SQRT(AWX_TMP**2+AWY_TMP**2)!zhangyang
      UDIR=ATAN2(AWY_TMP,AWX_TMP)!zhangyang
       if(UABS.lt.1e-8)then
          UABS=1
          UDIR=0
       end if 
!      DIMXP = ((NK+1)/2) * ((NTH-1)/2)
       DIMXP = (NK/2) * (NTH/2)
!DIMXP = NK * NTH
!
! -------------------------------------------------------------------- /
! 1.  Process input spectrum
! 1.a 2-D to 1-D spectrum
!
      DO ITH=1, NTH
        ZP(1+(ITH-1)*NK:ITH*NK) = SPEC(:,ITH)
        END DO
!
! 1.b Invert spectrum and 'digitize'
!
      ZMIN   = MINVAL ( ZP ) 
      ZMAX   = MAXVAL ( ZP ) 
      IF ( ZMAX-ZMIN .LT. 1.E-9 ) RETURN
!
      ZZZ      = ZMAX - ZP
!
      FACT   = REAL(IHMAX-1) / ( ZMAX - ZMIN )
      IMI    = MAX ( 1 , MIN ( IHMAX , NINT ( 1. + ZZZ*FACT ) ) )
! 1.c Sort digitized image
!
      CALL PTSORT ( IMI, IND, IHMAX,NSPEC )
! -------------------------------------------------------------------- /
! 2.  Perdom partitioning
! 2.a Update nearest neighbor info as needed.
!
      CALL PTNGHB(NK, NTH, NSPEC)
!
! 2.b Incremental flooding
!
      CALL PT_FLD ( IMI, IND, IMO, ZP, NP_MAX ,NSPEC)
!
! 2.c Compute parameters per partition
!     
! 
      MAX_IMO=maxval(IMO)    
      allocate (hs_parts(MAX_IMO));hs_parts=0.0
      allocate (dirdeg_parts(MAX_IMO));dirdeg_parts=0.0
      allocate (tpeak_parts(MAX_IMO));tpeak_parts=0.0
      allocate (parts(MAX_IMO,MDC,MSC));parts=0.0
      allocate (order(MAX_IMO));order=0
      allocate (energy_max(MAX_IMO));energy_max=0.0
      allocate (parts_swell(MAX_IMO,MDC,MSC));parts_swell=0.0
      allocate (fpeak_parts(MAX_IMO));fpeak_parts=0.0
      allocate (dtheta(MAX_IMO));dtheta=0.0
      allocate (fpeak_threshold(MAX_IMO));fpeak_threshold=0.0
         DO J=1,MDC
            DO K=1,MSC
               IMO_2D(J,K)=IMO(K+MSC*(J-1))
            END DO
         END DO
      DO I=1,MAX_IMO       
         DO J=1,MDC
            DO K=1,MSC
               IF (IMO_2D(J,K).eq.I)THEN
                  parts(I,J,K)=SPEC(K,J)
                  if (parts(I,J,K).lt.0.0)parts(I,J,K)=0
               END IF
            END DO
         END DO
      END DO
      DO I=1,MAX_IMO
         DO J=1,MDC
            DO K=1,MSC
               part(J,K)=parts(I,J,K)
            END DO
         END DO
            if(sum(part).gt.1e-4.and.sum(part).lt.1e+4)then
            hs_parto=0.0
            call part_info(part,SPCSIG,&
                 hs_parto,dirdeg_parto,tpeak_parto)   
            end if
            if(hs_parto.gt.0)then
               hs_parts(I)=hs_parto
               dirdeg_parts(I)=dirdeg_parto
               tpeak_parts(I)=tpeak_parto
               dirdeg_parts(I)=dirdeg_parts(I)/180*pi_w
            end if
      END DO
      DO I=1,MAX_IMO
         fpeak_parts(I)=1.0/tpeak_parts(I)
      END DO         
      part_wind=0.0
      DO I=1,MAX_IMO
         dtheta(I)=abs(UDIR-dirdeg_parts(I))
         if (dtheta(I).gt.pi_w)dtheta(I)=pi_w*2-dtheta(I)
         if (dtheta(I).le.pi_w/2)then
            fpeak_threshold(I)=9.8/2/pi_w/(gf*UABS*cos(dtheta(I)))
         if (fpeak_parts(I).ge.fpeak_threshold(I))then
            do j=1,MDC
               do k=1,MSC
                  part_wind(j,k)=part_wind(j,k)+parts(I,j,k)
                  parts(I,j,k)=0.0
               end do
            end do
            hs_parts(I)=0.0
            dirdeg_parts(I)=0.0
            tpeak_parts(I)=0.0
        do J=1,MDC
           do K=1,MSC
              if (IMO_2D(J,K).eq.I)then
                 IMO_2D(J,K)=0
              end if 
           end do
        end do
        end if 
        end if 
      END DO
! calculate_wind_wave_indomation
      if (sum(part_wind).gt.1e-4.and.sum(part).lt.1e+4)then
         call  part_info(part_wind,SPCSIG,hs_wind(NIP),dirdeg_wind(NIP),tpeak_wind(NIP))
         tpeak_energy=part_wind(1,1)
         etot_tmp=0.0
         do j=1,MDC
            do k=1,MSC
               etot_tmp=etot_tmp+part_wind(j,k)
               if (part_wind(j,k).gt.tpeak_energy)then
                  tpeak_energy=part_wind(j,k)
                  tpeak_wind_pos(NIP)=k
               end if 
            end do
         end do
      if (etot_tmp .lt. (AA/((SPCSIG(tpeak_wind_pos(NIP)))**4+BB)) )then
         hs_wind(NIP)=-999
         dirdeg_wind(NIP)=-999
         tpeak_wind(NIP)=-999
         tpeak_wind_pos(NIP)=-999
      elseif (hs_wind(NIP).eq.0.0) then
         hs_wind(NIP)=-999
         dirdeg_wind(NIP)=-999
         tpeak_wind(NIP)=-999
         tpeak_wind_pos(NIP)=-999
      end if 
      else
         hs_wind(NIP)=-999
         dirdeg_wind(NIP)=-999
         tpeak_wind(NIP)=-999
         tpeak_wind_pos(NIP)=-999
      end if 
      DO I=1,MAX_IMO
         part_tmp(:,:)=parts(I,:,:)
         energy_max(I)=sum(part_tmp)
         order(I)=I
      END DO
      do i=1,size(energy_max)-1
        do j=i+1,size(energy_max)
            if(energy_max(j)>energy_max(i))then
                tmp_r=energy_max(i)
                energy_max(i)=energy_max(j)
                energy_max(j)=tmp_r
                tmp_i=order(i)
                order(i)=order(j)
                order(j)=tmp_i
            end if
        end do
      end do
      DO I=1,MAX_IMO
         parts_swell(I,:,:)=parts(order(I),:,:) 
      END DO
      IMO_2D_SW=0
      DO I=1,MAX_IMO
         DO J=1,MDC
            DO K=1,MSC
               if (IMO_2D(J,K).eq.order(I)) IMO_2D_SW(J,K)=I
            END DO
         END DO
      END DO
      part_swell=0.0 
      DO I=1,MAX_IMO
         part_swell(:,:)=part_swell(:,:) + parts_swell(I,:,:)
      END DO 
      if (sum(part_swell).gt.1e-4.and.sum(part_swell).lt.1e+4)then
         call part_info(part_swell,SPCSIG,hs_swell,dirdeg_swell,tpeak_swell) 
         etot_swell=0.0
         tpeak_energy=part_swell(1,1)
         do j=1,MDC
            do k=1,MSC
               etot_swell=etot_swell+part_swell(j,k)
               if (part_swell(j,k).gt.tpeak_energy)then
                  tpeak_energy=part_swell(j,k)
                  tpeak_swell_pos=k
               end if
            end do
         end do  
      end if 
!         %%%%%%%%%%%%fy_mean%%%%%%%%%%%%%
         fy_mean=0;
         do id=1,MDC
            do is=1,MSC
               fy_mean=fy_mean+SPCSIG(is)*part_swell(id,is)*sin((id-0.5)/MDC*2*pi_w);
            end do
         end do

!         %%%%%%%%%%%%fy2_mean%%%%%%%%%%%%%
         fy2_mean=0;
         do id=1,MDC
            do is=1,MSC
               fy2_mean=fy2_mean+SPCSIG(is)**2*part_swell(id,is)*(sin((id-0.5)/MDC*2*pi_w))**2;
            end do
         end do

!         %%%%%%%%%%%%fx_mean%%%%%%%%%%%%%
         fx_mean=0;
         do id=1,MDC
            do is=1,MSC
               fx_mean=fx_mean+SPCSIG(is)*part_swell(id,is)*cos((id-0.5)/MDC*2*pi_w);
            end do
         end do

!         %%%%%%%%%%%%fx2_mean%%%%%%%%%%%%%
         fx2_mean=0;
         do id=1,MDC
            do is=1,MSC
               fx2_mean=fx2_mean+SPCSIG(is)**2*part_swell(id,is)*(cos((id-0.5)/MDC*2*pi_w))**2;
            end do
         end do
         gama_f2=fx2_mean+fy2_mean-fx_mean**2-fy_mean**2
!%%%%%%%%%%start_swell_merge%%%%%%%%%%%%%%%
         if (maxval(part_swell).gt.0)then
            try_time=MAX_IMO
            DO try_t=1,try_time
               DO I=1,try_t-1
                  DO J=(i+1),try_t
                  parts_swell_tmp1(:,:)=parts_swell(i,:,:)
                  parts_swell_tmp2(:,:)=parts_swell(j,:,:)
                  if((sum(parts_swell_tmp1).gt.1e-6).and.sum(parts_swell_tmp2).gt.1e-6)then
                  energy_min=min(maxval(parts_swell_tmp1),maxval(parts_swell_tmp2))
                  adjacent=0
                  mergee=0
                  energe_edge=10000
                     do ii=1,MDC
                        do oo=1,MSC
                           if (IMO_2D_SW(ii,oo).eq.I)then
                              adjacent_pos(1,:)=(/ii , oo+1/)
                              adjacent_pos(2,:)=(/ii-1 , oo+1/)
                              adjacent_pos(3,:)=(/ii-1 , oo/)
                              adjacent_pos(4,:)=(/ii-1 , oo-1/)
                              adjacent_pos(5,:)=(/ii , oo-1/)
                              adjacent_pos(6,:)=(/ii+1 , oo-1/)
                              adjacent_pos(7,:)=(/ii+1 , oo/)
                              adjacent_pos(8,:)=(/ii+1 , oo+1/)
                              do k=1,8
                                 if (adjacent_pos(k,1).eq.0)then
                                    adjacent_pos(k,1)=MDC
                                 elseif (adjacent_pos(k,1).eq.MDC+1)then
                                    adjacent_pos(k,1)=1
                                 end if
                                 if (adjacent_pos(k,2).eq.0)then
                                    adjacent_pos(k,2)=9999
                                 elseif (adjacent_pos(k,2).gt.MSC)then
                                    adjacent_pos(k,2)=9999
                                 end if
                              end do
                              do k=1,8
                                 if (adjacent_pos(k,2).lt.9999)then
                                    if (IMO_2D_SW(adjacent_pos(k,1),adjacent_pos(k,2)).eq.j)then
                                       adjacent=1
                                       energe_edge=min(energe_edge,parts_swell(i,ii,oo))
                                       energe_edge=min(energe_edge,parts_swell(j,adjacent_pos(k,1),adjacent_pos(k,2)))
                                    end if 
                                 end if
                              end do
                           end if 
                        end do
                     end do
                     if (adjacent.eq.1)then
                        part1(:,:)=parts(i,:,:)
                        part2(:,:)=parts(j,:,:)
                        theta1_pos=0
                        fp1_pos=0
                        theta2_pos=0
                        fp2_pos=0
                        peak_energy1=0
                        peak_energy2=0
                        do jj=1,MDC
                           do kk=1,MSC
                              if (part1(jj,kk).gt.peak_energy1)then
                                 peak_energy1=part1(jj,kk)
                                 fp1_pos=kk
                                 theta1_pos=jj
                              end if
                              if (part2(jj,kk).gt.peak_energy2)then
                                 peak_energy2=part2(jj,kk)
                                 fp2_pos=kk
                                 theta2_pos=jj
                               end if
                           end do
                        end do
                        theta1=(theta1_pos-0.5)/MDC*2*pi
                        theta2=(theta2_pos-0.5)/MDC*2*pi
                        fp1=spcsig(fp1_pos)
                        fp2=spcsig(fp2_pos)
                        fpx1=fp1*cos(theta1)
                        fpy1=fp1*sin(theta1)
                        fpx2=fp2*cos(theta2)
                        fpy2=fp2*sin(theta2)
                        d_f2=(fpx1-fpx2)**2+(fpy1-fpy2)**2
                        if (energe_edge.gt.0.75*energy_min.and.d_f2.le.0.5*gama_f2)then
                           mergee=1;
                        end if
                        if (mergee.eq.1)then
                           parts_swell(i,:,:)=parts_swell(i,:,:)+parts_swell(j,:,:)
                           if (j.lt.MAX_IMO)then
                           DO k=j,MAX_IMO-1
                           parts_swell(k,:,:)=parts_swell(k+1,:,:)
                           END DO
                           parts_swell(MAX_IMO,:,:)=0
                           else
                           parts_swell(MAX_IMO,:,:)=0
                           end if
                           exit
                        end if              
                     end if                       
                  end if
                  end do
                  if (mergee.eq.1)then
                  exit
                  end if   
               end do

            end do
            count_swells=MAX_IMO;
         end if
!%%%%%%%%%%%%%%end_swell_merge%%%%%%%%%%%%%%%%%  
       DO j=1,MAX_IMO
          part_swell(:,:)=parts_swell(j,:,:) 
           if (sum(part_swell).gt.1e-4.and.sum(part).lt.1e+4)then
            call  part_info(part_swell,SPCSIG,hs_swell,dirdeg_swell,tpeak_swell)
            tpeak_energy=part_swell(1,1)
            etot_tmp=0.0
            do i=1,MDC
               do k=1,MSC
                  etot_tmp=etot_tmp+part_swell(i,k)
                  if (part_swell(i,k).gt.tpeak_energy)then
                     tpeak_energy=part_swell(i,k)
                     tpeak_swell_pos=k
                  end if 
               end do
            end do
            if (etot_tmp.lt.(AA/((SPCSIG(tpeak_swell_pos))**4+BB)))then
               hs_swell=0!-999
               dirdeg_swell=0!-999
               tpeak_swell=0!-999
               tpeak_swell_pos=0!-999
            elseif (hs_swell.eq.0.0) then
               hs_swell=0!-999
               dirdeg_swell=0!-999
               tpeak_swell=0!-999
               tpeak_swell_pos=0!-999
            end if 
         else
            hs_swell=0!-999
            dirdeg_swell=0!-999
            tpeak_swell=0!-999
            tpeak_swell_pos=0!-999
         end if 
         hs_swell_all(NIP,j)=hs_swell
         dirdeg_swell_all(NIP,j)=dirdeg_swell
         tpeak_swell_all(NIP,j)=tpeak_swell  
         tpeak_swell_pos_all(NIP,j)=tpeak_swell_pos
       END DO 
      deallocate (hs_parts,dirdeg_parts,tpeak_parts,parts,order,energy_max,parts_swell,fpeak_parts)
      deallocate (dtheta,fpeak_threshold)
!
      RETURN
!/
!/ End of W3PART ----------------------------------------------------- /
!/
      END SUBROUTINE W3PART
!      END SUBROUTINE W3PART_NEW
!/ ------------------------------------------------------------------- /
      SUBROUTINE PTSORT ( IMI, IND, IHMAX ,NSPEC)
!/
!/                  +-----------------------------------+
!/                  | Original from WAVEWATCH III       |   
!/                  |               USACE/NOAA          |
!/                  |          Barbara  Tracy           |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         19-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    19-Oct-2006 : Origination.                       ( version 3.10 )
!/
!  1. Purpose :
!
!     This subroutine sorts the image data in ascending order.
!     This sort original to F.T.Tracy (2006)
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMI     I.A.   I   Input discretized spectrum.
!       IND     I.A.   O   Sorted data.
!       IHMAX   Int.   I   Number of integer levels.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
!/S      USE W3SERVMD, ONLY: STRACE
!
!      USE W3GDATMD, ONLY: NSPEC
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER                 ::NSPEC
      INTEGER, INTENT(IN)      :: IHMAX, IMI(NSPEC)
      INTEGER, INTENT(OUT)     :: IND(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, IN, IV
!/S      INTEGER, SAVE           :: IENT = 0
      INTEGER                 :: NUMV(IHMAX), IADDR(IHMAX),           &
                                 IORDER(NSPEC)
!/
!/S      CALL STRACE (IENT, 'PTSORT')
!
! -------------------------------------------------------------------- /
! 1.  Occurences per height
!
      NUMV   = 0
      DO I=1, NSPEC
        NUMV(IMI(I)) = NUMV(IMI(I)) + 1
        END DO
!
! -------------------------------------------------------------------- /
! 2.  Starting address per height
!
      IADDR(1) = 1
      DO I=1, IHMAX-1
        IADDR(I+1) = IADDR(I) + NUMV(I)
      END DO
!
! -------------------------------------------------------------------- /
! 3.  Order points
!
      DO I=1, NSPEC
        IV        = IMI(I)
        IN        = IADDR(IV)
        IORDER(I) = IN
        IADDR(IV) = IN + 1
        END DO
!
! -------------------------------------------------------------------- /
! 4.  Sort points
!
      DO I=1, NSPEC
        IND(IORDER(I)) = I
        END DO
!
      RETURN
!/
!/ End of PTSORT ----------------------------------------------------- /
!/
      END SUBROUTINE PTSORT
!/ ------------------------------------------------------------------- /
      SUBROUTINE PTNGHB(NK, NTH, NSPEC) 
!/
!/                  +-----------------------------------+
!/                  | Original from WAVEWATCH III       |   
!/                  |               USACE/NOAA          |
!/                  |          Barbara  Tracy           |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         20-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    20-Oct-2006 : Origination.                       ( version 3.10 )
!/
!  1. Purpose :
!
!     This subroutine computes the nearest neighbors do each grid
!     point. Wrapping of directional distribution (0 to 360)is taken
!     care of using the nearest neighbor system
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMI     I.A.   I   Input discretized spectrum.
!       IMD     I.A.   O   Sorted data.
!       IHMAX   Int.   I   Number of integer levels.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
!/S      USE W3SERVMD, ONLY: STRACE
!
!      USE W3GDATMD, ONLY: NK, NTH, NSPEC
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!     INTEGER, INTENT(IN)      :: IHMAX, IMI(NSPEC)
!     INTEGER, INTENT(IN)      :: IMD(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: N, J, I, K,NK, NTH, NSPEC
!/S      INTEGER, SAVE           :: IENT = 0
!/
!/S      CALL STRACE (IENT, 'PTNGHB')
!
! -------------------------------------------------------------------- /
! 1.  Check on need of processing
!
      IF ( MK.EQ.NK .AND. MTH.EQ.NTH ) RETURN
!
      IF ( MK.GT.0 ) DEALLOCATE ( NEIGH )
      ALLOCATE ( NEIGH(9,NSPEC) )
      MK     = NK
      MTH    = NTH
!
! -------------------------------------------------------------------- /
! 2.  Build map
!
      NEIGH  = 0
!
! ... Base loop
!
      DO N = 1, NSPEC
!
        J      = (N-1) / NK + 1
        I      = N - (J-1) * NK
        K      = 0
!
! ... Point at the left(1)
!
        IF ( I .NE. 1 ) THEN
            K           = K + 1
            NEIGH(K, N) = N - 1
          END IF
!
! ... Point at the right (2)
!
        IF ( I .NE. NK ) THEN 
            K           = K + 1
            NEIGH(K, N) = N + 1
          END IF
!
! ... Point at the bottom(3)
!
        IF ( J .NE. 1 ) THEN
            K           = K + 1
            NEIGH(K, N) = N - NK
          END IF
!
! ... ADD Point at bottom_wrap to top
!
        IF ( J .EQ. 1 ) THEN
            K          = K + 1
            NEIGH(K,N) = NSPEC - (NK-I)
          END IF
!
! ... Point at the top(4)
!
        IF ( J .NE. NTH ) THEN
            K           = K + 1
            NEIGH(K, N) = N + NK
          END IF
!
! ... ADD Point to top_wrap to bottom
!
         IF ( J .EQ. NTH ) THEN
             K          = K + 1
             NEIGH(K,N) = N - (NTH-1) * NK
            END IF
!
! ... Point at the bottom, left(5)
!
        IF ( (I.NE.1) .AND. (J.NE.1) ) THEN
            K           = K + 1
            NEIGH(K, N) = N - NK - 1
          END IF
!
! ... Point at the bottom, left with wrap.
!
         IF ( (I.NE.1) .AND. (J.EQ.1) ) THEN
             K          = K + 1
             NEIGH(K,N) = N - 1 + NK * (NTH-1)
           END IF
!
! ... Point at the bottom, right(6)
!
        IF ( (I.NE.NK) .AND. (J.NE.1) ) THEN
            K           = K + 1
            NEIGH(K, N) = N - NK + 1
          END IF
!
! ... Point at the bottom, right with wrap
!
        IF ( (I.NE.NK) .AND. (J.EQ.1) ) THEN
            K           = K + 1
            NEIGH(K,N) = N + 1 + NK * (NTH - 1)
          END  IF
!
! ... Point at the top, left(7)
!
        IF ( (I.NE.1) .AND. (J.NE.NTH) ) THEN
            K           = K + 1
            NEIGH(K, N) = N + NK - 1
          END IF
!
! ... Point at the top, left with wrap
!
         IF ( (I.NE.1) .AND. (J.EQ.NTH) ) THEN
             K           = K + 1
             NEIGH(K,N) = N - 1 - (NK) * (NTH-1)
           END IF
!
! ... Point at the top, right(8)
!
        IF ( (I.NE.NK) .AND. (J.NE.NTH) ) THEN
            K           = K + 1
            NEIGH(K, N) = N + NK + 1
          END IF
!
! ... Point at top, right with wrap
!
!
        IF ( (I.NE.NK) .AND. (J.EQ.NTH) ) THEN
            K           = K + 1
            NEIGH(K,N) = N + 1 - (NK) * (NTH-1)
          END IF
!
        NEIGH(9,N) = K
!
        END DO
!
      RETURN
!/
!/ End of PTNGHB ----------------------------------------------------- /
!/
      END SUBROUTINE PTNGHB
!/ ------------------------------------------------------------------- /
      SUBROUTINE PT_FLD ( IMI, IND, IMO, ZP, NPART ,NSPEC)
!/
!/                  +-----------------------------------+
!/                  | Original from WAVEWATCH III       |   
!/                  |               USACE/NOAA          |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Nov-2006 !
!/                  +-----------------------------------+
!/
!/    01-Nov-2006 : Origination.                       ( version 3.10 )
!/
!  1. Purpose :
!
!     This subroutine does incremental flooding of the image to
!     determine the watershed image.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMI     I.A.   I   Input discretized spectrum.
!       IND     I.A.   I   Sorted addresses.
!       IMO     I.A.   O   Output partitioned spectrum.
!       ZP      R.A.   I   Spectral array.
!       NPART   Int.   O   Number of partitions found.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
!
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMI(NSPEC), IND(NSPEC)
      INTEGER, INTENT(OUT)    :: IMO(NSPEC), NPART
      REAL, INTENT(IN)        :: ZP(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: MASK, INIT, IWSHED, IMD(NSPEC),      &
                                 IC_LABEL, IFICT_PIXEL, M, IH, MSAVE, &
                                 IP, I, IPP, IC_DIST, IEMPTY, IPPP,   &
                                 JL, JN, IPTT, J
      INTEGER                 :: IQ(NSPEC), IQ_START, IQ_END,NSPEC
      REAL                    :: ZPMAX, EP1, DIFF
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      MASK        = -2
      INIT        = -1
      IWSHED      =  0
      IMO         = INIT
      IC_LABEL    =  0
      IMD         =  0
      IFICT_PIXEL = -100
!
      IQ_START    =  1
      IQ_END      =  1
!
      ZPMAX       = MAXVAL ( ZP )
!
! -------------------------------------------------------------------- /
! 1.  Loop over levels
!
      M      =  1
!
      DO IH=1, IHMAX
        MSAVE  = M
!
! 1.a Pixels at level IH
!
        DO
          IP     = IND(M)
          IF ( IMI(IP) .NE. IH ) EXIT
!
!     Flag the point, if it stays flagge, it is a separate minimum.
!
          IMO(IP) = MASK
!
!     Consider neighbors. If there is neighbor, set distance and add
!     to queue.
!
          DO I=1, NEIGH(9,IP)
            IPP    = NEIGH(I,IP)
            IF ( (IMO(IPP).GT.0) .OR. (IMO(IPP).EQ.IWSHED) ) THEN
                IMD(IP) = 1
                CALL FIFO_ADD (IP)
                EXIT
              END IF
            END DO
!
          IF ( M+1 .GT. NSPEC ) THEN
              EXIT
            ELSE
              M = M + 1
            END IF
!
          END DO
!
! 1.b Process the queue
!
        IC_DIST = 1
        CALL FIFO_ADD (IFICT_PIXEL)
!
        DO
          CALL FIFO_FIRST (IP)
!
!     Check do end of processing
!
          IF ( IP .EQ. IFICT_PIXEL ) THEN
              CALL FIFO_EMPTY (IEMPTY)
              IF ( IEMPTY .EQ. 1 ) THEN
                  EXIT
                ELSE
                  CALL FIFO_ADD (IFICT_PIXEL)
                  IC_DIST = IC_DIST + 1
                  CALL FIFO_FIRST (IP)
                END IF
            END IF
!
!     Process queue
!
          DO I=1, NEIGH(9,IP)
            IPP = NEIGH(I,IP)
!
!     Check do labeled watersheds or basins
!
            IF ( (IMD(IPP).LT.IC_DIST) .AND. ( (IMO(IPP).GT.0) .OR.  &
                 (IMO(IPP).EQ.IWSHED))) THEN
!
                IF ( IMO(IPP) .GT. 0 ) THEN
!
                    IF ((IMO(IP) .EQ. MASK) .OR. (IMO(IP) .EQ. &
                        IWSHED)) THEN
                        IMO(IP) = IMO(IPP)
                      ELSE IF (IMO(IP) .NE. IMO(IPP)) THEN
                        IMO(IP) = IWSHED
                      END IF
!
                  ELSE IF (IMO(IP) .EQ. MASK) THEN
!
                    IMO(IP) = IWSHED
!
                  END IF
!
              ELSE IF ( (IMO(IPP).EQ.MASK) .AND. (IMD(IPP).EQ.0) ) THEN
!
                 IMD(IPP) = IC_DIST + 1
                 CALL FIFO_ADD (IPP)
!
              END IF
!
            END DO
!
          END DO
!
! 1.c Check do mask values in IMO to identify new basins
!
        M = MSAVE
!
        DO
          IP     = IND(M)
          IF ( IMI(IP) .NE. IH ) EXIT
          IMD(IP) = 0
!
          IF (IMO(IP) .EQ. MASK) THEN
!
! ... New label do pixel
!
              IC_LABEL = IC_LABEL + 1
              CALL FIFO_ADD (IP)
              IMO(IP) = IC_LABEL
!
! ... and all connected to it ...
!
              DO
                CALL FIFO_EMPTY (IEMPTY)
                IF ( IEMPTY .EQ. 1 ) EXIT
                CALL FIFO_FIRST (IPP)
!
                DO I=1, NEIGH(9,IPP)
                  IPPP   = NEIGH(I,IPP)
                  IF ( IMO(IPPP) .EQ. MASK ) THEN
                      CALL FIFO_ADD (IPPP)
                      IMO(IPPP) = IC_LABEL
                    END IF
                  END DO
!
                END DO
!
            END IF
!
          IF ( M + 1 .GT. NSPEC ) THEN
              EXIT
            ELSE
              M = M + 1
            END IF
!
          END DO
!
        END DO
!
! -------------------------------------------------------------------- /
! 2.  Find nearest neighbor of 0 watershed points and replace
!     use original input to check which group to affiliate with 0
!     Soring changes first in IMD to assure symetry in adjustment.
!
      DO J=1, 5
        IMD    = IMO
        DO JL=1 , NSPEC
          IPTT    = -1
          IF ( IMO(JL) .EQ. 0 ) THEN
              EP1    = ZPMAX
              DO JN=1, NEIGH (9,JL)
                DIFF   = ABS ( ZP(JL) - ZP(NEIGH(JN,JL)))
                IF ( (DIFF.LE.EP1) .AND. (IMO(NEIGH(JN,JL)).NE.0) ) THEN
                    EP1    = DIFF
                    IPTT    = JN
                  END IF
                END DO
              IF ( IPTT .GT. 0 ) IMD(JL) = IMO(NEIGH(IPTT,JL))
            END IF
          END DO
        IMO    = IMD
        IF ( MINVAL(IMO) .GT. 0 ) EXIT
        END DO
!
      NPART = IC_LABEL
!
      RETURN
!
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIFO_ADD ( IV )
!
!     Add point to FIFO queue.
!
      INTEGER, INTENT(IN)      :: IV
!
      IQ(IQ_END) = IV
!
      IQ_END = IQ_END + 1
      IF ( IQ_END .GT. NSPEC ) IQ_END = 1
!
      RETURN
      END SUBROUTINE FIFO_ADD
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIFO_EMPTY ( IEMPTY )
!
!     Check if queue is empty.
!
      INTEGER, INTENT(OUT)     :: IEMPTY
!
      IF ( IQ_START .NE. IQ_END ) THEN
        IEMPTY = 0
      ELSE
        IEMPTY = 1
      END IF
!
      RETURN
      END SUBROUTINE FIFO_EMPTY
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIFO_FIRST ( IV )
!
!     Get point out of queue.
!
      INTEGER, INTENT(OUT)     :: IV
!
      IV = IQ(IQ_START)
!
      IQ_START = IQ_START + 1
      IF ( IQ_START .GT. NSPEC ) IQ_START = 1
!
      RETURN
      END SUBROUTINE FIFO_FIRST
!/
!/ End of PT_FLD ----------------------------------------------------- /
!/
      END SUBROUTINE PT_FLD
!/ ------------------------------------------------------------------- /
!/ ------------------------------------------------------------------- /
!      SUBROUTINE part_info(part,ddir,frintf,SPCSIG1,MDC,MSC,hs_part,dirdeg_part,tpeak_part) 
      SUBROUTINE part_info(part,SPCSIG1,hs_part,dirdeg_part,tpeak_part) 
!/
!/                  +-----------------------------------+
!/                  | wave parameters        yang zhang !
!/                  +-----------------------------------+
!/
!/    28-Oct-2006 : Origination.                       ( version 3.10 )
!/    02-Dec-2010 : Adding a mapping PMAP between      ( version 3.14 )
!/                  original and combined partitions
!/                  ( M. Szyszka )
!/
!  1. Purpose :
!
!     Calculater the wave parameters
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       part          R.A.   I   2-D spectrum E(f,theta).
!       hs_part       R      O   Significant Wave height of this part
!       dirdeg_part   R      O   wave direction of this part
!       tpeak_part    R      O   peak period of this part
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     ----------------------------------------------------------------
      USE MOD_PREC
      USE SWCOMM3
      IMPLICIT NONE
      INTEGER      :: is,id
      INTEGER      ::isigm
!      REAL(SP)     ::part(MDC,MSC),ddir,frintf,SPCSIG1(MSC),hs_part,dirdeg_part,tpeak_part
      REAL(SP)     ::part(MDC,MSC),hs_part,dirdeg_part,tpeak_part
      REAL         :: SPCSIG1(MSC)
      REAL         ::emax,etd,etot,spcdir(MDC,3),eex,eey,ead,ds,edi,eftail,ehfr
!      %%%%%%%%%%%tpeak_part%%%%%%%%%%%%%%%%%
       emax=0;
       isigm=-1;
       do is=1,MSC
       etd=0;
       do id=1,MDC
        etd=etd+SPCSIG1(is)*part(id,is)*ddir;
       end do
       if (etd.gt.emax)then
           emax=etd;
           isigm=is;
       end if
       if (isigm.gt.0)then
          tpeak_part=2*pi_w/SPCSIG1(isigm);
       else
          tpeak_part=0;
       end if
       end do

!     %%%%%%%%%%%%hs_part%%%%%%%%%%%%%
      etot=0.0;
      do id=1,MDC
          do is=1,MSC
              etot=etot+SPCSIG1(is)**2*part(id,is);
          end do
      end do
      !print *, "etot",etot,"frintf",frintf,"ddir",ddir
      if (etot.gt.1e-8.and.etot.le.1e8)then
          hs_part=4.0*sqrt(etot*frintf*ddir);
      else
      hs_part=0.0
      end if


!     %%%%%%%%%%%%dir_part%%%%%%%%%%%%%%%
      etot=0;
      eex=0;
      eey=0;
      spcdir=0.0;
      do id=1,MDC
          spcdir(id,1)=(id-0.5)*360/MDC;
          spcdir(id,2)=cos(spcdir(id,1)/360*2*pi_w);
          spcdir(id,3)=sin(spcdir(id,1)/360*2*pi_w);
      end do
      do id=1,MDC
          ead=0;
          do is=2,MSC
              ds=SPCSIG1(is)-SPCSIG1(is-1);
              edi=0.5*(SPCSIG1(is)*part(id,is)+SPCSIG1(is-1)*part(id,is-1))*ds;
              ead=ead+edi;
          end do
          eftail=1./(5.-1.);
          ehfr=part(id,MSC)*SPCSIG1(MSC);
          ead=ead+ehfr*SPCSIG1(MSC)*eftail;
          ead=ead*ddir;
          etot=etot+ead;
          eex=eex+ead*spcdir(id,2);
          eey=eey+ead*spcdir(id,3);
      end do
      if (etot.gt.0)then
          dirdeg_part=atan2(eey,eex)*180/pi_w;
          if (dirdeg_part.lt.0)then
              dirdeg_part=dirdeg_part+360;
          end if
      end if
      RETURN
      END SUBROUTINE part_info
!/ ------------------------------------------------------------------- /
!/ End of module W3PARTMD -------------------------------------------- /
!/
#     endif
      END MODULE W3PARTMD
