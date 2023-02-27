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

!==============================================================================|
!     this subroutine is used to calculate the true temperature                !
!     and salinity by including vertical diffusion implicitly.                 !
!==============================================================================|
! NOTE:  This subroutine has a lot of room for optimization.
!        Suggestions: 
!        Remove extra if/then statements inside double do loops
!        Switch loop order to only check if(wetdry) once per node            
!        Remove extra do loops- use array assignment!

! REMOVED DEPTH CHECK FOR NONE WET DRY CASE. THIS IS A LEGACY FROM
! ECOM-SI WHICH HAS LAND IN THE DOMAIN. IT IS NOT NEEDED IN FVCOM

SUBROUTINE VDIF_TS(NCON1,F)                

  !------------------------------------------------------------------------------|

  USE ALL_VARS
  USE MOD_UTILS
  USE BCS
#  if defined (WET_DRY)
  USE MOD_WD
#  endif

# if defined (THIN_DAM)
  USE MOD_DAM,only : NODE_DAM1_N,NODE_DAM2_N,NODE_DAM3_N&
                      &,I_NODE_DAM1_N,I_NODE_DAM2_N,I_NODE_DAM3_N,KDAM
# endif

  IMPLICIT NONE
  INTEGER :: I,K,NCON1,J,KI,hutemp
  !   REAL(SP) :: TMP,TMP1,TMP2,TMP3,QTMP,GW,ZDEP,FKH,UMOLPR,WETFAC
  !   REAL(SP), DIMENSION(0:MT,KB)  :: F
  !   REAL(SP), DIMENSION(M,KB)     :: FF,AF,CF,VHF,VHPF,RAD
  !   REAL(SP), DIMENSION(M)        :: KHBOTTOM,WFSURF,SWRADF
  REAL(DP) :: TMP,TMP1,TMP2,TMP3,QTMP,GW,ZDEP,FKH,UMOLPR,WETFAC
  REAL(SP), DIMENSION(0:MT,KB)  :: F
  REAL(DP), DIMENSION(M,KB)     :: FF,AF,CF,VHF,VHPF,RAD
  REAL(DP), DIMENSION(M)        :: KHBOTTOM,WFSURF,SWRADF,SASURF
  REAL(DP), DIMENSION(M)        :: COSGAMA1,COSGAMA2

# if defined (THIN_DAM)
  REAL(DP) :: ftmp,stmp
# endif

  REAL(SP) :: WFTMP1, WFTMP2, WFTMP3

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: vdif_ts :", NCON1

  UMOLPR = UMOL*1.E0_SP
  SASURF = 0.0_SP
  GW = 0.0_DP
  !
  !------------------------------------------------------------------------------!
  !                                                                              !
  !        the following section solves the equation                             !
  !         dti*(kh*f')'-f=-fb                                                   !
  !                                                                              !
  !------------------------------------------------------------------------------!

  DO K = 2, KBM1
     DO I = 1, M
#  if !defined (WET_DRY)
        !       IF (D(I) > 0.0_SP) THEN
#  else
        IF(ISWETN(I) == 1)THEN
#  endif
           FKH = KH(I,K)

           IF(K == KBM1) THEN
              KHBOTTOM(I)=FKH
           END IF

           AF(I,K-1)=-DTI*(FKH+UMOLPR)/(DZ(I,K-1)*DZZ(I,K-1)*D(I)*D(I))
           CF(I,K)=-DTI*(FKH+UMOLPR)/(DZ(I,K)*DZZ(I,K-1)*D(I)*D(I))
# if defined(WET_DRY)
        END IF
# endif
     END DO
  END DO


  !------------------------------------------------------------------------------!
  !     the net heat flux input.                                                 !
  !     the method shown below can be used when we turn off the                  !
  !     body force in subroutine advt. be sure this method could                 !
  !     cause the surface overheated if the vertical resolution                  !
  !     is not high enough.                                                      !
  !------------------------------------------------------------------------------!

  SELECT CASE(NCON1)
  CASE(1)

     RAD    = 0.0_SP
     WFSURF = 0.0_SP
     SWRADF = 0.0_SP
     IF(HEATING_TYPE == 'flux')THEN
        !     DO I=1,M
        !       WFSURF(I)=WTSURF(I)
        !       SWRADF(I)=SWRAD(I)
        !       SASURF(I)=0.0_SP
        !     END DO
        WFSURF(1:M)=WTSURF(1:M)
        SWRADF(1:M)=SWRAD(1:M)
        SASURF=0.0_SP

#    if !defined (ICE)
!
!-----------------------------------------------------------------------
!  If net heat flux is cooling and SST is at freezing point or below
!  then suppress further cooling. Note: net heat flux sign convention is that 
!  positive (in FVCOM negitive - J. Qi) means heating the ocean (J Wilkin - ROMS).
!-----------------------------------------------------------------------
! 
!  Below the surface heat flux WFSURF is ZERO if cooling AND
!  the SST is cooler than the threshold.  The value is retained if
!  warming.
!
!    WFTMP3 = 0      if SST warmer than threshold (WFTMP1) - change nothing
!    WFTMP3 = 1      if SST colder than threshold (WFTMP1)
!
!    0.5*(WFTMP2-ABS(WFTMP2)) = 0                        if flux is warming in ROMS
!    0.5*(WFTMP2+ABS(WFTMP2)) = 0                        if flux is warming
!                             = WFSURF                   if flux is cooling 
!
        WFTMP1=-2.0_SP              ! nominal SST threshold to cease cooling
        DO I=1,M
          WFTMP2=WFSURF(I)
          WFTMP3=0.5_SP*(1.0_SP+SIGN(1.0_SP,WFTMP1-F(I,1)))
          WFSURF(I) = WFTMP2-WFTMP3*0.5_SP*(WFTMP2+ABS(WFTMP2))
        END DO
#    endif

        DO K = 1, KB
           DO I = 1, M
#    if !defined (WET_DRY)
              !         IF(D(I) > 0.0_SP) THEN
#    else
              IF(ISWETN(I) == 1)THEN
#    endif
                 ZDEP = Z(I,K)*D(I)
                 RAD(I,K)=SWRADF(I)*(RHEAT*EXP(ZDEP/ZETA1)+(1-RHEAT)*EXP(ZDEP/ZETA2))
#    if defined (WET_DRY)
              END IF
# endif
           END DO
        END DO
     END IF

  CASE(2)

     SWRADF   = 0.0_SP
     WFSURF   =0.0_SP
     COSGAMA1 =1.0_SP
!---Considering the salinity conservation, the sea surface salinity flux-----!
!---can be set as zero, that is----------------------------------------------!
     SASURF   = 0.0_SP
     RAD      =0.0_SP

  CASE DEFAULT
     CALL FATAL_ERROR('NCON NOT CORRECT IN VDIF_TS')
  END SELECT


  !------------------------------------------------------------------------------!
  !   surface bcs; wfsurf                                                        !
  !------------------------------------------------------------------------------!

  DO I = 1, M
#  if !defined (WET_DRY)
     !     IF (D(I) > 0.0_SP) THEN
#  else
     IF(ISWETN(I) == 1)THEN
#  endif
        VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
        VHPF(I,1) = -DTI *(SASURF(I)+WFSURF(I)-SWRADF(I) &
             +RAD(I,1)-RAD(I,2)) / (-DZ(I,1)*D(I)) - F(I,1)
        VHPF(I,1) = VHPF(I,1) / (AF(I,1)-1.)
#  if defined (WET_DRY)
     END IF
# endif
  END DO

  DO K = 2, KBM2
     DO I = 1, M
#  if !defined (WET_DRY)
        !     IF (D(I) > 0.0_SP) THEN
#  else
        IF(ISWETN(I) == 1)THEN
#  endif
           VHPF(I,K)=1./ (AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
           VHF(I,K) = AF(I,K) * VHPF(I,K)
           VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-DBLE(F(I,K)) &
                +DTI*(RAD(I,K)-RAD(I,K+1))/(D(I)*DZ(I,K)))*VHPF(I,K)
#  if defined (WET_DRY)
        END IF
# endif
     END DO
  END DO


!!$   DO  K = 1, KBM1
!!$     DO  I = 1, M
!!$#  if !defined (WET_DRY)
!!$     IF (D(I) > 0.0_SP) THEN
!!$#  else
!!$     IF(ISWETN(I) == 1)THEN
!!$#  endif
!!$         FF(I,K) = F(I,K)
!!$       END IF
!!$     END DO
!!$   END DO


#  if !defined (WET_DRY)
  FF(1:M,1:KBM1) = F(1:M,1:KBM1)
#  else
  DO  K = 1, KBM1
     DO  I = 1, M
        IF(ISWETN(I) == 1)THEN
           FF(I,K) = F(I,K)
        END IF
     END DO
  END DO
#  endif

  ! THIS PIECE OF CODE DESPERATELY NEEDS TO BE CLARIFIED
  ! AND STREAMLINED

  DO I = 1, M
     IF (ISONB(I) /= 2) THEN
#  if defined (WET_DRY)
        IF(ISWETN(I) == 1)THEN
#  endif
           TMP1=PFPXB(I)*COS(SITA_GD(I))+PFPYB(I)*SIN(SITA_GD(I))
           TMP2=AH_BOTTOM(I)*PHPN(I)
           TMP3=KHBOTTOM(I)+UMOLPR+AH_BOTTOM(I)*PHPN(I)*PHPN(I)
           TMP=TMP1*TMP2/TMP3*(KHBOTTOM(I)+UMOLPR)

           ! -------------------------------------------------------------------
   IF(NOFLUX_BOT_CONDITION)THEN
     SELECT CASE(NCON1)
     CASE(1)
       IF (TMP1 > 0.0_SP) TMP=0.0_SP
     CASE(2)
!!$       IF (TMP1 > 0.0_SP) TMP=0.0_SP
       IF (TMP1 < 0.0_SP) TMP=0.0_SP
       TMP = -TMP
     CASE DEFAULT
       CALL FATAL_ERROR('NCON NOT CORRECT IN VDIF_TS')
     END SELECT
   ELSE	   
           TMP = 0.0_SP
   END IF	   
           ! -------------------------------------------------------------------

!!$ THIS IS FOR AN OLDER VERSION OF GROUNDWATER
!!$           GW=0.0_SP
!!$           IF(NCON1 == 2) THEN ! MODIFIED FOR NEW BFWDIS
!!$              IF (BFWDIS(I) /= 0.0_SP) THEN
!!$
!!$                 QTMP=-(F(I,KBM1)*D(I)*DZ(I,KBM1)*BFWDIS(I))/ &
!!$                      (D(I)*DZ(I,KBM1)*ART1(I)+BFWDIS(I))
!!$                 GW=DTI/D(I)/DZ(I,KBM1)*QTMP
!!$                 
!!$                 TMP=0.0_SP
!!$
!!$              END IF
!!$           END IF


           FF(I,KBM1) = (CF(I,KBM1)*VHPF(I,KBM2)-FF(I,KBM1)-GW+DTI*(RAD(I,KBM1)-RAD(I,KB)-TMP)/(D(I)*DZ(I,KBM1))) &
                /(CF(I,KBM1)*(1._SP-VHF(I,KBM2))-1._SP)


#  if defined (WET_DRY)
        END IF
# endif
     END IF
  END DO

  DO  K = 2, KBM1
     KI = KB - K
     DO  I = 1, M
        IF(ISONB(I) /= 2) THEN
#  if defined (WET_DRY)
           IF(ISWETN(I) == 1)THEN
#  endif
              FF(I,KI) = (VHF(I,KI)*FF(I,KI+1)+VHPF(I,KI))
#  if defined (WET_DRY)
           END IF
# endif
        END IF


     END DO
  END DO

  DO I = 1, M
#  if defined (WET_DRY)
     IF(ISWETN(I)*ISWETNT(I) == 1 )then
#  endif
        DO K = 1, KBM1
           F(I,K) = FF(I,K)
        END DO
#  if defined (WET_DRY)
     END IF
#  endif
  END DO

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

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: vdif_ts"
  RETURN
END SUBROUTINE VDIF_TS
   !==============================================================================|
