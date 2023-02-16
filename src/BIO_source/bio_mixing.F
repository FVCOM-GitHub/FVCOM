SUBROUTINE BIO_MIXING 
!==============================================================================|
!     vertical mixing adopted from FVCOM                                       !
!==============================================================================|
   USE MOD_1D
   USE MOD_PHYTOPLANKTON
   USE MOD_ZOOPLANKTON
   USE MOD_BACTERIA
   USE MOD_DETRITUS
   USE MOD_DOM
   USE MOD_NUTRIENT
   !USE MOD_BIO_3D
   IMPLICIT NONE
   INTEGER :: I,K,J,J1,KI, KVM1,KVM2
   REAL(SPP) :: GW,FKH
   REAL(SPP), DIMENSION(KBV,NTT)     :: F
   REAL(SPP), DIMENSION(KBV,NTT)     :: FF,VHF,VHPF
   REAL(SPP), DIMENSION(KBV)         :: AF,CF 

!------------------------------------------------------------------------------!
!        the following section solves the equation                             !
!         dti*(kh*f')'-f=-fb                                                   !
!------------------------------------------------------------------------------!
   KVM1=KBV-1
   KVM2=KBV-2
   F(1:KBV,1:NTT)=BIO_VAR(1:KBV,1:NTT)

   DO K = 2, KVM1
     FKH=KM_BIO(K)
     AF(K-1)=-T_STEP*FKH/(DELTA_D(K-1)*DELTA_Z(K-1))
     CF(K)=-T_STEP*FKH/(DELTA_D(K)*DELTA_Z(K-1))
   END DO

   DO I=1,NTT
     VHF(1,I)  = AF(1) / (AF(1)-1.)
     VHPF(1,I) = -F(1,I)
     VHPF(1,I) = VHPF(1,I) / (AF(1)-1.)
   END DO !I=1.NTT

   DO I=1,NTT
     DO K = 2, KVM2
       VHPF(K,I) = 1./ (AF(K)+CF(K)*(1.-VHF(K-1,I))-1.)
       VHF(K,I)  = AF(K) * VHPF(K,I)
       VHPF(K,I) = (CF(K)*VHPF(K-1,I)-DBLE(F(K,I)))*VHPF(K,I)
     END DO
   END DO !I=1,NTT

   DO I=1,NTT
     DO  K = 1, KVM1
       FF(K,I) = F(K,I) 
     END DO
   END DO

   DO I=1,NTT  !BOTTOM KBV-1
     GW=0.0_SPP
     FF(KVM1,I) = ((CF(KVM1)*VHPF(KVM2,I)-FF(KVM1,I)-GW) &
                  /(CF(KVM1)*(1.-VHF(KVM2,I))-1._SPP)) !GW IS GROUND FLUX
   END DO

   DO I=1,NTT
     DO  K = 2, KVM1
       KI = KBV - K
       FF(KI,I) = (VHF(KI,I)*FF(KI+1,I)+VHPF(KI,I))
     END DO
   END DO

   DO I = 1, NTT
     DO K = 1, KVM1
       BIO_VAR(K,I) = FF(K,I)
!      if (BIO_VAR(K,I) < 0.0_SPP) BIO_VAR(K,I)=0.0001	 
     END DO
   END DO

   RETURN
END SUBROUTINE BIO_MIXING
