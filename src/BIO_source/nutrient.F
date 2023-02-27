SUBROUTINE NUTRIENT
!==========================================================================
!  This routine computes source and sink terms of phytoplankton in 1D   ===
!  without vertical mixing                                              ===
!==========================================================================
!  Define global data.
   USE MOD_1D
   USE MOD_NUTRIENT
   USE MOD_PHYTOPLANKTON
   USE MOD_ZOOPLANKTON
   USE MOD_DETRITUS
   USE MOD_BACTERIA
   USE MOD_DOM
   
   IMPLICIT NONE
   INTEGER  :: I,K,J,I1,I2
   REAL(SPP) :: SOURCE
!=======================================================================
!  BEGINNING EXECUTABLE
!=======================================================================
!***    ACTIVE RESPIRATION ELEMENTAL RATIO, LANDRY et al., 1993   ***
      DO J=1,NNN
        DO I1=1,NNZ
          DO I2=1,NNP
            N2CZP(J,I2,I1)=(N2CP(J,I2)  &
              -EFFIP(I2,I1)*N2CZ(J,I1))/(1.-EFFIP(I2,I1))
            IF(N2CZ(J,I1).EQ.0) N2CZP(J,I2,I1)=0.0
          END DO
          DO I2=1,NNB
            N2CZB(J,I2,I1)=(N2CB(J,I2)  & 
              -EFFIB(I2,I1)*N2CZ(J,I1))/(1.-EFFIB(I2,I1))
          END DO          
          DO I2=1,NNZ
            N2CZZ(J,I2,I1)=(N2CZ(J,I2)  & 
              -EFFIZ(I2,I1)*N2CZ(J,I1))/(1.-EFFIZ(I2,I1))
          END DO
          DO I2=1,NND
            N2CZD(J,I2,I1)=(N2CD(J,I2)  & 
              -EFFID(I2,I1)*N2CZ(J,I1))/(1.-EFFID(I2,I1))
          END DO
        END DO
      END DO
!*********************  NUTRIENT SOURCES AND SINKS   **********************
      DO 10 J=1,NNN
        DO 20 K=1,KBVM1
!        DO 20 J=1,NNN
          SOURCE=0.
          DO I1=1,NNP
            IF(BIO_NAME(INN+J-1,1).EQ.'SiO3') SOURCE=SOURCE+N2CP(J,I1) &
              *(P_DOM(K,I1)+D_DOM(I1)*U_P(K,I1)) !PASSIVE+ACTIVE EXUDATION 
            SOURCE=SOURCE+N2CP(J,I1)*P_N(K,I1)              &!PHYTO RESPIRATION
                  -UPTAKE_PN(K,J,I1)                         !PHYTO UPTAKE
          END DO 
          DO I1=1,NNZ
            SOURCE=SOURCE+N2CZ(J,I1)*Z_N(K,I1)               !ZOO RESPIRATION
          END DO
          DO I1=1,NNB
            SOURCE=SOURCE+N2CB(J,I1)*RESPIRATION_B(K,I1)     &!BACTERIA RESPIRATION 
                  -UPTAKE_BN(K,J,I1)*EFFIN(J,I1)              !BACTERIA UPTAKE
            DO I2=1,NNM
              SOURCE=SOURCE+N2CDOM(J,I2) &
                *UPTAKE_DOM(K,I2,I1)*(1-EFFIDOM(I2,I1))      !BAC ACTIVE RESPIRATION
            END DO
            DO I2=1,NND
              SOURCE=SOURCE+N2CD(J,I2) &
                *UPTAKE_BD(K,I2,I1)*(1-EFFIBD(I2,I1))        !BAC ACTIVE RESPIRATION
            END DO
          END DO
          DO I1=1,NND
            SOURCE=SOURCE+N2CD(J,I1)*REMINERALIZATION(K,I1)  !DETRITUS REMINERALIZATION
          END DO
          DO I1=1,NNZ
            DO I2=1,NNP
              SOURCE=SOURCE+ACTIVE_R(I1)*N2CZP(J,I2,I1)     &!ACTIVE RESPIRATION
                    *G_P(K,I2,I1)
            END DO
            DO I2=1,NND
              SOURCE=SOURCE+ACTIVE_R(I1)*N2CZD(J,I2,I1)     &!ACTIVE RESPIRATION
                    *G_D(K,I2,I1)
            END DO
            DO I2=1,NNZ
              SOURCE=SOURCE+ACTIVE_R(I1)*N2CZZ(J,I2,I1)     &!ACTIVE RESPIRATION
                    *G_Z(K,I2,I1)
            END DO
            DO I2=1,NNB
              SOURCE=SOURCE+ACTIVE_R(I1)*N2CZB(J,I2,I1)     &!ACTIVE RESPIRATION
                    *G_B(K,I2,I1)
            END DO
          END DO
          BIO_N(K,J)=BIO_N(K,J)+SOURCE*T_STEP                !UPDATE
 20     CONTINUE   !KBVM1
          DO I1=1,NND
             BIO_N(KBVM1,J)=BIO_N(KBVM1,J)+ N2CD(J,I1)*REMINERALIZATION(KBV,I1)*T_STEP   !BENTIHC REMINERALIZATION
          END DO
 10   CONTINUE   !NNN
          IF (NO3_ON) THEN
             DO K=1,KBVM1
                BIO_N(K,1)=BIO_N(K,1)-NH4_NO3(K)*T_STEP
                BIO_N(K,2)=BIO_N(K,2)+NH4_NO3(K)*T_STEP
             END DO
          END IF
      RETURN
END SUBROUTINE NUTRIENT
