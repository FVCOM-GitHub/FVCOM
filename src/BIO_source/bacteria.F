SUBROUTINE BACTERIA
!==========================================================================
!  This routine computes source and sink terms of bacteria in 1D        ===
!  without vertical mixing                                              ===
!==========================================================================
!
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
   REAL(SPP):: TEMPP, TRANS,TRANS1,DOC,DON,N2CDOMM,EFFIDOM_MEAN,PSI_DON_NH4

!=======================================================================
!  BEGINNING EXECUTABLE
!=======================================================================
        DO 10 K=1,KBVM1  !vertical
!******  Determine availability of all biological pool *********
      TEMPP=T_BIO(K)
      DO J=1,NNB
        U_B(J)=UBMAX(J)*EXP(-A_TB(J)*ABS(TEMPP-T_OPTB(J)))*BIO_B(K,J) ! Temperature forcing
        IF (DOCDON) THEN
          DOC=SIGMA_DOM(1,J)*FVDOM(K,1)
          DON=SIGMA_DOM(2,J)*FVDOM(K,2)
          N2CDOMM=DON/DOC              !N:C ratio in DOM
          EFFIDOM_MEAN=0.
          DO I1=1,NNM
            EFFIDOM_MEAN=EFFIDOM_MEAN+EFFIDOM(I1,J)
          END DO
          EFFIDOM_MEAN=EFFIDOM_MEAN/NNM
          TRANS=EFFIDOM_MEAN*N2CB(1,J)/(EFFIN(1,J)*N2CDOMM)
          IF (TRANS > 1.) THEN
            DELTA_B(J)=TRANS-1
          ELSE
            DELTA_B(J)=0.
          END IF
          PSI_DON_NH4=MIN(FVNN(K,1),DELTA_B(J)*DON)
          UPTAKE_BN(K,1,J)=U_B(J)*PSI_DON_NH4/(1+DON+PSI_DON_NH4)
          UPTAKE_DOM(K,2,J)=U_B(J)*DON/(1+DON+PSI_DON_NH4)
          UPTAKE_DOM(K,1,J)=UPTAKE_DOM(K,2,J)*N2CDOMM
        ELSEIF (PO4_ON) THEN
          TRANS=1.
          DO I1=1,NNM
            TRANS=TRANS+SIGMA_DOM(I1,J)*FVDOM(K,I1)
          END DO
          DO I1=1,NND
            TRANS=TRANS+SIGMA_BD(I1,J)*FVD(K,I1)
          END DO
          DO I1=1,NNN
            IF (N2CB(I1,J).GT.0.) TRANS=TRANS+SIGMA_N(I1,J)*FVNN(K,I1)/N2CB(I1,J)
          END DO
          DO I1=1,NNM
            UPTAKE_DOM(K,I1,J)=U_B(J)*SIGMA_DOM(I1,J)*FVDOM(K,I1)/TRANS
          END DO
          DO I1=1,NND
            UPTAKE_BD(K,I1,J)=U_B(J)*SIGMA_BD(I1,J)*FVD(K,I1)/TRANS
          END DO
          DO I1=1,NNN
!           IF (N2CB(I1,J).GT.0.)UPTAKE_BN(K,I1,J)=U_B(J)*SIGMA_N(I1,J)*FVNN(K,I1)/TRANS/N2CB(I1,J)
           IF (N2CB(I1,J).GT.0.)UPTAKE_BN(K,I1,J)=U_B(J)*SIGMA_N(I1,J)*FVNN(K,I1)/TRANS
          END DO
        ELSE    !DOC or DON only
          TRANS=1.
          TRANS1=0.
          DO I1=1,NNM
            TRANS=TRANS+SIGMA_DOM(I1,J)*FVDOM(K,I1)*N2CDOM(1,I1)
            TRANS1=TRANS1+FVDOM(K,I1)*N2CDOM(1,I1)
          END DO
          PSI_DON_NH4=MIN(SIGMA_N(1,J)*FVNN(K,1),DELTA_B(J)*TRANS1*EFFIDOM(1,J))
          TRANS=TRANS+PSI_DON_NH4
          UPTAKE_BN(K,1,J)=U_B(J)*PSI_DON_NH4/TRANS
          DO I1=1,NNM
            UPTAKE_DOM(K,I1,J)=U_B(J)*SIGMA_DOM(I1,J)*FVDOM(K,I1)*N2CDOM(1,I1)/TRANS            
          END DO
        END IF  !IF (DOCDON) 
        RESPIRATION_B(K,J)=R_B(J)*FVB(K,J)*EXP(RP_T*(TEMPP-T_OPTB(J)))
!**********************   BACTERIA SOURCE AND SINK TERMS   ***************
        BIO_B(K,J)=BIO_B(K,J)-RESPIRATION_B(K,J)*T_STEP        !RESPIRATION
        DO I1=1,NNM           
          BIO_B(K,J)=BIO_B(K,J) &                              !DOM UPTAKE
                  +UPTAKE_DOM(K,I1,J)*EFFIDOM(I1,J)*T_STEP
        END DO
        DO I1=1,NNN
          IF(N2CB(I1,J).GT.0.0) THEN
             BIO_B(K,J)=BIO_B(K,J) &                           !DOM UPTAKE
                  +UPTAKE_BN(K,I1,J)*EFFIN(J,I1)*T_STEP/N2CB(I1,J)
          END IF
        END DO
        DO I1=1,NND           
          BIO_B(K,J)=BIO_B(K,J) &                              !DOM UPTAKE
                  +UPTAKE_BD(K,I1,J)*EFFIBD(I1,J)*T_STEP
        END DO
        DO I1=1,NNZ
          BIO_B(K,J)=BIO_B(K,J) &                              !GRAZING LOSSES
                  -G_B(K,J,I1)*T_STEP
        END DO
      END DO !J=1,NNB
 10   CONTINUE
      RETURN
END SUBROUTINE BACTERIA
