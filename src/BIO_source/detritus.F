subroutine detRITUS
!==========================================================================
!  This routine computes source and sink terms of detritus in 1D        ===
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
   INTEGER :: I,K,J,I1,I2, N_SINK
   REAL(DPP) :: DELTA_DMIN,W_D1,SOURCE,TEMPP,RIEN
!=======================================================================
!  BEGINNING EXECUTABLE
!=======================================================================
!*****************   SINKING   ************************
      DO J=1,NND
!         DO k=1,KBVM1
         DO k=1,KBV
           D_SINK(K,J)=0.0_DPP
           DETRITE(K,J)=FVD(K,J)
         END DO
!********  AGGREGATION, DISSOLUTION AND REMINERLIZATION   *********
        DO K=1,KBVM1
          TEMPP=T_BIO(K)
          AGGREGATION(K,J)=ALPHA_DAG(J)*DETRITE(K,J)**2
          DISAGGREGATION(K,J)=ALPHA_DDG(J)*DETRITE(K,J)
          DISSOLUTION(K,J)=D_D(J)*DETRITE(K,J)
 !         REMINERALIZATION(K,J)=D_RN(J)*DETRITE(K,J)
         REMINERALIZATION(K,J)=D_RN(J)*DETRITE(K,J)*EXP(RP_T*TEMPP)
!         REMINERALIZATION(K,J)=D_RN(J)*DETRITE(K,J)*(RP_T)**((TEMPP-T_OPTP(1))/10.)
        END DO
         REMINERALIZATION(KBV,J)=D_RN(J)*BIO_D(KBV,J)*EXP(RP_T*TEMPP)  !BENTHIC LAYER 	 	 
      END DO

      DELTA_DMIN=1.0E+5_SPP    !CFL REQUIRES W*DTI<DELTA_DMIN, THE THINNEAST LAYER
      DO K=1,KBVM1      
        DELTA_DMIN=MIN(DELTA_DMIN,DELTA_D(K))  
      END DO

      DO 10 J=1,NND
        W_D_TOTAL(J)=W_D(J)*T_STEP
        IF (W_D_TOTAL(J) > DELTA_DMIN) THEN
          N_SINK=INT(W_D_TOTAL(J)/DELTA_DMIN)+1
        ELSE
          N_SINK=1
        END IF
        W_D1=W_D_TOTAL(J)/N_SINK

        DO I1=1,N_SINK
!          D_SINK1(1)=-0.5*W_D1*(DETRITE(1,J)+DETRITE(2,J))/DELTA_D(1)    !SURFACE 
          D_SINK1(1)=-1._DPP*W_D1*(DETRITE(1,J))/DELTA_D(1)    !SURFACE         
          DO K=2,KBV-2
!            D_SINK1(K)=0.5*W_D1*(DETRITE(K-1,J)-DETRITE(K+1,J))/DELTA_D(K)
            D_SINK1(K)=W_D1*(DETRITE(K-1,J)-DETRITE(K,J))/DELTA_D(K)
          END DO
          D_SINK1(KBVM1)=W_D1*(DETRITE(KBVM1-1,J)-DETRITE(KBVM1,J))/DELTA_D(KBVM1) !Sink out
!          D_SINK1(KBVM1)=W_D1*(DETRITE(KBVM1-1,J))/DELTA_D(KBVM1) !No sink out 
          D_SINK1(KBV)=W_D1*(DETRITE(KBVM1,J))/DELTA_D(KBVM1) !FLUX TO BENTHIC
!          DO K=1,KBVM1
          DO K=1,KBV
!RIEN=RIEN+D_SINK1(K)*DELTA_D(K)
            DETRITE(K,J)=DETRITE(K,J)+D_SINK1(K)
            D_SINK(K,J)=D_SINK(K,J)+D_SINK1(K)
          END DO
!print*,'de 72 ',rien,D_SINK1(:)
        END DO
        DO K=1,KBV
          D_SINK(K,J)=D_SINK(K,J)/T_STEP
        END DO

 10    CONTINUE
!*****************   DETRITE SOURCES AND SINKS   *******************
      DO 20 J=1,NND
         DO K=1,KBVM1
           SOURCE=0.
           SOURCE=SOURCE                                   &
                     +D_SINK(K,J)                          & ! SINKING
                     +AGGREGATION(K,J-1)                   & 
                     +DISAGGREGATION(K,J+1)                &
                     -AGGREGATION(K,J)                     &
                     -DISAGGREGATION(K,J)                  &
                     -DISSOLUTION(K,J)                     &
                     -REMINERALIZATION(K,J)            
           DO I1=1,NNP
             SOURCE=SOURCE+EPSILON_PD(J,I1)*P_D(K,I1)       !PHYTO MORTALITY
           END DO
           DO I1=1,NNZ
             SOURCE=SOURCE+EPSILON_ZD(J,I1)*Z_D(K,I1)      &!ZOO MORTALITY
                          -G_D(K,J,I1)                      !GRAZING LOSS
           END DO 
           DO I1=1,NNB
           SOURCE=SOURCE-UPTAKE_BD(K,J,I1)                  !BACTERIA CONSUMPTION
           END DO
                    
           DO I1=1,NNZ                                      ! FEEDING LOSS
             DO I2=1,NNP
               SOURCE=SOURCE+ALPHA_PD(J,I2,I1)*G_P(K,I2,I1) 
             END DO
             DO I2=1,NND
               SOURCE=SOURCE+ALPHA_DD(J,I2,I1)*G_D(K,I2,I1) 
             END DO
             DO I2=1,NNZ
               SOURCE=SOURCE+ALPHA_ZD(J,I2,I1)*G_Z(K,I2,I1) 
             END DO
             DO I2=1,NNB
               SOURCE=SOURCE+ALPHA_BD(J,I2,I1)*G_B(K,I2,I1) 
             END DO
           END DO
  	   BIO_D(K,J)=BIO_D(K,J)+SOURCE*T_STEP   
         END DO !K=1,KBVM!
         BIO_D(KBV,J)=BIO_D(KBV,J)+(D_SINK(KBV,J)-REMINERALIZATION(KBV,J))*T_STEP         
 20    CONTINUE
      RETURN
END SUBROUTINE DETRITUS
