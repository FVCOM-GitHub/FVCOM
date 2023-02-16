SUBROUTINE DOM
!==========================================================================
!  This routine computes source and sink terms of DOM in 1D             ===
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
   REAL(SPP) :: SOURCE 
!=======================================================================
!  BEGINNING EXECUTABLE
!=======================================================================

      DO 10 K=1,KBVM1
        DO J=1,NNM
          AGE_DOM(J)=ALPHA_DOM(J)*FVDOM(K,J)
        END DO 
        DO J=1,NNM
          SOURCE=0.
          SOURCE=SOURCE+AGE_DOM(J-1)-AGE_DOM(J)             !AGEING GAIN
          DO I1=1,NNP
            SOURCE=SOURCE+ALPHA_PDOM(J,I1)*P_DOM(K,I1)     &!PASSIVE EXUDATION
                  +ALPHA_PDOM(J,I1)*D_DOM(I1)*U_P(K,I1)     !ACTIVE EXUDATION
          END DO
          DO I1=1,NND
            SOURCE=SOURCE+ALPHA_DDOM(J,I1)*DISSOLUTION(K,I1)
          END DO
          DO I1=1,NNZ
            DO I2=1,NNP
              SOURCE=SOURCE+ALPHA_ZPDOM(J,I2,I1)*G_P(K,I2,I1)
            END DO
            DO I2=1,NND
              SOURCE=SOURCE+ALPHA_ZDDOM(J,I2,I1)*G_D(K,I2,I1)
            END DO
            DO I2=1,NNZ
              SOURCE=SOURCE+ALPHA_ZZDOM(J,I2,I1)*G_Z(K,I2,I1)
            END DO
            DO I2=1,NNB
              SOURCE=SOURCE+ALPHA_ZBDOM(J,I2,I1)*G_B(K,I2,I1)
            END DO         
          END DO !I1=1,NNZ
          DO I1=1,NNB
            SOURCE=SOURCE-UPTAKE_DOM(K,J,I1)
          END DO
          BIO_DOM(K,J)=BIO_DOM(K,J)+SOURCE*T_STEP
        END DO
 10   CONTINUE
      RETURN
END SUBROUTINE DOM
