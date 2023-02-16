SUBROUTINE ZOOPLANKTON
!==========================================================================
!                                                                       ===
!  This routine computes source and sink terms of zooplankton in 1D     ===
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
   INTEGER :: I,K,J,I1,I2
   REAL(SPP) :: SOURCE, TEMPP, TRANS, R_TOTAL, R_TOTAL2 
   REAL(SPP),DIMENSION(NTT) :: SINKK
   REAL(SPP) :: RATIOT

!=======================================================================
!  BEGINNING EXECUTABLE
!=======================================================================
        do 10 k=1,KBV
!******  Determine availability of all biological pool *********
       DO J=1,NNN
        FVNN(K,J)=BIO_N(K,J)-N_0(J)
        IF (FVNN(K,J).LT.0.0) FVNN(K,J)=0.0
       END DO 
       DO J=1, NNP
       FVP(K,J)=BIO_P(K,J)-P_0(J)
        IF (FVP(K,J).LT.0.0) FVP(K,J)=0.0
       END DO
       DO J=1,NNZ
        FVZ(K,J)=BIO_Z(K,J)-Z_0(J)
        IF (FVZ(K,J).LT.0.0) FVZ(K,J)=0.0
       END DO
       DO J=1,NND
        FVD(K,J)=BIO_D(K,J)-D_0(J)
        IF (FVD(K,J).LT.0.0) FVD(K,J)=0.0
       END DO
       DO J=1,NNB
        FVB(K,J)=BIO_B(K,J)-B_0(J)
        IF (FVB(K,J).LT.0.0) FVB(K,J)=0.0
       END DO
       DO J=1,NNM
        FVDOM(K,J)=BIO_DOM(K,J)-DOM_0(J)
        IF (FVDOM(K,J).LT.0.0) FVDOM(K,J)=0.0
       END DO
!*************       ZOOPLANKTON SINKS       *********************
      TEMPP=T_BIO(K)
      DO J=1,NNZ
        TRANS=EXP(RP_T*(TEMPP-T_OPTZ(J)))          !Temporal for the lake model
        Z_D(K,J)=MZD(J)*FVZ(K,J)**M_Z(J)           !Mortality to detritus

!       Z_D(K,J)=MZD(J)*FVZ(K,J)**M_Z(J)*(A_TZ(J))**((TEMPP-T_OPTZ(J))/10.)  !Ji
       	
        Z_N(K,J)=R_Z(J)*FVZ(K,J)*EXP(RP_T*TEMPP)   !Basic respiration to nutrient
!JQI        Z_N(K,J)=R_Z(J)*( FVZ(K,J)**M_Z(J) )*EXP(RP_T*TEMPP)
!JQIJQIif(k == 1)write(499,*) Z_D(K,J),Z_N(K,J)	
      END DO

!*****************     ZOOPLANKTON GRAZING     ********************
      SELECT CASE (TRIM(G_FUNCTION))
        CASE('IVLE1_G')                                  !Ivlev, 1955
        DO 30 I1=1,NNZ
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
          *(A_TZ(I1))**((TEMPP-T_OPTZ(I1))/10.)
!              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))                  
          DO I2=1,NNP                                    ! Grazing on phyto
            G_P(K,I2,I1)=G_RATE(I1)*(1-EXP(-SIGMA_P(I2,I1)*FVP(K,I2)))
          END DO
          DO I2=1,NNZ                                    ! Predation on zoo
            G_Z(K,I2,I1)=G_RATE(I1)*(1-EXP(-SIGMA_Z(I2,I1)*FVZ(K,I2)))
          END DO
          DO I2=1,NND                                    ! Grazing on detritus
            G_D(K,I2,I1)=G_RATE(I1)*(1-EXP(-SIGMA_D(I2,I1)*FVD(K,I2)))
          END DO
          DO I2=1,NNB                                    ! Predation on bacteria
            G_B(K,I2,I1)=G_RATE(I1)*(1-EXP(-SIGMA_B(I2,I1)*FVB(K,I2)))
          END DO
 30   CONTINUE
        CASE('RECTI_G')                                  !Armstrong et al., 1999
        DO 40 I1=1,NNZ
          R_TOTAL=0.                                     ! Total food
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  

          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+SIGMA_P(I2,I1)*FVP(K,I2)
          END DO
          DO I2=1,NNZ                                    
            R_TOTAL=R_TOTAL+SIGMA_Z(I2,I1)*FVZ(K,I2)
          END DO
          DO I2=1,NND                                    
            R_TOTAL=R_TOTAL+SIGMA_D(I2,I1)*FVD(K,I2)
          END DO
          DO I2=1,NNB 
            R_TOTAL=R_TOTAL+SIGMA_B(I2,I1)*FVB(K,I2)
          END DO
          
          IF (R_TOTAL <= K_ZG(I1)) THEN              
            DO I2=1,NNP                                    ! Grazing on phyto
              G_P(K,I2,I1)=G_RATE(I1)*SIGMA_P(I2,I1)*FVP(K,I2)/K_ZG(I1)
            END DO
            DO I2=1,NNZ                                    ! Predation on zoo
              G_Z(K,I2,I1)=G_RATE(I1)*SIGMA_Z(I2,I1)*FVZ(K,I2)/K_ZG(I1)
            END DO
            DO I2=1,NND                                    ! Grazing on detritus
              G_D(K,I2,I1)=G_RATE(I1)*SIGMA_D(I2,I1)*FVD(K,I2)/K_ZG(I1)
            END DO
            DO I2=1,NNB                                    ! Predation on bacteria
              G_B(K,I2,I1)=G_RATE(I1)*SIGMA_B(I2,I1)*FVB(K,I2)/K_ZG(I1)
            END DO
          ELSE
            DO I2=1,NNP                                    ! Grazing on phyto
              G_P(K,I2,I1)=G_RATE(I1)*SIGMA_P(I2,I1)*FVP(K,I2)/R_TOTAL
            END DO
            DO I2=1,NNZ                                    ! Predation on zoo
              G_Z(K,I2,I1)=G_RATE(I1)*SIGMA_Z(I2,I1)*FVZ(K,I2)/R_TOTAL
            END DO
            DO I2=1,NND                                    ! Grazing on detritus
              G_D(K,I2,I1)=G_RATE(I1)*SIGMA_D(I2,I1)*FVD(K,I2)/R_TOTAL
            END DO
            DO I2=1,NNB                                    ! Predation on bacteria
              G_B(K,I2,I1)=G_RATE(I1)*SIGMA_B(I2,I1)*FVB(K,I2)/R_TOTAL
            END DO
          END IF
 40   CONTINUE          
        CASE('CLI_G')                                    !Leonard et al., 1999
        DO 50 I1=1,NNZ
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))                  
          DO I2=1,NNP                                    ! Grazing on phyto
            G_P(K,I2,I1)=G_RATE(I1)*SIGMA_P(I2,I1)*FVP(K,I2) &
                        *EXP(-SIGMA_P(I2,I1)*FVP(K,I2))
          END DO
          DO I2=1,NNZ                                    ! Predation on zoo
            G_Z(K,I2,I1)=G_RATE(I1)*SIGMA_Z(I2,I1)*FVZ(K,I2) &
                        *EXP(-SIGMA_Z(I2,I1)*FVZ(K,I2))
          END DO
          DO I2=1,NND                                    ! Grazing on detritus
            G_D(K,I2,I1)=G_RATE(I1)*SIGMA_D(I2,I1)*FVD(K,I2) &
                        *EXP(-SIGMA_D(I2,I1)*FVD(K,I2))
          END DO
          DO I2=1,NNB                                    ! Predation on bacteria
            G_B(K,I2,I1)=G_RATE(I1)*SIGMA_B(I2,I1)*FVB(K,I2) &
                        *EXP(-SIGMA_B(I2,I1)*FVB(K,I2))
          END DO
 50   CONTINUE
        CASE('IVLE2_G')                                  ! Hofmann & Amber, 88
        DO 60 I1=1,NNZ
          R_TOTAL=0.                                     ! Total food
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  

          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+SIGMA_P(I2,I1)*FVP(K,I2)
          END DO
          DO I2=1,NNZ                                    
            R_TOTAL=R_TOTAL+SIGMA_Z(I2,I1)*FVZ(K,I2)
          END DO
          DO I2=1,NND                                    
            R_TOTAL=R_TOTAL+SIGMA_D(I2,I1)*FVD(K,I2)
          END DO
          DO I2=1,NNB 
            R_TOTAL=R_TOTAL+SIGMA_B(I2,I1)*FVB(K,I2)
          END DO
                       
            DO I2=1,NNP                                    ! Grazing on phyto
              G_P(K,I2,I1)=G_RATE(I1)*SIGMA_P(I2,I1)*FVP(K,I2)/R_TOTAL &
                          *(1-EXP(-SIGMA_P(I2,I1)*FVP(K,I2)))
            END DO
            DO I2=1,NNZ                                    ! Predation on zoo
              G_Z(K,I2,I1)=G_RATE(I1)*SIGMA_Z(I2,I1)*FVZ(K,I2)/R_TOTAL &
                          *(1-EXP(-SIGMA_Z(I2,I1)*FVZ(K,I2)))
            END DO
            DO I2=1,NND                                    ! Grazing on detritus
              G_D(K,I2,I1)=G_RATE(I1)*SIGMA_D(I2,I1)*FVD(K,I2)/R_TOTAL &
                          *(1-EXP(-SIGMA_D(I2,I1)*FVD(K,I2)))
            END DO
            DO I2=1,NNB                                    ! Predation on bacteria
              G_B(K,I2,I1)=G_RATE(I1)*SIGMA_B(I2,I1)*FVB(K,I2)/R_TOTAL &
                          *(1-EXP(-SIGMA_B(I2,I1)*FVB(K,I2)))
            END DO
 60   CONTINUE
        CASE('MM1_G')                                    ! Molonet & Field, 1991
        DO 70 I1=1,NNZ
          R_TOTAL=0.                                     ! Total food
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  
          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+SIGMA_P(I2,I1)*FVP(K,I2)
          END DO
          DO I2=1,NNZ                                    
            R_TOTAL=R_TOTAL+SIGMA_Z(I2,I1)*FVZ(K,I2)
          END DO
          DO I2=1,NND                                    
            R_TOTAL=R_TOTAL+SIGMA_D(I2,I1)*FVD(K,I2)
          END DO
          DO I2=1,NNB 
            R_TOTAL=R_TOTAL+SIGMA_B(I2,I1)*FVB(K,I2)
          END DO
          
            DO I2=1,NNP                                    ! Grazing on phyto
              G_P(K,I2,I1)=G_RATE(I1)*SIGMA_P(I2,I1)*FVP(K,I2) &
              /(K_ZG(I1)+R_TOTAL)
            END DO
            DO I2=1,NNZ                                    ! Predation on zoo
              G_Z(K,I2,I1)=G_RATE(I1)*SIGMA_Z(I2,I1)*FVZ(K,I2) &
                          /(K_ZG(I1)+R_TOTAL)
            END DO
            DO I2=1,NND                                    ! Grazing on detritus
              G_D(K,I2,I1)=G_RATE(I1)*SIGMA_D(I2,I1)*FVD(K,I2) &
                          /(K_ZG(I1)+R_TOTAL)
            END DO
            DO I2=1,NNB                                    ! Predation on bacteria
              G_B(K,I2,I1)=G_RATE(I1)*SIGMA_B(I2,I1)*FVB(K,I2) &
                          /(K_ZG(I1)+R_TOTAL)
            END DO
 70   CONTINUE   
         CASE('DM_G')                                    ! Denman 1991
        DO 75 I1=1,NNZ
          R_TOTAL=0.                                     ! Total food
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))
!          *(A_TZ(I1))**((TEMPP-T_OPTZ(I1))/10.)      !Ji 

!              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  


          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+FVP(K,I2)**2
          END DO
 !         DO I2=1,NNZ                                    
 !           R_TOTAL=R_TOTAL+FVZ(K,I2)**2
 !         END DO
 !         DO I2=1,NND                                    
 !           R_TOTAL=R_TOTAL+FVD(K,I2)**2
 !         END DO
 !         DO I2=1,NNB 
 !           R_TOTAL=R_TOTAL+FVB(K,I2)**2
 !         END DO
          
            DO I2=1,NNP                                    ! Grazing on phyto
              G_P(K,I2,I1)=G_RATE(I1)*FVP(K,I2)**2 &
              /(K_ZG(I1)**2+R_TOTAL)
!JQIJQIif(k == 1)write(200,*) A_TZ(I1),TEMPP,T_OPTZ(I1),G_P(K,I2,I1),FVP(K,I2)
!JQIJQIif(k == 1)write(201,*) BIO_Z(K,I1),G_P(K,I2,I1),G_MAX(I1),FVP(K,I2),(A_TZ(I1))**((TEMPP-T_OPTZ(I1))/10.)
            END DO
 !           DO I2=1,NNZ                                    ! Predation on zoo
 !             G_Z(K,I2,I1)=G_RATE(I1)*FVZ(K,I2)**2 &
 !                         /(K_ZG(I1)**2+R_TOTAL)
 !           END DO
 !           DO I2=1,NND                                    ! Grazing on detritus
 !             G_D(K,I2,I1)=G_RATE(I1)*FVD(K,I2)**2 &
 !                         /(K_ZG(I1)**2+R_TOTAL)
 !           END DO
 !           DO I2=1,NNB                                    ! Predation on bacteria
 !             G_B(K,I2,I1)=G_RATE(I1)*FVB(K,I2)**2 &
 !                         /(K_ZG(I1)**2+R_TOTAL)
 !           END DO
 75   CONTINUE 
	            
        CASE('MM2_G')                                    ! Evans 1988
        DO 80 I1=1,NNZ
          R_TOTAL=0.                                     ! Total food
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  

          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+SIGMA_P(I2,I1)*FVP(K,I2)
          END DO
          DO I2=1,NNZ                                    
            R_TOTAL=R_TOTAL+SIGMA_Z(I2,I1)*FVZ(K,I2)
          END DO
          DO I2=1,NND                                    
            R_TOTAL=R_TOTAL+SIGMA_D(I2,I1)*FVD(K,I2)
          END DO
          DO I2=1,NNB 
            R_TOTAL=R_TOTAL+SIGMA_B(I2,I1)*FVB(K,I2)
          END DO
          
          DO I2=1,NNP                                    ! Grazing on phyto
            G_P(K,I2,I1)=G_RATE(I1)*SIGMA_P(I2,I1)*FVP(K,I2)/(K_ZG(I1)+R_TOTAL) &
              *(R_TOTAL-P_C(I1))/(R_TOTAL+K_ZG(I1)-P_C(I1))
          END DO
          DO I2=1,NNZ                                    ! Predation on zoo
            G_Z(K,I2,I1)=G_RATE(I1)*SIGMA_Z(I2,I1)*FVZ(K,I2)/(K_ZG(I1)+R_TOTAL) &
              *(R_TOTAL-P_C(I1))/(R_TOTAL+K_ZG(I1)-P_C(I1))
          END DO
          DO I2=1,NND                                    ! Grazing on detritus
            G_D(K,I2,I1)=G_RATE(I1)*SIGMA_D(I2,I1)*FVD(K,I2)/(K_ZG(I1)+R_TOTAL) &
              *(R_TOTAL-P_C(I1))/(R_TOTAL+K_ZG(I1)-P_C(I1))
          END DO
          DO I2=1,NNB                                    ! Predation on bacteria
            G_B(K,I2,I1)=G_RATE(I1)*SIGMA_B(I2,I1)*FVB(K,I2)/(K_ZG(I1)+R_TOTAL) &
              *(R_TOTAL-P_C(I1))/(R_TOTAL+K_ZG(I1)-P_C(I1))
          END DO
 80   CONTINUE          
        CASE('MM3_G')                                    ! Verity, 1991
        DO 90 I1=1,NNZ
          R_TOTAL=1.0                                    ! Total food
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  

          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+SIGMA_P(I2,I1)*FVP(K,I2)
          END DO
          DO I2=1,NNZ                                    
            R_TOTAL=R_TOTAL+SIGMA_Z(I2,I1)*FVZ(K,I2)
          END DO
          DO I2=1,NND                                    
            R_TOTAL=R_TOTAL+SIGMA_D(I2,I1)*FVD(K,I2)
          END DO
          DO I2=1,NNB 
            R_TOTAL=R_TOTAL+SIGMA_B(I2,I1)*FVB(K,I2)
          END DO
          
          DO I2=1,NNP                                    ! Grazing on phyto
            G_P(K,I2,I1)=G_RATE(I1)*SIGMA_P(I2,I1)*FVP(K,I2)/R_TOTAL
          END DO
          DO I2=1,NNZ                                    ! Predation on zoo
            G_Z(K,I2,I1)=G_RATE(I1)*SIGMA_Z(I2,I1)*FVZ(K,I2)/R_TOTAL
          END DO
          DO I2=1,NND                                    ! Grazing on detritus
            G_D(K,I2,I1)=G_RATE(I1)*SIGMA_D(I2,I1)*FVD(K,I2)/R_TOTAL
          END DO
          DO I2=1,NNB                                    ! Predation on bacteria
            G_B(K,I2,I1)=G_RATE(I1)*SIGMA_B(I2,I1)*FVB(K,I2)/R_TOTAL
          END DO
 90   CONTINUE          
        CASE('SMM_G')                                    ! Fasham et al., 1990
        DO 100 I1=1,NNZ
          R_TOTAL=0.                                     ! Total food
          R_TOTAL2=0.
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  

          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+SIGMA_P(I2,I1)*FVP(K,I2)
            R_TOTAL2=R_TOTAL+SIGMA_P(I2,I1)*FVP(K,I2)**2
          END DO
          DO I2=1,NNZ                                    
            R_TOTAL=R_TOTAL+SIGMA_Z(I2,I1)*FVZ(K,I2)
            R_TOTAL2=R_TOTAL+SIGMA_Z(I2,I1)*FVZ(K,I2)**2     
          END DO
          DO I2=1,NND                                    
            R_TOTAL=R_TOTAL+SIGMA_D(I2,I1)*FVD(K,I2)
            R_TOTAL2=R_TOTAL+SIGMA_D(I2,I1)*FVD(K,I2)**2
          END DO
          DO I2=1,NNB 
            R_TOTAL=R_TOTAL+SIGMA_B(I2,I1)*FVB(K,I2)
            R_TOTAL2=R_TOTAL+SIGMA_B(I2,I1)*FVB(K,I2)**2
          END DO
          
          DO I2=1,NNP                                    ! Grazing on phyto
            G_P(K,I2,I1)=G_RATE(I1)*SIGMA_P(I2,I1)*FVP(K,I2)**2 &
              /(K_ZG(I1)*R_TOTAL*R_TOTAL2)
          END DO
          DO I2=1,NNZ                                    ! Predation on zoo
            G_Z(K,I2,I1)=G_RATE(I1)*SIGMA_Z(I2,I1)*FVZ(K,I2)**2 &
              /(K_ZG(I1)*R_TOTAL*R_TOTAL2)
          END DO
          DO I2=1,NND                                    ! Grazing on detritus
            G_D(K,I2,I1)=G_RATE(I1)*SIGMA_D(I2,I1)*FVD(K,I2)**2 &
              /(K_ZG(I1)*R_TOTAL*R_TOTAL2)
          END DO
          DO I2=1,NNB                                    ! Predation on bacteria
            G_B(K,I2,I1)=G_RATE(I1)*SIGMA_B(I2,I1)*FVB(K,I2)**2 &
              /(K_ZG(I1)*R_TOTAL*R_TOTAL2)
          END DO
 100  CONTINUE 
        CASE('GSF1_G')                                   ! Tansky, 1978
        DO 110 I1=1,NNZ
          R_TOTAL=0.                                     ! Total food
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  

          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+(SIGMA_P(I2,I1)*FVP(K,I2))**M_G(I1)
          END DO
          DO I2=1,NNZ                                    
            R_TOTAL=R_TOTAL+(SIGMA_Z(I2,I1)*FVZ(K,I2))**M_G(I1)
          END DO
          DO I2=1,NND                                    
            R_TOTAL=R_TOTAL+(SIGMA_D(I2,I1)*FVD(K,I2))**M_G(I1)
          END DO
          DO I2=1,NNB 
            R_TOTAL=R_TOTAL+(SIGMA_B(I2,I1)*FVB(K,I2))**M_G(I1)
          END DO
          
          DO I2=1,NNP                                    ! Grazing on phyto
            G_P(K,I2,I1)=G_RATE(I1)*(SIGMA_P(I2,I1)*FVP(K,I2)) &
                        **(M_G(I1))/R_TOTAL
          END DO
          DO I2=1,NNZ                                    ! Predation on zoo
            G_Z(K,I2,I1)=G_RATE(I1)*(SIGMA_Z(I2,I1)*FVZ(K,I2)) &
                        **(M_G(I1))/R_TOTAL
          END DO
          DO I2=1,NND                                    ! Grazing on detritus
            G_D(K,I2,I1)=G_RATE(I1)*(SIGMA_D(I2,I1)*FVD(K,I2)) &
                        **(M_G(I1))/R_TOTAL
          END DO
          DO I2=1,NNB                                    ! Predation on bacteria
            G_B(K,I2,I1)=G_RATE(I1)*(SIGMA_B(I2,I1)*FVB(K,I2)) &
                        **(M_G(I1))/R_TOTAL
          END DO
 110  CONTINUE 
        CASE('GSF2_G')                                   ! Vance, 1978
        DO 120 I1=1,NNZ
          R_TOTAL=0.                                     ! Total food
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  

          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+(SIGMA_P(I2,I1)*FVP(K,I2))
          END DO
          DO I2=1,NNZ                                    
            R_TOTAL=R_TOTAL+(SIGMA_Z(I2,I1)*FVZ(K,I2))
          END DO
          DO I2=1,NND                                    
            R_TOTAL=R_TOTAL+(SIGMA_D(I2,I1)*FVD(K,I2))
          END DO
          DO I2=1,NNB 
            R_TOTAL=R_TOTAL+(SIGMA_B(I2,I1)*FVB(K,I2))
          END DO
          
          DO I2=1,NNP                                    ! Grazing on phyto
            G_P(K,I2,I1)=G_RATE(I1)*(SIGMA_P(I2,I1)*FVP(K,I2)) &
                        **(M_G(I1))/(R_TOTAL**M_G(I1))
          END DO
          DO I2=1,NNZ                                    ! Predation on zoo
            G_Z(K,I2,I1)=G_RATE(I1)*(SIGMA_Z(I2,I1)*FVZ(K,I2)) &
                        **(M_G(I1))/(R_TOTAL**M_G(I1))
          END DO
          DO I2=1,NND                                    ! Grazing on detritus
            G_D(K,I2,I1)=G_RATE(I1)*(SIGMA_D(I2,I1)*FVD(K,I2)) &
                        **(M_G(I1))/(R_TOTAL**M_G(I1))
          END DO
          DO I2=1,NNB                                    ! Predation on bacteria
            G_B(K,I2,I1)=G_RATE(I1)*(SIGMA_B(I2,I1)*FVB(K,I2)) &
                        **(M_G(I1))/(R_TOTAL**M_G(I1))
          END DO
 120  CONTINUE 
        CASE('GSMM_G')                                   ! Gismervik&Andersen, 97
        DO 130 I1=1,NNZ
          R_TOTAL=1.                                     ! Total food
          G_RATE(I1)=G_MAX(I1)*BIO_Z(K,I1)  &            ! T forcing and biomass
              *EXP(-A_TZ(I1)*ABS(TEMPP-T_OPTZ(I1)))  

          DO I2=1,NNP                                    ! Total food
            R_TOTAL=R_TOTAL+SIGMA_P(I2,I1)*FVP(K,I2)**M_G(I1)
          END DO
          DO I2=1,NNZ                                    
            R_TOTAL=R_TOTAL+SIGMA_Z(I2,I1)*FVZ(K,I2)**M_G(I1)
          END DO
          DO I2=1,NND                                    
            R_TOTAL=R_TOTAL+SIGMA_D(I2,I1)*FVD(K,I2)**M_G(I1)
          END DO
          DO I2=1,NNB 
            R_TOTAL=R_TOTAL+SIGMA_B(I2,I1)*FVB(K,I2)**M_G(I1)
          END DO
          
          DO I2=1,NNP                                    ! Grazing on phyto
            G_P(K,I2,I1)=G_RATE(I1)*SIGMA_P(I2,I1)*FVP(K,I2) &
                        **M_G(I1)/R_TOTAL
          END DO
          DO I2=1,NNZ                                    ! Predation on zoo
            G_Z(K,I2,I1)=G_RATE(I1)*SIGMA_Z(I2,I1)*FVZ(K,I2) &
                        **M_G(I1)/R_TOTAL
          END DO
          DO I2=1,NND                                    ! Grazing on detritus
            G_D(K,I2,I1)=G_RATE(I1)*SIGMA_D(I2,I1)*FVD(K,I2) &
                        **M_G(I1)/R_TOTAL
          END DO
          DO I2=1,NNB                                    ! Predation on bacteria
            G_B(K,I2,I1)=G_RATE(I1)*SIGMA_B(I2,I1)*FVB(K,I2) &
                        **M_G(I1)/R_TOTAL
          END DO
 130  CONTINUE      
      END SELECT 
      GO TO 8    !No check in Ji's program
!************   CHECK FOOD AVAILABILITY   *******
      DO I1=1,NNP
      RATIOT=1.
!        SINKK(I1)= P_DOM(K,I1)+P_D(K,I1)+P_N(K,I1)
        SINKK(I1)= DPDOM(I1)*FVP(K,I1) + MPD(I1)*FVP(K,I1)**M_P(I1) + R_P(I1)*FVP(K,I1)*EXP(RP_T*TEMPP)
        DO I2=1,NNZ
          SINKK(I1)=SINKK(I1)+G_P(K,I1,I2)
        END DO
        SINKK(I1)=SINKK(I1)*T_STEP         
        RATIOT=FVP(K,I1)/(SINKK(I1)+1.E-30)
!JQI        IF(SINKK(I1).LT.0.0001) RATIOT=0.0
        IF(RATIOT.LT.1.) THEN
          DO I2=1,NNZ
            G_P(K,I1,I2)=G_P(K,I1,I2)*RATIOT
          END DO
        END IF
      END DO

      DO I1=1,NND
        RATIOT=1.
        SINKK(I1)=0.
        DO I2=1,NNZ
          SINKK(I1)=SINKK(I1)+G_D(K,I1,I2)
        END DO
        SINKK(I1)=SINKK(I1)*T_STEP
        RATIOT=FVD(K,I1)/(SINKK(I1)+1.E-30)
        IF(RATIOT.LT.1. .and. SINKK(I1).GT. 0.0) THEN
          DO I2=1,NNZ
            G_D(K,I1,I2)=G_D(K,I1,I2)*RATIOT
          END DO
        END IF
      END DO

      DO I1=1,NNZ
        RATIOT=1.
        SINKK(I1)= Z_D(K,I1)+Z_N(K,I1)
        DO I2=1,NNZ
          SINKK(I1)=SINKK(I1)+G_Z(K,I1,I2)
        END DO
        SINKK(I1)=SINKK(I1)*T_STEP
        RATIOT=FVZ(K,I1)/(SINKK(I1)+1.E-30)
        IF(RATIOT.LT.1.0 .and. SINKK(I1).GT. 0.0) THEN
          DO I2=1,NNZ
            G_Z(K,I1,I2)=G_Z(K,I1,I2)*RATIOT
          END DO
          Z_D(K,I1)=Z_D(K,I1)*RATIOT
          Z_N(K,I1)=Z_N(K,I1)*RATIOT
        END IF
      END DO
 8    CONTINUE
!********************  ZOOPLANKTON SOURCE AND SINKS   ***********************
      DO J=1,NNZ
        SOURCE=-Z_D(K,J)-Z_N(K,J)                   !Losses
        DO I1=1,NNP        
          SOURCE=SOURCE+EFFIP(I1,J)*G_P(K,I1,J)     !Grazing on phyto
!          FVP(K,I1)=FVP(K,I1)-G_P(K,I1,J)*T_STEP
        END DO
        DO I1=1,NND        
          SOURCE=SOURCE+EFFID(I1,J)*G_D(K,I1,J)     !Grazing on detritus
!          FVD(K,I1)=FVD(K,I1)-G_D(K,I1,J)*T_STEP
        END DO
        DO I1=1,NNB        
          SOURCE=SOURCE+EFFIB(I1,J)*G_B(K,I1,J)     !Grazing on bacteria
!          FVB(K,I1)=FVB(K,I1)-G_B(K,I1,J)*T_STEP
        END DO 
        DO I1=1,NNZ        
          SOURCE=SOURCE+EFFIZ(I1,J)*G_Z(K,I1,J)   & !Predation on zoo
                       -G_Z(K,J,I1)                 !Predation losses 
!          FVZ(K,I1)=FVZ(K,I1)-G_Z(K,I1,J)*T_STEP
        END DO 
        RECRUIT_LOSS(J)=R_RECRUIT(J)*SOURCE  
        BIO_Z(K,J)=BIO_Z(K,J)+T_STEP*SOURCE       
      END DO !J=1,NNZ
      DO J=1,NNZ
        BIO_Z(K,J)=BIO_Z(K,J)+T_STEP*(            &
                  +RECRUIT_LOSS(J-1)              & !Recruitment gain
                  -RECRUIT_LOSS(J))                 !Recruitment loss
      END DO
 10   CONTINUE  !depth
      RETURN
END SUBROUTINE ZOOPLANKTON
