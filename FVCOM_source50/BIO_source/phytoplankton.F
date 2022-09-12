SUBROUTINE PHYTOPLANKTON(index_id,i_id,iint_id)
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
   INTEGER :: I,K,J,I1,I2
   REAL(SPP) :: TEMPP, TRANS, TATANUL, IRRAD, IRRAD1
   REAL(SPP) :: SINKK(10),RATIO(10),RATIOT,SMALL
   
   integer :: index_id,i_id,iint_id

!=======================================================================
!  BEGINNING EXECUTABLE
!=======================================================================
      SMALL = 1.E-30

        DO 10 K=1,KBVM1  !vertical
!**************  PAR-PHYTOPLANKTON PRODUCTION    **************************
      TEMPP=T_BIO(K)
      TATANUL=1.0E-10                                  !Total attanuation
      TATANUL=TATANUL+ATANU_W*DELTA_D(K)               !Water attanuation
      TRANS=0.0_SPP
      DO J=1,NNP
        TRANS=TRANS+BIO_P(K,J)*CHL2C(J)                !Total chlorophyll
      END DO
      IF(TRANS.GT.0.0) TATANUL=TATANUL+TRANS*ATANU_C*DELTA_D(K)
      TRANS=0.0_SPP
      DO J=1,NND
        TRANS=TRANS+BIO_D(K,J)                         !Total detritus
      END DO
      IF(TRANS.GT.0.0) TATANUL=TATANUL+TRANS*ATANU_D*DELTA_D(K)     ! Detrital attanuation
      IRRAD1=IRRAD0*EXP(-TATANUL)                      ! Irradiance at the layer bottom
      IRRAD=(IRRAD0-IRRAD1)/TATANUL                    ! Average irradiance 
      DO 20 J=1,NNP
        UTR(J)=EXP(-A_TP(J)*ABS(TEMPP-T_OPTP(J)))      ! Temperature forcing on phyto
     ! 	UTR(J)=A_TP(J)**((TEMPP-T_OPTP(J))/10.)
      SELECT CASE (TRIM(L_FUNCTION))                   ! Light function choice
        CASE('EXP_LIGHT')                              ! Exponential (Frank et al. ,86)
!          ULR(J)=EXP(TATANUL*DEPTH_Z(K))
          ULR(J)=EXP(ATANU_W*DEPTH_Z(K))
        CASE('SL62_LIGHT')                             ! Exponential (Steele, 62)
          ULR(J)=IRRAD/I_OPT(J)*EXP(1-IRRAD/I_OPT(J))
        CASE('MM_LIGHT')                               ! Michaealis (Baly 35)
          ULR(J)=ALPHAP(J)*IRRAD/(K_LIGHT(J)+ALPHAP(J)*IRRAD)
        CASE('LB_LIGHT')                               ! Michaelis (Bannister, 79)
          ULR(J)=ALPHAP(J)*IRRAD/((K_LIGHT(J))**N_P(J) &
                +(ALPHAP(J)*IRRAD)**BETAP(J))**(1/N_P(J))
        CASE('V65_LIGHT')                              ! Michaelis (Vollenweider 65)
          ULR(J)=ALPHAP(J)*IRRAD/(SQRT(I_OPT(J)**2+(ALPHAP(J)*IRRAD)**2) &
                *(1+(BETAP(J)*IRRAD/I_OPT(J))**2)**(N_P(J)/2))
        CASE('PE78_LIGHT')                             ! Michaelis (Peeters & Eiler, 78)
          ULR(J)=IRRAD/I_OPT(J)*(2+ALPHAP(J))/(1+ALPHAP(J)*IRRAD/I_OPT(J) &
                +(IRRAD/I_OPT(J))**2)
        CASE('WNS74_LIGHT')                            ! Exponentail (Webb, 1978)
          ULR(J)=1-EXP(-ALPHAP(J)*IRRAD/(UMAX(J)+SMALL))
        CASE('PGH80_LIGHT')                            ! Exponential (Platt et al., 80)
          ULR(J)=(1-EXP(-ALPHAP(J)*IRRAD/(UMAX(J)+SMALL))) &
                *EXP(-BETAP(J)*IRRAD/(UMAX(J)+SMALL))
	CASE('SH92_LIGHT')                             ! Exponential (STEEL & Henderson, 1992) 
	  ULR(J)=(1-EXP(-ALPHAP(J)*IRRAD)*EXP(-BETAP(J)*IRRAD))
        CASE('JI_LIGHT')                            ! Exponential (Platt et al., 80)
          ULR(J)=(1-EXP(-ALPHAP(J)*IRRAD))*EXP(-BETAP(J)*IRRAD)
        CASE('JP76_LIGHT')                             ! Tangent (Jessy & Platt, 76)
          ULR(J)=TANH(ALPHAP(J)*IRRAD/(UMAX(J)+SMALL))
        CASE('BWDC99_LIGHT')                           ! Tangent (Bissett et al., 99)
          ULR(J)=TANH(ALPHAP(J)*(IRRAD-I_OPT(J))/(UMAX(J)+SMALL)) &
                *EXP(BETAP(J)*(I_OPT(J)-IRRAD))
      END SELECT       
!*************   NUTRIENT-PHYTOPLANKTON GORWTH   ***************
      UNRMIN(J)=1.0_SPP
      DO I1=1,NNN
        UNR(J,I1)=FVNN(K,I1)/(FVNN(K,I1)+KSN(I1,J))
      END DO
      IF (.NOT. NO3_ON) THEN
         DO I1=1,NNN
           IF(N2CP(I1,J).GT.0.) UNRMIN(J)=MIN(UNRMIN(J),UNR(J,I1))
         END DO    
      ELSE
        UNR(J,2)=FVNN(K,2)/(FVNN(K,2)+KSN(2,J)) &   ! NO3 must be nutrient 2
                *KSN(J,1)/(FVNN(K,1)+KSN(1,J))      ! Parker NH4 inhibition
        UNR(J,NNN+1)=UNR(J,1)+UNR(J,2)                ! Total nitrogen factor
         DO I1=3,NNN
           IF(N2CP(I1,J).GT.0.) UNRMIN(J)=MIN(UNRMIN(J),UNR(J,I1)) 
         END DO
         UNRMIN(J)=MIN(UNRMIN(J),UNR(J,NNN+1))        ! Combined NH4 and NO3
      END IF
      U_P(K,J)=UMAX(J)*UTR(J)*(ALPHA_U(J)*ULR(J)*UNRMIN(J) &
              +(1-ALPHA_U(J))*MIN(ULR(J),UNRMIN(J)))       &
              *BIO_P(K,J)*CHL2C(J)                            !Chl concentration
!JQIJQI
!      if(index_id == 1 .and. mod(iint_id,30) == 0)then
!        if(i_id == 153)write(101,'(2i6,i3,9e12.5)')   iint_id,i_id,k,u_p(k,j),utr(j),ulr(j),unrmin(j),bio_p(k,j),irrad,TATANUL,irrad0,irrad1
!        if(i_id == 5681)write(102,'(2i6,i3,9e12.5)')  iint_id,i_id,k,u_p(k,j),utr(j),ulr(j),unrmin(j),bio_p(k,j),irrad,TATANUL,irrad0,irrad1
!        if(i_id == 10244)write(103,'(2i6,i3,9e12.5)') iint_id,i_id,k,u_p(k,j),utr(j),ulr(j),unrmin(j),bio_p(k,j),irrad,TATANUL,irrad0,irrad1
!      end if	

      DO I1=1,NNN
        UPTAKE_PN(K,I1,J)=U_P(K,J)*N2CP(I1,J)                    !Nutrient uptake
      END DO
      IF (NO3_ON) THEN
        if(unr(j,nnn+1).ne.0) then
          UPTAKE_PN(K,1,J)=U_P(K,J)*N2CP(1,J)*UNR(J,1)/UNR(J,NNN+1)    !NH4
          UPTAKE_PN(K,2,J)=U_P(K,J)*N2CP(1,J)*UNR(J,2)/UNR(J,NNN+1)    !NO3
        else
          uptake_pn(k,1,j)=0.
          uptake_pn(k,2,j)=0.
        endif
      END IF
!***************     SINK TERMS    ****************
      P_DOM(K,J)=DPDOM(J)*FVP(K,J)                       !DOM EXUDATION
      P_D(K,J)=UTR(J)*MPD(J)*FVP(K,J)**M_P(J)
      P_N(K,J)=R_P(J)*FVP(K,J)*EXP(RP_T*TEMPP)         !Phyto respiration with T forcing
     ! RATIOT=FVP(K,J)/((P_DOM(K,J)+P_D(K,J)+P_N(K,J))*T_STEP+1.E-30)
     ! IF(RATIOT.LT.1.) THEN
     !   P_DOM(K,J)=0.0
     !   P_D(K,J)=0.0
     !   P_N(K,J)=0.0
     ! END IF
!      FVP(K,J)=FVP(K,J)-(P_DOM(K,J)+P_D(K,J)+P_N(K,J))*T_STEP
 20   CONTINUE
      GO TO 21  !TO CONFORM JI'S
!**************     CHECK NUTRIENT AVAILABILITY  ***********
      RATIOT=1.0
      SINKK=0.
      RATIO=1.0
      DO I1=1,NNN
        DO I2=1,NNP
           SINKK(I1)=SINKK(I1)+UPTAKE_PN(K,I1,I2)
        END DO
        SINKK(I1)=SINKK(I1)*T_STEP
      END DO
      DO I1=1,NNN
        IF (SINKK(I1).GT.FVNN(K,I1).AND.SINKK(I1).GT.0.0) THEN
          RATIO(I1)=FVNN(K,I1)/(SINKK(I1)+1.E-30)
          RATIOT=MIN(RATIOT,RATIO(I1))
        END IF
      END DO  
      IF (RATIOT .LT. 1.0) THEN
        DO I2=1,NNP
          U_P(K,I2)=U_P(K,I2)*RATIOT
          DO I1=1,NNN
            UPTAKE_PN(K,I1,I2)=UPTAKE_PN(K,I1,I2)*RATIOT
          END DO
        END DO
      END IF
 21   CONTINUE   
!****************   LIGHT-CONTROLLED NITRIFICATION    *************
      IF (NO3_ON) THEN
        IF (IRRAD0.GT.0.) THEN
           TRANS=L_NH4N-IRRAD1
           IF (TRANS > 0.) THEN
             NH4_NO3(K)=R_AN*TRANS*FVNN(K,1)/L_NH4N
           ELSE
             NH4_NO3(K)=0.
           END IF
        ELSE 
            NH4_NO3(K)=R_AN*FVNN(K,1)
        END IF
     END IF
     
      IRRAD0=IRRAD1     !For next level
 10   CONTINUE
!*******************   CELLS SINKING    ***********
      DO I1=1,NNP             
!         WSNK_P(1,I1)=(-0.5)*W_P(I1)*(FVP(1,I1)+FVP(2,I1))/DELTA_D(1) 
         WSNK_P(1,I1)=(-1)*W_P(I1)*(FVP(1,I1))/DELTA_D(1) 
!         WSNK_P(KBVM1,I1)=W_P(I1)*(FVP(KBV-2,I1)-FVP(KBVM1,I1))/DELTA_D(KBVM1) !Sink out
         WSNK_P(KBVM1,I1)=W_P(I1)*(FVP(KBV-2,I1))/DELTA_D(KBVM1) !No sink out
      END DO 
      DO K=2,KBV-2
        DO I1=1,NNP
!         WSNK_P(K,I1)=W_P(I1)*0.5*(FVP(K-1,I1)-FVP(K+1,I1))/DELTA_D(K)
         WSNK_P(K,I1)=W_P(I1)*(FVP(K-1,I1)-FVP(K,I1))/DELTA_D(K)
        END DO
      END DO
!**************   PHYTOPLANKTON SOURCES AND SINKS   **************
      DO K=1,KBVM1
        DO J=1,NNP
          BIO_P(K,J)=BIO_P(K,J)+((1-D_DOM(J))*U_P(K,J)+WSNK_P(K,J) & !Growth & Sinking 
               -P_DOM(K,J)-P_D(K,J)-P_N(K,J))*T_STEP       !Exudaiton, mortality, respiration
          DO I2=1,NNZ
            BIO_P(K,J)=BIO_P(K,J)-G_P(K,J,I2)*T_STEP       !GRAZING LOSSES
          END DO
        END DO
      END DO
      RETURN
END SUBROUTINE PHYTOPLANKTON
