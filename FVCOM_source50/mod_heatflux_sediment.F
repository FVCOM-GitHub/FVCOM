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

  MODULE MOD_HEATFLUX_SEDIMENT
  
  USE MOD_PREC
  USE MOD_UTILS
  USE CONTROL,ONLY  : casename,NMLUNIT,IINT,DTI,IREPORT
  USE ALL_VARS,ONLY : MGL,M,KBM1,D,MSR,PAR,SERIAL,TF1,WTSURF
  USE MOD_PAR,ONLY  : NGID
  
  IMPLICIT NONE
  
  ! =================== Added code(heat flux at tidal flat) ykcho ==============|   
   REAL(SP), ALLOCATABLE  :: FMST1(:,:) !! MUD TEMPERATURE IN PREVIOUS TIME
   REAL(SP), ALLOCATABLE  :: FMST2(:,:) !! MUD TEMPERATURE IN PRESENT TIME
   REAL(SP), ALLOCATABLE  :: NMST(:,:) !! CALCULATED MUD TEMERATURE
   REAL(SP), ALLOCATABLE  :: WHEAT1(:) !! HEAT FLUX at tidal flat (positive/negative means heat gain/loss at tidal flat) 
   REAL(SP), ALLOCATABLE  :: WHEAT(:)  !! HEAT FLUX at bottom layer of  sea water (positive/negative means heat gain/loss at bottm layer of sea water)
!==============================================================================
  REAL(SP) :: MVHCDEP(8),MTDCDEP(8),MVHC(8),MTDC(8)

!  Parameters in NameList NML_HEATFLUX_SEDIMENT
   LOGICAL HEATFLUX_SEDIMENT_ON
   CHARACTER(LEN=80) MUD_INITIAL_TEMP_FILE
   REAL(SP) :: VOLUMETRIC_HEAT_CAPACITY
   REAL(SP) :: THERMAL_DIFFUSIVITY
   REAL(SP) :: EFFECTIVE_THICKNESS
   REAL(SP) :: CRITICAL_DEPTH
   REAL(SP) :: SEDIMENT_TEMPERATURE
   
   NAMELIST /NML_HEATFLUX_SEDIMENT/      &
        & HEATFLUX_SEDIMENT_ON,          &
	& MUD_INITIAL_TEMP_FILE,         &       !doesn't need
	& VOLUMETRIC_HEAT_CAPACITY,      &
	& THERMAL_DIFFUSIVITY,           &
	& EFFECTIVE_THICKNESS,           &
	& CRITICAL_DEPTH,                &
	& SEDIMENT_TEMPERATURE
	
   CONTAINS	
  !==============================================================================!
  !
  !==============================================================================!
   SUBROUTINE NAME_LIST_INITIALIZE_HEATFLUX_SEDIMENT
   
   IMPLICIT NONE
   
   HEATFLUX_SEDIMENT_ON     = .FALSE.
   MUD_INITIAL_TEMP_FILE    = "none"
   VOLUMETRIC_HEAT_CAPACITY = 3.65e+6_SP     !unit = J/(m**3*K)
   THERMAL_DIFFUSIVITY      = 0.38e-6_SP     !unit = m**2/s
   EFFECTIVE_THICKNESS      = 0.01           !unit = m
   CRITICAL_DEPTH           = 1.0            !unit = m   
   SEDIMENT_TEMPERATURE     = 10.0           !unit = degree
   
   RETURN
   END SUBROUTINE NAME_LIST_INITIALIZE_HEATFLUX_SEDIMENT
  !==============================================================================!
  !
  !==============================================================================!
   SUBROUTINE NAME_LIST_PRINT_HEATFLUX_SEDIMENT
   
   IMPLICIT NONE
   
   WRITE(UNIT=IPT,NML=NML_HEATFLUX_SEDIMENT)
   
   RETURN
   END SUBROUTINE NAME_LIST_PRINT_HEATFLUX_SEDIMENT
  !==============================================================================!
  !
  !==============================================================================!
   SUBROUTINE NAME_LIST_READ_HEATFLUX_SEDIMENT
   
   IMPLICIT NONE

   INTEGER :: ios
   CHARACTER(LEN=120) :: FNAME

   IF(DBG_SET(dbg_sbr)) write(IPT,*) "Subroutine Begins: name_list_read_heatflux_sediment;"

    ios = 0
    FNAME = "./"//trim(casename)//"_run.nml"
    IF(DBG_SET(dbg_io)) write(IPT,*) "Get_nestpar: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    !READ NESTING FLAG
    READ(UNIT=NMLUNIT, NML=NML_HEATFLUX_SEDIMENT,IOSTAT=ios)  
    IF(ios /= 0)THEN
       IF(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_HEATFLUX_SEDIMENT)
       CALL FATAL_ERROR("Can Not Read NameList NML_HEATFLUX_SEDIMENT from file: "//trim(FNAME))
    END IF

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:NML_HEATFLUX_SEDIMENT"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_HEATFLUX_SEDIMENT)

    CLOSE(NMLUNIT)

    MTDC(1) = VOLUMETRIC_HEAT_CAPACITY
    MVHC(1) = THERMAL_DIFFUSIVITY
!    EFFECTIVE_THICKNESS
!    CRITICAL_DEPTH
    FMST2(:,1) = SEDIMENT_TEMPERATURE
   
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Ends: name_list_read_heatflux_sediment;"    

  END SUBROUTINE NAME_LIST_READ_HEATFLUX_SEDIMENT
  !==============================================================================!
  !
  !==============================================================================!
  SUBROUTINE READ_MUD_INITIAL_TEMPERATURE
  
  IMPLICIT NONE  
   
  REAL(SP),ALLOCATABLE :: RTEMP1(:,:),RTEMP2(:,:)
  INTEGER :: I,J,IERR
  
!
! =================== Added code(heat flux at tidal flat) ykcho===============|   
!-------MUD INITIAL TEMPERATURE------------------------------------------------!
!

  ALLOCATE(RTEMP1(MGL,50),RTEMP2(MGL,50))
!  OPEN(111,FILE='MSST.DAT',STATUS='OLD')
  OPEN(111,FILE=trim(MUD_INITIAL_TEMP_FILE),STATUS='OLD')

  IF(MSR)THEN

    DO I=1,MGL
      READ(111,*)(RTEMP1(I,J),J=1,50),(RTEMP2(I,J),J=1,50)
    ENDDO
  ENDIF 
 
  CLOSE(111)

# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    CALL MPI_BCAST(RTEMP1,50*MGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(RTEMP2,50*MGL,  MPI_F,0,MPI_COMM_WORLD,IERR)
  END IF
# endif
  ALLOCATE(FMST1(M,50),FMST2(M,50),NMST(M,50))
  IF(SERIAL)THEN
    FMST1(1:MGL,:) = RTEMP1(1:MGL,:)
    FMST2(1:MGL,:) = RTEMP2(1:MGL,:)
  END IF
# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
    DO I=1,M
      FMST1(I,:) = RTEMP1(NGID(I),:)
      FMST2(I,:) = RTEMP2(NGID(I),:)
    END DO
  END IF
# endif
  DEALLOCATE(RTEMP1,RTEMP2)
  
  
!  mod_main.F:   REAL(SP), ALLOCATABLE  :: FMST1(:,:) !! MUD TEMPERATURE IN PREVIOUS TIME
!  mod_main.F:   REAL(SP), ALLOCATABLE  :: FMST2(:,:) !! MUD TEMPERATURE IN PRESENT TIME

  
  RETURN
  END SUBROUTINE READ_MUD_INITIAL_TEMPERATURE
  !==============================================================================|   
  !
  !==============================================================================!
  SUBROUTINE CALCULATE_HEATFLUX_SEDIMENT
  
  IMPLICIT NONE  

! =================== Added code(heat flux at tidal flat) ykcho===============|   
!  REAL(SP) :: MVHCDEP(8),MTDCDEP(8),MVHC(8),MTDC(8)
  REAL(SP) :: MHCDZ,MWCONT,MHEAT,MHEAT1,MHEAT2,STMUDZ,MUDHC,MUDDF,GRTS,GRTB,ALBEDO
  REAL(SP)  ,ALLOCATABLE :: NMTMP(:,:), WHTMP(:),WHTMP1(:),WTTMP(:)
  INTEGER  :: MHCMDZ,IMUD,I,K
  REAL(SP) :: AAA,BBB,KTEM,KAIRT,KSST,EVAPS,EVAPS1,EVAPS2,EVAPS3,EVAPS4,EVAPA,DHAIRTE,DHRETHU, &
              CLOUDF,DHCLOUD,LONGHEAT,STBOC,EMIS,CONDHEAT,DENAIR,CONHC,ZDIMM
! ============================================================================|   


!
!--- Heat Flux and Short Wave Radiation----------------------------------------!
!
!     DHAIRTE = UFACT*UHAIRTE(L1) + FACT*UHAIRTE(L2)
!     DHRETHU = UFACT*UHRETHU(L1) + FACT*UHRETHU(L2)
!     DHAIRPR = UFACT*UHAIRPR(L1) + FACT*UHAIRPR(L2)
!     DHCLOUD = UFACT*UHCLOUD(L1) + FACT*UHCLOUD(L2)
!          TX = UFACT*UWIND(L1)   + FACT*UWIND(L2)
!          TY = UFACT*VWIND(L1)   + FACT*VWIND(L2)
!     WINDSP=(TX**2.+TY**2.)**0.5
!     SWRAD(:)= UFACT*UHSHORT(L1) + FACT*UHSHORT(L2)
!     ALBEDO=0.92
!     SWRAD=SWRAD*ALBEDO
!     SPCP=4.2174E3_SP
!     ROSEA = 1.023E3_SP
!     SPRO = SPCP*ROSEA
!     DO I=1,M
!     DSST=T1(I,1)
!     STBOC = 5.67/(10**8)
!     EMIS = 0.96
!     DENAIR = 1.2929
!     CONHC = 1003.0
!     ZDIMM = 0.0014
!     KTEM = 273.16
!     KAIRT = KTEM + DHAIRTE
!     KSST = KTEM + DSST
!     EVAPS1=10.79574*(1.-(KTEM/KAIRT))
!     EVAPS2=5.02800*ALOG10(KAIRT/KTEM)
!     EVAPS3=(1.50475E-4)*(1.-10.**(-8.2969*((KAIRT/KTEM)-1.)))
!     EVAPS4=(0.42873E-3)*(10.**(+4.76955*(1.-(KTEM/KAIRT)))-1.)
!     EVAPS=10.**(EVAPS1-EVAPS2+EVAPS3+EVAPS4+0.78614)
!     EVAPA=DHRETHU*EVAPS/100.
!     CLOUDF=DHCLOUD/10.
!     LONGHEAT=STBOC*EMIS*(KAIRT**4)*(0.4-(0.05*(EVAPA**0.5)))
!     LONGHEAT=LONGHEAT+4*STBOC*EMIS*(KAIRT**3)*(KSST-KAIRT)
!     LONGHEAT=LONGHEAT*(1-0.75*CLOUDF**3.4)
!     CONDHEAT=DENAIR*CONHC*ZDIMM*(1.+WINDSP)*(KSST-KAIRT)
!     LVEV=(2500.84-2.35*DSST)*(10.**3)
!     ABH=(0.621*EVAPA)/(DHAIRPR-(1.-0.621)*EVAPA)
!     EVAPS1=10.79574*(1.-(KTEM/KSST))
!     EVAPS2=5.02800*ALOG10(KSST/KTEM)
!     EVAPS3=(1.50475E-4)*(1.-10.**(-8.2969*((KSST/KTEM)-1.)))
!     EVAPS4=(0.42873E-3)*(10.**(+4.76955*(1.-(KTEM/KSST)))-1.)
!     EVAPS=10.**(EVAPS1-EVAPS2+EVAPS3+EVAPS4+0.78614)
!     SPH=(0.621*EVAPS)/(DHAIRPR-(1.-0.621)*EVAPS)
!     EVAPHEAT=DENAIR*LVEV*ZDIMM*(1.+WINDSP)*(SPH-ABH)
!     WTSURF(I)=SWRAD(I)-LONGHEAT-CONDHEAT-EVAPHEAT
!     END DO
! =================== Added code(heat flux at tidal flat) ykcho===============|   
!----------TIDAL FLAT HEAT FLUX -------------------------------------------------     
!  MVHCDEP(1)  =  0.
!  MVHCDEP(2)  =  0.03
!  MVHCDEP(3)  =  0.07
!  MVHCDEP(4)  =  0.10
!  MVHCDEP(5)  =  0.15
!  MVHCDEP(6)  =  0.25
!  MVHCDEP(7)  =  0.45
!  MVHCDEP(8)  =  1.00
!  MTDCDEP     =  MVHCDEP
!  MVHC(1)     =  3.84
!  MVHC(2)     =  3.65
!  MVHC(3)     =  3.48
!  MVHC(4)     =  3.31
!  MVHC(5)     =  3.14
!  MVHC(6)     =  2.96
!  MVHC(7)     =  2.65
!  MVHC(8)     =  2.50
!  MTDC(1)     =  0.30
!  MTDC(2)     =  0.40
!  MTDC(3)     =  0.50
!  MTDC(4)     =  0.60
!  MTDC(5)     =  0.70
!  MTDC(6)     =  0.80
!  MTDC(7)     =  0.90
!  MTDC(8)     =  1.00
!  MHCDZ       =  0.02
!  MHCMDZ      =  50
!  MWCONT      =  0.7
   
  WHEAT  =0.
  WHEAT1 =0.
  DO I=1,M
    STMUDZ      =  0.1
    MHEAT  =0.
    MHEAT1 =0.
    MHEAT2 =0.
!    DO IMUD=2,MHCMDZ-1
!      DO K=1,7
!        IF(STMUDZ > MVHCDEP(K) .AND. STMUDZ <= MVHCDEP(K+1))THEN
!          AAA=(MVHC(K)-MVHC(K+1))/(MVHCDEP(K)-MVHCDEP(K+1))
!          BBB=MVHC(K)-AAA*MVHCDEP(K)
!          MUDHC=AAA*STMUDZ+BBB
!          AAA=(MTDC(K)-MTDC(K+1))/(MTDCDEP(K)-MTDCDEP(K+1))
!          BBB=MTDC(K)-AAA*MTDCDEP(K)
!          MUDDF=AAA*STMUDZ+BBB
!        END IF
!      ENDDO
!      GRTS = FMST2(I,IMUD-1) -  FMST2(I,IMUD)
!      GRTB = FMST2(I,IMUD)   -  FMST2(I,IMUD+1)
!      NMST(I,IMUD) = (MUDDF*(10.0**(-6.))*((GRTS-GRTB)/(MHCDZ**2)))*DTI
!      NMST(I,IMUD) = FMST2(I,IMUD) + NMST(I,IMUD)
!      MHEAT1=MHEAT1+MUDHC*(10.0**6.)*((NMST(I,IMUD)-FMST1(I,IMUD))/(2.*DTI))*MHCDZ
!      STMUDZ=STMUDZ+MHCDZ
!    END DO

    IF (D(I) > 0.201) THEN
      FMST2(I,1) = SEDIMENT_TEMPERATURE
      MHEAT2=MVHC(1)*MTDC(1)*200.*(TF1(I,KBM1)-FMST2(I,1))
      WHEAT(I)=MHEAT2
      WHEAT1(I)=WHEAT(I)*(-1.)
    ELSE
!      KAIRT = KTEM + DHAIRTE
!      KSST  = KTEM + FMST2(I,1)
!      EVAPS1=10.79574*(1.-(KTEM/KAIRT))
!      EVAPS2=5.02800*ALOG10(KAIRT/KTEM)
!      EVAPS3=(1.50475E-4)*(1.-10.**(-8.2969*((KAIRT/KTEM)-1.)))
!      EVAPS4=(0.42873E-3)*(10.**(+4.76955*(1.-(KTEM/KAIRT)))-1.)
!      EVAPS=10.**(EVAPS1-EVAPS2+EVAPS3+EVAPS4+0.78614)
!      EVAPA=DHRETHU*EVAPS/100.
!      CLOUDF=DHCLOUD/10.
!      LONGHEAT=STBOC*EMIS*(KAIRT**4)*(0.4-(0.05*(EVAPA**0.5)))
!      LONGHEAT=LONGHEAT+4*STBOC*EMIS*(KAIRT**3)*(KSST-KAIRT)
!      LONGHEAT=LONGHEAT*(1-0.75*CLOUDF**3.4)
!      CONDHEAT=DENAIR*CONHC*ZDIMM*(1.+WINDSP)*(KSST-KAIRT)
!      LVEV=(2500.84-2.35*FMST2(I,1))*(10.**3)
!      ABH=(0.621*EVAPA)/(DHAIRPR-(1.-0.621)*EVAPA)
!      EVAPS1=10.79574*(1.-(KTEM/KSST))
!      EVAPS2=5.02800*ALOG10(KSST/KTEM)
!      EVAPS3=(1.50475E-4)*(1.-10.**(-8.2969*((KSST/KTEM)-1.)))
!      EVAPS4=(0.42873E-3)*(10.**(+4.76955*(1.-(KTEM/KSST)))-1.)
!      EVAPS=10.**(EVAPS1-EVAPS2+EVAPS3+EVAPS4+0.78614)
!      SPH=(0.621*EVAPS)/(DHAIRPR-(1.-0.621)*EVAPS)
!      EVAPHEAT=DENAIR*LVEV*ZDIMM*(1.+WINDSP)*(SPH-ABH)*MWCONT
!      MHEAT2=SWRAD(I)*(ALBEDO*2.-0.99)-LONGHEAT-CONDHEAT-EVAPHEAT
!      WHEAT(I)=MHEAT2
    END IF
!    NMST(I,1)=FMST1(I,1)+(2.*DTI)*((MHEAT2-MHEAT1)/(MVHC(1)*(10.0**6.)*MHCDZ))
!    NMST(I,1)=NMST(I,1)*0.5+FMST2(I,1)*0.5
!-------------------initial tidal sediment temperature------------------------------          
!    NMST(I,MHCMDZ)=15.-COS((IINT*DTI/86400.+141.-51.)/365.*3.141592)*7.
!--------------------initial tidal sediment temperature------------------------------          
  ENDDO  
!  IF(SERIAL)THEN 
!    IF(MOD(IINT,IREPORT) == 0)THEN
!      WRITE(CMTO,'(I5.5)')IINT/IREPORT
!      MOUTF='MUD_SM'//CMTO//'.DAT'
!      OPEN(2,FILE=MOUTF,STATUS='UNKNOWN')
!      DO I=1,MGL
!        WRITE(2,'(4f10.5)')NMST(I,1),WHEAT(I),WHEAT1(I),WTSURF(I)
!      ENDDO
!      CLOSE(2)
!    ENDIF
!  ENDIF
!# if defined (MULTIPROCESSOR)
!  IF(PAR)THEN
!    ALLOCATE(WHTMP(MGL))
!    ALLOCATE(WHTMP1(MGL))
!    ALLOCATE(WTTMP(MGL))
!    ALLOCATE(NMTMP(MGL,50))
!    CALL GATHER(LBOUND(WHEAT,1),  UBOUND(WHEAT,1),  M,MGL, 1,MYID,NPROCS,NMAP,WHEAT,  WHTMP)
!    CALL GATHER(LBOUND(WHEAT1,1), UBOUND(WHEAT1,1), M,MGL, 1,MYID,NPROCS,NMAP,WHEAT1, WHTMP1)
!    CALL GATHER(LBOUND(WTSURF,1), UBOUND(WTSURF,1), M,MGL, 1,MYID,NPROCS,NMAP,WTSURF, WTTMP)
!    CALL GATHER(LBOUND(NMST,1),   UBOUND(NMST,1),   M,MGL,50,MYID,NPROCS,NMAP,NMST,   NMTMP)
!    IF(MSR)THEN
!      IF(MOD(IINT,IREPORT) == 0)THEN
!        WRITE(CMTO,'(I5.5)')IINT/IREPORT
!        MOUTF='MUD_SM'//CMTO//'.DAT'
!        OPEN(2,FILE=MOUTF,STATUS='UNKNOWN')
!        DO I=1,MGL
!          WRITE(2,'(4f10.5)')NMTMP(I,1),WHTMP(I),WHTMP1(I),WTTMP(I)
!        ENDDO
!        CLOSE(2)
!      ENDIF
!    ENDIF
!    DEALLOCATE(WHTMP,WHTMP1,WTTMP,NMTMP)
!  ENDIF
!# endif

!  WTSURF    = -WTSURF/SPRO*RAMP
!  SWRAD     = -SWRAD/SPRO*RAMP*0.
!  WHEAT1    = -WHEAT1/SPRO*RAMP
!  FMST1=FMST2
!  FMST2=NMST
  
  RETURN
  END SUBROUTINE CALCULATE_HEATFLUX_SEdIMENT
! ============================================================================|   

END MODULE MOD_HEATFLUX_SEDIMENT

