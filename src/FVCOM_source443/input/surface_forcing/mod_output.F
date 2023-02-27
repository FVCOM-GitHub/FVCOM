MODULE binout 
  USE MOD_SURFACEFORCE
  USE MOD_BINARY


CONTAINS



! -------------------------------------------------------------------  
  SUBROUTINE DUMP_HFX(NOW)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP) :: hour, mean
    INTEGER :: STATUS, I

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_HFX "

    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,SWR,SWRGL)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,NHF,NHFGL)
# endif
    END IF

    IF (.NOT. MSR) RETURN

    hour = DAYS(NOW - ZEROTIME) * 24.0_SP

    write(ipt,*) "time", hour

    mean = sum(swrgl(1:mGL))/real(NGL,SP)
    if(isnan(mean)) call fatal_error("swrgl is nan")
    write(ipt,*) "SWRGL: min/max/mean",minval(swrgl(1:MGL)),maxval(swrgl(1:MGL)),mean

    mean = sum(nhfgl(1:mGL))/real(NGL,SP)
    if(isnan(mean)) call fatal_error("nhfgl is nan")
    write(ipt,*) "NHFGL: min/max/mean",minval(nhfgl(1:MGL)),maxval(nhfgl(1:MGL)),mean

    WRITE(hfxunit) hour
    write(hfxunit)(NHFGL(i),SWRGL(i),i=1,mgl)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_HFX "

  END SUBROUTINE DUMP_HFX


  SUBROUTINE DUMP_EVP(NOW)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP) :: hour, mean
    INTEGER :: STATUS, I

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_EVP "

    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,PRC,prcgl)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,EVP,evpgl)
# endif
    END IF

    IF (.NOT. MSR) RETURN

    hour = DAYS(NOW - ZEROTIME) * 24.0_SP

    write(ipt,*) "time", hour

    mean = sum(evpgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Evpgl is nan")
    write(ipt,*) "EVPGL: min/max/mean",minval(Evpgl(1:MGL)),maxval(Evpgl(1:MGL)),mean

    mean = sum(Prcgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Prcgl is nan")
    write(ipt,*) "PRCGL: min/max/mean",minval(Prcgl(1:MGL)),maxval(Prcgl(1:MGL)),mean

    WRITE(evpunit) hour
    write(evpunit)(evpgl(i),prcgl(i),i=1,mgl)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_EVP "

  END SUBROUTINE DUMP_EVP


  SUBROUTINE DUMP_WND(NOW)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP) :: hour, mean
    INTEGER :: STATUS, I

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_WND "

    IF (PAR)THEN
# if defined(MULTIPROCESSOR)

       CALL PCOLLECT(MYID,MSRID,NPROCS,EMAP,STRX,WNDXGL)
       CALL PCOLLECT(MYID,MSRID,NPROCS,EMAP,STRY,WNDYGL)

# endif
    END IF

    IF (.NOT. MSR) RETURN

    hour = DAYS(NOW - ZEROTIME) * 24.0_SP

    write(ipt,*) "time", hour

    mean = sum(wndxgl(1:NGL))/real(NGL,SP)
    if(isnan(mean)) call fatal_error("wndxgl is nan")
    write(ipt,*) "WNDXGL: min/max/mean",minval(wndxgl(1:NGL)),maxval(wndxgl(1:NGL)),mean


    mean = sum(wndygl(1:NGL))/real(NGL,SP)
    if(isnan(mean)) call fatal_error("wndygl is nan")
    write(ipt,*) "WNDYGL: min/max/mean",minval(wndygl(1:NGL)),maxval(wndygl(1:NGL)),mean

    WRITE(wndunit) hour
    write(wndunit)(WNDXGL(i),WNDYGL(i),i=1,ngl)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_WND "

  END SUBROUTINE DUMP_WND

  SUBROUTINE DUMP_AIP(NOW)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP) :: hour, mean
    INTEGER :: STATUS, I

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_AIP "

    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,AIP,aipgl)
# endif
    END IF

    IF (.NOT. MSR) RETURN

    hour = DAYS(NOW - ZEROTIME) * 24.0_SP

    write(ipt,*) "time", hour

    mean = sum(aipgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Aipgl is nan")
    write(ipt,*) "AIPGL: min/max/mean",minval(Aipgl(1:MGL)),maxval(Aipgl(1:MGL)),mean

    WRITE(aipunit) hour
    write(aipunit)(aipgl(i),i=1,mgl)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_AIP "

  END SUBROUTINE DUMP_AIP

  SUBROUTINE DUMP_ICE(NOW)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP) :: hour, mean
    INTEGER :: STATUS, I

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_ICE "

    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,SAT,satgl)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,SPQ,spqgl)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,CLD,cldgl)
# endif
    END IF

    IF (.NOT. MSR) RETURN

    hour = DAYS(NOW - ZEROTIME) * 24.0_SP

    write(ipt,*) "time", hour

    mean = sum(satgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Satgl is nan")
    write(ipt,*) "SATGL: min/max/mean",minval(Satgl(1:MGL)),maxval(Satgl(1:MGL)),mean

    mean = sum(spqgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Spqgl is nan")
    write(ipt,*) "SPQGL: min/max/mean",minval(Spqgl(1:MGL)),maxval(Spqgl(1:MGL)),mean

    mean = sum(cldgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Cldgl is nan")
    write(ipt,*) "CLDGL: min/max/mean",minval(Cldgl(1:MGL)),maxval(Cldgl(1:MGL)),mean

    WRITE(satunit) hour
    write(satunit)(satgl(i),i=1,mgl)

    WRITE(spqunit) hour
    write(spqunit)(spqgl(i),i=1,mgl)

    WRITE(cldunit) hour
    write(cldunit)(cldgl(i),i=1,mgl)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_ICE "

  END SUBROUTINE DUMP_ICE

! EJA <
  SUBROUTINE DUMP_SOLAR(NOW)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP) :: hour, mean
    INTEGER :: STATUS, I

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_SOLAR "

    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,SAT,satgl)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,DPT,dptgl)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,CLD,cldgl)
# endif
    END IF

    IF (.NOT. MSR) RETURN

    hour = DAYS(NOW - ZEROTIME) * 24.0_SP

    write(ipt,*) "time", hour

    mean = sum(satgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Satgl is nan")
    write(ipt,*) "SATGL: min/max/mean",minval(Satgl(1:MGL)),maxval(Satgl(1:MGL)),mean

    mean = sum(dptgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Dptgl is nan")
    write(ipt,*) "DPTGL: min/max/mean",minval(Dptgl(1:MGL)),maxval(Dptgl(1:MGL)),mean

    mean = sum(cldgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Cldgl is nan")
    write(ipt,*) "CLDGL: min/max/mean",minval(Cldgl(1:MGL)),maxval(Cldgl(1:MGL)),mean

    WRITE(satunit) hour
    write(satunit)(satgl(i),i=1,mgl)

    WRITE(dptunit) hour
    write(dptunit)(dptgl(i),i=1,mgl)

    WRITE(cldunit) hour
    write(cldunit)(cldgl(i),i=1,mgl)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_SOLAR "

  END SUBROUTINE DUMP_SOLAR
! EJA >
  
    
END MODULE binout


MODULE MOD_output 
  USE MOD_SURFACEFORCE
  USE MOD_MM5
  USE MOD_WRF
  USE MOD_BINARY

  TYPE(GRID), SAVE :: MYGRID

  CONTAINS

  SUBROUTINE DUMP_BINARY(NOW)
    USE BINOUT
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_BINARY"

    IF(PRECIPITATION_ON) THEN
       CALL DUMP_EVP(NOW)
    END IF

    IF(WIND_ON) THEN
       CALL DUMP_WND(NOW)
    END IF

    IF(HEATING_ON) THEN
       CALL DUMP_HFX(NOW)
    END IF

    IF(AIRPRESSURE_ON) THEN
       CALL DUMP_AIP(NOW)
    END IF

    IF(ICE_MODEL) THEN
       CALL DUMP_ICE(NOW)
    END IF
! EJA <
    IF(HEATING_SOLAR_ON) THEN
       CALL DUMP_SOLAR(NOW)
    END IF
! EJA >
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_BINARY"

  END SUBROUTINE DUMP_BINARY
  

!/////////////////////////////////////////////////////////////
  Subroutine MAKE_BINARY_OUTFILE(WNDNM,HFXNM,EVPNM,AIPNM,SATNM,SPQNM,CLDNM,DPTNM)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WNDNM,HFXNM,EVPNM,AIPNM,SATNM,SPQNM,CLDNM,DPTNM
    integer status
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:MAKE_BINARY_OUTFILE "

    IF (.NOT. MSR) RETURN

    IF(WIND_ON) THEN
       CALL FOPEN(WNDUNIT,WNDNM,'our')


       IF(PAR) THEN
          ALLOCATE(WNDXGL(0:NGL))
          ALLOCATE(WNDYGL(0:NGL))
       ELSE
          WNDXGL => STRX
          WNDYGL => STRY
       END IF
    END IF

    IF (HEATING_ON) THEN
       CALL FOPEN(HFXUNIT,HFXNM,'our')

       IF(PAR) THEN
          ALLOCATE(SWRGL(0:MGL))
          ALLOCATE(NHFGL(0:MGL))
       ELSE
          swrgl  => SWR
          nhfgl  => NHF
       END IF
    END IF

    IF(PRECIPITATION_ON) THEN
       CALL FOPEN(EVPUNIT,EVPNM,'our')

       IF(PAR) THEN
          ALLOCATE(evpgl(0:MGL))
          ALLOCATE(prcgl(0:MGL))
       ELSE
          evpgl => EVP
          prcgl => PRC
       END IF
    END IF

    IF(AIRPRESSURE_ON) THEN
       CALL FOPEN(AIPUNIT,AIPNM,'our')

       IF(PAR) THEN
          ALLOCATE(aipgl(0:MGL))
       ELSE
          aipgl => AIP
       END IF
    END IF

    IF(SAT_ON) THEN
       CALL FOPEN(SATUNIT,SATNM,'our')

       IF(PAR) THEN
          ALLOCATE(satgl(0:MGL))
       ELSE
          satgl => SAT
       END IF
    END IF

    IF(SPQ_ON) THEN
       CALL FOPEN(SPQUNIT,SPQNM,'our')

       IF(PAR) THEN
          ALLOCATE(spqgl(0:MGL))
       ELSE
          spqgl => SPQ
       END IF
    END IF

    IF(CLD_ON) THEN
       CALL FOPEN(CLDUNIT,CLDNM,'our')

       IF(PAR) THEN
          ALLOCATE(cldgl(0:MGL))
       ELSE
          cldgl => CLD
       END IF
    END IF
! EJA <
    IF(DPT_ON) THEN
       CALL FOPEN(DPTUNIT,DPTNM,'our')

       IF(PAR) THEN
          ALLOCATE(DPTgl(0:MGL))
       ELSE
          dptgl => DPT
       END IF
    END IF
! EJA >


    ! CHECK TO MAKE SURE ATLEAST ONE IS ACTIVE
    IF (.not. PRECIPITATION_ON .and. .not. HEATING_ON .and. &
      & .not. WIND_ON .and. .not. AIRPRESSURE_ON .and. .not. ICE_MODEL .and. .not. HEATING_SOLAR_ON)  &
         CALL FATAL_ERROR("FOUND NO BINARY FORCING INPUT FILES?")

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:MAKE_BINARY_OUTFILE "

  END SUBROUTINE MAKE_BINARY_OUTFILE



!
! OUTPUT TYPES
!
!/////////////////////////////////////////////////////////////
  SUBROUTINE MAKE_UNSTRUCTURED_OUTFILE(FNAME)
    USE MOD_NCDIO
    USE MOD_FORCE, only : fvcom_cap_grid_source
    IMPLICIT NONE
    CHARACTER(LEN=*) :: FNAME
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    LOGICAL :: FOUND

    Call SET_FVCOM_GRID(MYGRID)    

    CALL DEFINE_DIMENSIONS(MYGRID)

    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    NC_OUT => NCF

    NCF%FTIME => new_ftime()

    NCF%FNAME = TRIM(FNAME)

    NCF => ADD(NCF,GRID_FILE_OBJECT(MYGRID) )
    
    NCF => ADD(NCF,TIME_FILE_OBJECT() )

    ATT => FIND_ATT(NCF,'source',FOUND)
    IF(.NOT.FOUND) CALL FATAL_ERROR("LOOKING FOR 'source' ATTRIBUTE: NOT FOUND")
    
    ATT%CHR = fvcom_cap_grid_source

    NCF => ADD(NCF,FILE_DATE_OBJECT() )

    IF (HEATING_ON) THEN
       ! SWRAD
       VAR  => NC_MAKE_PVAR(name='short_wave',&
            & values=SWR, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Short Wave Radiation') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='W m-2') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

       ! WTSURF - NET HEAT FLUX
       VAR  => NC_MAKE_PVAR(name='net_heat_flux',&
            & values=NHF, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Surface Net Heat Flux') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='W m-2') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

!       IF(ASSOCIATED(LNGWV)) THEN
          ! WTSURF - NET HEAT FLUX
!          VAR  => NC_MAKE_PVAR(name='lngwv',&
          VAR  => NC_MAKE_PVAR(name='long_wave',&
               & values=LWR, DIM1= DIM_node, DIM2= DIM_time)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='long wave radiation') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='W m-2') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='type',values='data') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)
!       END IF

       IF(ASSOCIATED(SENSE)) THEN
          ! WTSURF - NET HEAT FLUX
          VAR  => NC_MAKE_PVAR(name='sensible',&
               & values=sense, DIM1= DIM_node, DIM2= DIM_time)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='sensible heat flux') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='W m-2') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='type',values='data') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)
       END IF

       IF(ASSOCIATED(latent)) THEN
          ! WTSURF - NET HEAT FLUX
          VAR  => NC_MAKE_PVAR(name='latent',&
               & values=latent, DIM1= DIM_node, DIM2= DIM_time)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='latent heat flux') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='W m-2') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='type',values='data') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)
       END IF

    END IF

    IF(WIND_ON) THEN
       ! UUWIND
       VAR  => NC_MAKE_PVAR(name='uwind_stress',&
            & values=STRX, DIM1= DIM_nele, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Wind Stress') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Stress') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='Pa') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

       ! VVWIND
       VAR  => NC_MAKE_PVAR(name='vwind_stress',&
            & values=STRY, DIM1= DIM_nele, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Northward Wind Stress') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Stress') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='Pa') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)


       IF(ASSOCIATED(SPDX)) THEN
          ! WTSURF - NET HEAT FLUX
          VAR  => NC_MAKE_PVAR(name='U10',&
               & values=spdx, DIM1= DIM_nele, DIM2= DIM_time)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Wind Speed') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='m/s') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='type',values='data') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)
       END IF

       IF(ASSOCIATED(SPDY)) THEN
          ! WTSURF - NET HEAT FLUX
          VAR  => NC_MAKE_PVAR(name='V10',&
               & values=spdy, DIM1= DIM_nele, DIM2= DIM_time)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='Northward Wind Speed') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='m/s') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='type',values='data') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)
       END IF


    END IF

    IF(PRECIPITATION_ON) THEN
       VAR  => NC_MAKE_PVAR(name='precip',&
            & values=PRC, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Precipitation') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Precipitation, ocean &
            &lose water is negative') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

       ! EVAPORATION
       VAR  => NC_MAKE_PVAR(name='evap',&
            & values=EVP, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Evaporation') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Evaporation, ocean &
            &lose water is negative') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

       IF(ASSOCIATED(SST)) THEN
          ! WTSURF - NET HEAT FLUX
          VAR  => NC_MAKE_PVAR(name='sst',&
               & values=sst, DIM1= DIM_node, DIM2= DIM_time)
          
          ATT  => NC_MAKE_ATT(name='long_name',values='sea surface temperature') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='units',values='degrees_C"') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
          VAR  => ADD(VAR,ATT)
          
          ATT  => NC_MAKE_ATT(name='type',values='data') 
          VAR  => ADD(VAR,ATT)
          
          NCF  => ADD(NCF,VAR)
       END IF

!!$       IF(ASSOCIATED(SAT)) THEN
!!$          ! WTSURF - NET HEAT FLUX
!!$          VAR  => NC_MAKE_PVAR(name='SAT',&
!!$               & values=sat, DIM1= DIM_node, DIM2= DIM_time)
!!$          
!!$          ATT  => NC_MAKE_ATT(name='long_name',values='surface air temperature') 
!!$          VAR  => ADD(VAR,ATT)
!!$          
!!$          ATT  => NC_MAKE_ATT(name='units',values='degrees_C"') 
!!$          VAR  => ADD(VAR,ATT)
!!$          
!!$          ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
!!$          VAR  => ADD(VAR,ATT)
!!$          
!!$          ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
!!$          VAR  => ADD(VAR,ATT)
!!$          
!!$          ATT  => NC_MAKE_ATT(name='type',values='data') 
!!$          VAR  => ADD(VAR,ATT)
!!$          
!!$          NCF  => ADD(NCF,VAR)
!!$       END IF


    END IF

    IF(AIRPRESSURE_ON) THEN
       VAR  => NC_MAKE_PVAR(name='air_pressure',&
            & values=AIP, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Air Pressure') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Sea surface air pressure') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='Pa') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)
    END IF

    IF(ICE_MODEL) THEN
!       VAR  => NC_MAKE_PVAR(name='SAT',&
       VAR  => NC_MAKE_PVAR(name='air_temperature',&
            & values=SAT, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Air Temperature') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Sea surface air temperature') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='degree (C)') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)


       VAR  => NC_MAKE_PVAR(name='SPQ',&
            & values=SPQ, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Specific Humidity') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Sea surface specific humidity') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='kg/kg') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

       VAR  => NC_MAKE_PVAR(name='relative_humidity',&
            & values=RH, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Relative Humidity') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Sea surface relative humidity') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='kg/kg') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)


       VAR  => NC_MAKE_PVAR(name='cloud_cover',&
            & values=CLD, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Cloud Cover') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Cloud cover, 0 - 1') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values=' ') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

    END IF

! EJA <
    IF(HEATING_SOLAR_ON) THEN
!       VAR  => NC_MAKE_PVAR(name='SAT',&
       VAR  => NC_MAKE_PVAR(name='air_temperature',&
            & values=SAT, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Air Temperature') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Sea surface air temperature') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='degree (C)') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)


       VAR  => NC_MAKE_PVAR(name='dew_point',&
            & values=DPT, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Dew Point') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Dew Point Temperature') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='degree (C)') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)


       VAR  => NC_MAKE_PVAR(name='cloud_cover',&
            & values=CLD, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Cloud Cover') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='description',values='Cloud cover, 0 - 1') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values=' ') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

    END IF

! EJA >


    IF (WIND_TYPE=='speed' .and. WIND_ON) THEN

       ! UUWIND
       VAR  => NC_MAKE_PVAR(name='uwind_speed',&
            & values=SPDX, DIM1= DIM_nele, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Wind Speed') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Speed') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='m/s') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

       ! VVWIND
       VAR  => NC_MAKE_PVAR(name='vwind_speed',&
            & values=SPDY, DIM1= DIM_nele, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Northward Wind Speed') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Speed') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='m/s') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)

       NCF  => ADD(NCF,VAR)

    END IF

!!$    ! FOR TESTING THE BILINEAR METHOD UNDER A MASKED WRF CASE
!!$    VAR  => NC_MAKE_PVAR(name='NFOUND',&
!!$            & values=WRF_intp_n%N_FOUND, DIM1= DIM_node)
!!$    NCF  => ADD(NCF,VAR)
!!$
!!$    VAR  => NC_MAKE_PVAR(name='CFOUND',&
!!$            & values=WRF_intp_c%N_FOUND, DIM1= DIM_nele)
!!$    NCF  => ADD(NCF,VAR)
!!$
!!$    VAR  => NC_MAKE_PVAR(name='CCNT',&
!!$            & values=WRF_intp_c%N_CNT, DIM1= DIM_nele)
!!$    NCF  => ADD(NCF,VAR)
!!$
!!$    VAR  => NC_MAKE_PVAR(name='NCNT',&
!!$            & values=WRF_intp_n%N_CNT, DIM1= DIM_node)
!!$    NCF  => ADD(NCF,VAR)



    NCF%FTIME%NEXT_STKCNT = 0
    CALL NC_WRITE_FILE(NC_OUT)


  END SUBROUTINE MAKE_UNSTRUCTURED_OUTFILE
!/////////////////////////////////////////////////////////////
  SUBROUTINE MAKE_STRUCTURED_OUTFILE(FNAME)
    USE MOD_NCDIO
    IMPLICIT NONE
    CHARACTER(LEN=*) :: FNAME
    TYPE(NCFILE), POINTER :: NCF
    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()
    
    NC_OUT => NCF
    
    NCF%FTIME => new_ftime()
    
    NCF%FNAME = TRIM(FNAME)
    
    ! STRUCTURED OUPUT
    
    
    NCF%FTIME%NEXT_STKCNT = 0
    CALL NC_WRITE_FILE(NC_OUT)

  END SUBROUTINE MAKE_STRUCTURED_OUTFILE


  END MODULE MOD_output
