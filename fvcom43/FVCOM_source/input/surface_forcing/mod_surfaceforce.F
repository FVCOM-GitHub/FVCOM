MODULE MOD_SURFACEFORCE
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_INPUT
  USE MOD_TIME
  USE MOD_SOLAR
  USE CONTROL
  USE LIMS
  USE ALL_VARS
  IMPLICIT NONE


  CHARACTER(len=80) :: RUN_FILE

  CHARACTER(len=80) :: Input_file
  CHARACTER(len=80) :: Input_file_type
  CHARACTER(len=80) :: output_file
  CHARACTER(len=80) :: OUTPUT_FILE_TYPE
  CHARACTER(len=80) :: OUT_INTERVAL
  CHARACTER(len=80) :: ZERO_TIME
!  TYPE(TIME)        :: ZEROTIME


  INTEGER :: input_type, output_type
  INTEGER, PARAMETER :: wrf=1
  INTEGER, PARAMETER :: mm5=2
  INTEGER, PARAMETER :: binary=3
  INTEGER, PARAMETER :: unstruct=4
  INTEGER, PARAMETER :: struct=5
  INTEGER, PARAMETER :: prototypical_hurricane_model=6
  


  CHARACTER(LEN=160) :: WRF_FILE_NAME
  CHARACTER(LEN=160) :: MM5HFX_FILE_NAME
  CHARACTER(LEN=160) :: MM5WND_FILE_NAME
  CHARACTER(LEN=160) :: MM5EVP_FILE_NAME

  CHARACTER(LEN=160) :: BWND_FILE_NAME
  CHARACTER(LEN=160) :: BHFX_FILE_NAME
  CHARACTER(LEN=160) :: BEVP_FILE_NAME

  CHARACTER(LEN=160) :: BAIP_FILE_NAME
  CHARACTER(LEN=160) :: BSAT_FILE_NAME
  CHARACTER(LEN=160) :: BSPQ_FILE_NAME
  CHARACTER(LEN=160) :: BCLD_FILE_NAME
  CHARACTER(LEN=160) :: BDPT_FILE_NAME

  CHARACTER(LEN=160) :: STRUCT_FILE_NAME
  CHARACTER(LEN=160) :: UNSTRUCT_FILE_NAME


  LOGICAL :: FEXIST

  TYPE(TIME), SAVE :: INTERVAL

  TYPE(NCFILE), POINTER :: NC_OUT


  ! FOR FVCOM GRID OUTPUT
  REAL(SP), POINTER :: SPDY(:), SPDX(:)
  REAL(SP), POINTER :: STRY(:), STRX(:)
  REAL(SP), POINTER :: SWR(:),  NHF(:),  LWR(:)
  REAL(SP), POINTER :: PRC(:),  EVP(:)

  REAL(SP), POINTER :: PRS(:),  SST(:),LNGWV(:),SENSE(:),LATENT(:)

  ! for ice model code & Heating-SOLAR
  REAL(SP), POINTER :: AIP(:),  SAT(:), SPQ(:), CLD(:), RH(:), DPT(:)

  ! FOR REGULAR GRID OUTPUT
  REAL(SP), POINTER :: SPDY_rg(:,:), SPDX_rg(:,:)
  REAL(SP), POINTER :: STRY_rg(:,:), STRX_rg(:,:)
  REAL(SP), POINTER :: SWR_rg(:,:),  NHF_rg(:,:)
  REAL(SP), POINTER :: PRC_rg(:,:),  EVP_rg(:,:)

  REAL(SP), POINTER :: PRS_rg(:,:),  SST_rg(:,:),LNGWV_rg(:,:),SENSE_rg(:,:),LATENT_rg(:,:)

  REAL(SP), POINTER :: AIP_rg(:,:),  SAT_rg(:,:), SPQ_rg(:,:), CLD_rg(:,:), DPT_rg(:,:)


  INTEGER :: DIMX,DIMY  ! X/Y dimension for regular grid data


  NAMELIST /NML_SURFACE_FORCE/           &
       & INPUT_DIR,                      &
       & OUTPUT_DIR,                     &
       & INPUT_FILE,                     &
       & INPUT_FILE_TYPE,                &
       & OUTPUT_FILE,                    &
       & OUTPUT_FILE_TYPE,               &
       & GRID_FILE,                      &
       & GRID_FILE_UNITS,                &
       & PROJECTION_REFERENCE,           &
       & DATE_FORMAT,                    &
       & START_DATE,                     &
       & TIMEZONE,                       &
       & END_DATE,                       &
       & OUT_INTERVAL,                   &
       & ZERO_TIME,                      &
       & DATE_REFERENCE


  ! FILE UNITS FOR OUTPUT
  INTEGER, PARAMETER :: wndunit=118 
  INTEGER, PARAMETER :: hfxunit=119
  INTEGER, PARAMETER :: evpunit=120

  INTEGER, PARAMETER :: aipunit=121 
  INTEGER, PARAMETER :: satunit=122
  INTEGER, PARAMETER :: spqunit=123
  INTEGER, PARAMETER :: cldunit=124
  INTEGER, PARAMETER :: dptunit=126

  INTEGER, PARAMETER :: mm5unit=125

CONTAINS


  SUBROUTINE GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)
    use mod_sng


    character(len=*), INTENT(IN)::CVS_Id  ! [sng] CVS Identification
    character(len=*), INTENT(IN)::CVS_Date ! [sng] Date string
    character(len=*), INTENT(IN)::CVS_Name ! [sng] File name string
    character(len=*), INTENT(IN)::CVS_Revision ! [sng] File revision string

    character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
    ! Command-line parsing
    character(80)::arg_val ! [sng] command-line argument value
    character(200)::cmd_ln ! [sng] command-line
    character(80)::opt_sng ! [sng] Option string
    character(2)::dsh_key ! [sng] command-line dash and switch
    character(200)::prg_ID ! [sng] Program ID

    integer::arg_idx ! [idx] Counting index
    integer::arg_nbr ! [nbr] Number of command-line arguments
    integer::opt_lng ! [nbr] Length of option

    ! Main code
    call ftn_strini(cmd_ln) ! [sng] sng(1:len)=NUL

    call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
    call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID

    arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments

    if (arg_nbr .LE. 0 ) then
       if(MSR) WRITE(IPT,*) "You must specify: '--filename=<namelist>' "
       if(MSR) Call MYHelpTxt
       call PSHUTDOWN
    end if

    arg_idx=1 ! [idx] Counting index
    do while (arg_idx <= arg_nbr)
       call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg, increment arg_idx
       dsh_key=arg_val(1:2) ! [sng] First two characters of option
       if (dsh_key == "--") then
          opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
          if (opt_lng <= 0) then
             if(MSR) write(IPT,*) "Long option has no name"
             call PSHUTDOWN
          end if

          opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
          if (dbg_lvl >= dbg_io) then
             if(MSR) write (6,"(5a,i3)") prg_nm(1:ftn_strlen(prg_nm)), &
                  ": DEBUG Double hyphen indicates multi-character option: ", &
                  "opt_sng = ",opt_sng(1:ftn_strlen(opt_sng)),", opt_lng = ",opt_lng
          end if
          if (opt_sng == "dbg" .or. opt_sng == "dbg_lvl" ) then
             call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level

          else if (opt_sng == "dbg_par" .or.opt_sng == "Dbg_Par"&
               & .or.opt_sng == "DBG_PAR") then

             dbg_par = .true.

          else if (opt_sng == "FileName" .or.opt_sng == "FILENAME"&
               & .or.opt_sng == "filename") then

             call ftn_arg_get(arg_idx,arg_val,RUN_FILE) ! [sng] Input file
             RUN_FILE=RUN_FILE(1:ftn_strlen(RUN_FILE))
             ! Convert back to a fortran string!

             !          else if (opt_sng == "Create_NameList" .or.opt_sng == "create_namelist"&
             !               & .or.opt_sng == "CREATE_NAMELIST") then
             !
             !             call ftn_arg_get(arg_idx,arg_val,NAMELIST_NAME)
             !             NAMELIST_NAME = NAMELIST_NAME(1:ftn_strlen(NAMELIST_NAME))
             !
             !             BLANK_NAMELIST = .true.


          else if (opt_sng == "help" .or.opt_sng == "HELP" .or. opt_sng&
               & == "Help") then

             if(MSR) call MYHelpTxt


             call PSHUTDOWN
!!$   THIS DOES NOT SEEM PRACTICAL - MODIFY THE RUN FILE INSTEAD
!!$          else if (opt_sng == "CrashRestart") then
!!$             call ftn_arg_get(arg_idx,arg_val,CrashRestart) ! [lgc] Logical

          else ! Option not recognized
             arg_idx=arg_idx-1 ! [idx] Counting index
             if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
          endif ! endif option is recognized
          ! Jump to top of while loop
          cycle ! C, F77, and F90 use "continue", "goto", and "cycle"
       endif ! endif long option

       if (dsh_key == "-V" .or.dsh_key == "-v" ) then

          if(MSR) write(IPT,*) prg_id
          call PSHUTDOWN

       else if (dsh_key == "-H" .or.dsh_key == "-h" ) then

          if(MSR) call MYHelpTxt
          Call PSHUTDOWN

       else ! Option not recognized
          arg_idx=arg_idx-1 ! [idx] Counting index
          if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
       endif ! endif arg_val


    end do ! end while (arg_idx <= arg_nbr)

    CALL dbg_init(IPT_BASE,.false.)

  END SUBROUTINE GET_COMMANDLINE

  SUBROUTINE MYHELPTXT
    IMPLICIT NONE

    WRITE(IPT,*) "Add better help here!"
    WRITE(IPT,*) "! OPTIONS:"
    WRITE(IPT,*) "! --filename=<a namelist file>"
    WRITE(IPT,*) "!   The namelist runfile for the program! "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   Namelist OPTIONS: "
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!   INPUT_DIR              # input path"
    WRITE(IPT,*) "!   OUTPUT_DIR             # output path"
    WRITE(IPT,*) "!   INPUT_FILE             # filename"
    WRITE(IPT,*) "!   INPUT_FILE_TYPE        # 'mm5', 'wrf', 'binary', 'prototypical hurricane model'"
    WRITE(IPT,*) "!   OUT_FILE_TYPE          # 'binary', 'unstructured', 'structured'" 
    WRITE(IPT,*) "!   GRID_FILE              # filename"
    WRITE(IPT,*) "!   GRID_FILE_UNITS        # 'meters' or 'degrees'"
    WRITE(IPT,*) "!   PROJECTION_REFERENCE   # see proj for help" 
    WRITE(IPT,*) "!   DATE_FORMAT            # 'ymd' or 'dmy' ?"
    WRITE(IPT,*) "!   START_DATE             # Start Here"
    WRITE(IPT,*) "!   TIMEZONE               # The Time zone or none"
    WRITE(IPT,*) "!   END_DATE               # End Here"
    WRITE(IPT,*) "!   OUT_INTERVAL           # 'seconds = n.n' or 'days=n.n'"
    WRITE(IPT,*) "!   ZERO_TIME              # For binary output, when is zero?"

    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    Exmaple Namelist"
    write(UNIT=IPT,NML=NML_SURFACE_FORCE)


    WRITE(IPT,*) "! NOTES: This program is parallel!"
  END SUBROUTINE MYHELPTXT


  SUBROUTINE GET_RUN_FILE
    IMPLICIT NONE
    integer :: ios, i
    Character(Len=120):: FNAME
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: GET_RUNFILE;"


    FNAME = RUN_FILE

    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read IO Information
    READ(UNIT=NMLUNIT, NML=NML_SURFACE_FORCE,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) &
            & write(UNIT=IPT,NML=NML_SURFACE_FORCE)
       CALL FATAL_ERROR("Can Not Read NameList NML_SURFACE_FORCE from file: "//trim(FNAME))
    end if
    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_SURFACE_FORCE)


  END SUBROUTINE GET_RUN_FILE

  SUBROUTINE SET_PARAMETERS
    IMPLICIT NONE

    


  ! SELECT INPUT TYPE  :: TRANSLATE TO A INTEGER TYPE
  SELECT CASE(INPUT_FILE_TYPE)
  CASE("wrf","WRF")
     input_type=wrf
  CASE('mm5','MM5')
     input_type=mm5
  CASE('binary stress','BINARY STRESS')
     input_type=binary
     WIND_TYPE = "stress"
  CASE('binary','BINARY')
     input_type=binary
     WIND_TYPE = "none"
  CASE('binary speed','BINARY SPEED')
     input_type=binary
     WIND_TYPE="speed"
  CASE('prototypical_hurricane_model','PROTOTYPICAL_HURRICANE_MODEL')
     input_type=prototypical_hurricane_model
     WIND_TYPE="speed"
  CASE DEFAULT
     CALL FATAL_ERROR("ILLEGAL INPUT FILE TYPE:"//TRIM(INPUT_FILE_TYPE))
  END SELECT


  ! SET OUTPUT TYPE :: TRANSLATE TO A INTEGER TYPE
  SELECT CASE(OUTPUT_FILE_TYPE)
  CASE('unstructured','UNSTRUCTURED')
     if(len_Trim(OUTPUT_FILE) == 0 ) CALL FATAL_ERROR&
          & ("You must specify an ouput_file name in the name list")
     output_type=unstruct
  CASE('structured','STRUCTURED')     
     if(len_Trim(OUTPUT_FILE) == 0 ) CALL FATAL_ERROR&
          & ("You must specify an ouput_file name in the name list")
     output_type=struct
  CASE('binary','BINARY')
     output_type=binary
  CASE DEFAULT     
     CALL FATAL_ERROR("ILLEGAL OUTPUT FILE TYPE:"//TRIM(OUTPUT_FILE_TYPE))
  END SELECT
  

  ! SETUP EACH TYPE OF INPUT AND OUTPUT
  SELECT CASE(input_type)
  CASE(WRF)

     WRF_FILE_NAME=TRIM(INPUT_DIR)//TRIM(INPUT_FILE)

     wrf_out: SELECT CASE(output_type)
     CASE(struct)

        STRUCT_FILE_NAME = TRIM(output_dir)//trim(output_file)
        


        CALL FATAL_ERROR("NOT SET UP FOR WRF INPUT, STRUCTURED OUTPUT YET...")



     CASE(unstruct) 
        IF (GRID_FILE_UNITS == 'meters' .and. .not. USE_PROJ) CALL FATAL_ERROR &
             & ('Proj must be active to convert a WRF file to a cartesian grid!')

        UNSTRUCT_FILE_NAME = TRIM(output_dir)//trim(output_file)



        ALLOCATE(WUSURF(0:NT),WVSURF(0:NT))
        STRX => WUSURF
        STRY => WVSURF
        
        ALLOCATE(Qprec(0:MT),Qevap(0:MT))
        PRC => Qprec
        EVP => Qevap
        
        ALLOCATE(SWRAD_WATTS(0:MT), WTSURF_WATTS(0:MT))
        SWR => SWRAD_WATTS
        NHF => WTSURF_WATTS
        
# if defined (AIR_PRESSURE)
        ALLOCATE(PA_AIR(0:MT))
        AIP => PA_AIR
        AIRPRESSURE_ON = .TRUE.
# endif
        
# if defined(ICE)
        ALLOCATE(T_AIR(0:MT), QA_AIR(0:MT), CLOUD(0:MT))
        SPQ => QA_AIR
        CLD => CLOUD
        SAT => T_AIR
# endif

# if defined(HEATING_SOLAR)
        ALLOCATE(T_AIR(0:MT), CLOUD(0:MT), TDEW(0:MT))
        CLD => CLOUD
        SAT => T_AIR
        DPT => TDEW
# endif



        CALL RUN_TGE


     CASE(binary)
        IF (GRID_FILE_UNITS == 'meters' .and. .not. USE_PROJ) CALL FATAL_ERROR &
             & ('Proj must be active to convert a WRF file to a cartesian grid!')

        
        BWND_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_wnd.dat"
        BHFX_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_hfx.dat"
        BEVP_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_evp.dat"

        BAIP_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_aip.dat"
        BSAT_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_sat.dat"
        BSPQ_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_spq.dat"
        BCLD_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_cld.dat"
        BDPT_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_dpt.dat"



        ALLOCATE(WUSURF(0:NT),WVSURF(0:NT))
        STRX => WUSURF
        STRY => WVSURF

        ALLOCATE(Qprec(0:MT),Qevap(0:MT))
        PRC => Qprec
        EVP => Qevap

        ALLOCATE(SWRAD_WATTS(0:MT), WTSURF_WATTS(0:MT))
        SWR => SWRAD_WATTS
        NHF => WTSURF_WATTS

# if defined (AIR_PRESSURE)
        ALLOCATE(PA_AIR(0:MT))
        AIP => PA_AIR
        AIRPRESSURE_ON = .TRUE.
# endif

# if defined(ICE)
        ALLOCATE(T_AIR(0:MT), QA_AIR(0:MT), CLOUD(0:MT))
        SPQ => QA_AIR
        CLD => CLOUD
        SAT => T_AIR
# endif

# if defined(HEATING_SOLAR)
        ALLOCATE(T_AIR(0:MT), CLOUD(0:MT), TDEW(0:MT))
        CLD => CLOUD
        SAT => T_AIR
        DPT => TDEW
# endif

        CALL RUN_TGE



     CASE DEFAULT
        CALL FATAL_ERROR("THIS SHOULD NEVER HAPPEN?")
     END SELECT wrf_out

     

  CASE(mm5)

     MM5HFX_FILE_NAME = TRIM(INPUT_DIR)//"heat."//TRIM(INPUT_FILE)

     MM5EVP_FILE_NAME = TRIM(INPUT_DIR)//"pre_evap."//TRIM(INPUT_FILE)

     MM5WND_FILE_NAME = TRIM(INPUT_DIR)//"stress."//TRIM(INPUT_FILE)


     mm5_out: SELECT CASE(output_type)
     CASE(struct)

        STRUCT_FILE_NAME = TRIM(output_dir)//trim(output_file)

        ALLOCATE(STRX_RG(DIMX,DIMY),STRY_RG(DIMX,DIMY))
        ALLOCATE(SPDX_RG(DIMX,DIMY),SPDY_RG(DIMX,DIMY))
        ALLOCATE(PRC_rg(dimx,dimy),EVP_rg(dimx,dimy))
        ALLOCATE(SWR_rg(dimx,dimy), NHF_rg(dimx,dimy))
        
        ALLOCATE(PRS_rg(dimx,dimy), SAT_rg(dimx,dimy))
        ALLOCATE(SST_rg(dimx,dimy),LNGWV_rg(dimx,dimy))
        ALLOCATE(SENSE_rg(dimx,dimy),LATENT_rg(dimx,dimy))


        
     CASE(unstruct)
        IF (GRID_FILE_UNITS == 'meters' .and. .not. USE_PROJ) CALL FATAL_ERROR &
             & ('Proj must be active to convert a MM5 file to a cartesian grid!')

        UNSTRUCT_FILE_NAME = TRIM(output_dir)//trim(output_file)

        
        ALLOCATE(STRX(0:NT),STRY(0:NT))
        ALLOCATE(SPDX(0:NT),SPDY(0:NT))
        ALLOCATE(PRC(0:MT),EVP(0:MT))
        ALLOCATE(SWR(0:MT), NHF(0:MT))
        
        ALLOCATE(PRS(0:MT), SAT(0:MT),SST(0:MT),LNGWV(0:MT),SENSE(0:MT),LATENT(0:MT))
     
        AIP => PRS

        AIRPRESSURE_ON = .TRUE.

        CALL RUN_TGE

     CASE(binary)
        IF (GRID_FILE_UNITS == 'meters' .and. .not. USE_PROJ) CALL FATAL_ERROR &
             & ('Proj must be active to convert a MM5 file to a cartesian grid!')


       
        BWND_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_wnd.dat"
        BHFX_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_hfx.dat"
        BEVP_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_evp.dat"

        BAIP_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_aip.dat"
        BSAT_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_sat.dat"
        BSPQ_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_spq.dat"
        BCLD_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_cld.dat"
        BDPT_FILE_NAME= TRIM(output_dir)//trim(output_file)//"_dpt.dat"




        ALLOCATE(STRX(0:NT),STRY(0:NT))
        ALLOCATE(SPDX(0:NT),SPDY(0:NT))
        ALLOCATE(PRC(0:MT),EVP(0:MT))
        ALLOCATE(SWR(0:MT), NHF(0:MT))
        
        ALLOCATE(PRS(0:MT), SAT(0:MT),SST(0:MT),LNGWV(0:MT),SENSE(0:MT),LATENT(0:MT))


        CALL RUN_TGE


     CASE DEFAULT
        CALL FATAL_ERROR("THIS SHOULD NEVER HAPPEN?")
     END SELECT mm5_out


  CASE(binary)

     
     BWND_FILE_NAME= TRIM(input_dir)//trim(input_file)//"_wnd.dat"
     BHFX_FILE_NAME= TRIM(input_dir)//trim(input_file)//"_hfx.dat"
     BEVP_FILE_NAME= TRIM(input_dir)//trim(input_file)//"_evp.dat"

     BAIP_FILE_NAME= TRIM(input_dir)//trim(input_file)//"_aip.dat"
     BSAT_FILE_NAME= TRIM(input_dir)//trim(input_file)//"_sat.dat"
     BSPQ_FILE_NAME= TRIM(input_dir)//trim(input_file)//"_spq.dat"
     BCLD_FILE_NAME= TRIM(input_dir)//trim(input_file)//"_cld.dat"
     BDPT_FILE_NAME= TRIM(input_dir)//trim(input_file)//"_dpt.dat"
     

     ALLOCATE(STRX(0:NT),STRY(0:NT))
     ALLOCATE(PRC(0:MT),EVP(0:MT))
     ALLOCATE(SWR(0:MT), NHF(0:MT), LWR(0:MT))

     IF(WIND_TYPE == 'speed') THEN
        ALLOCATE(SPDX(0:NT),SPDY(0:NT))
     END IF


# if defined (AIR_PRESSURE)
     ALLOCATE(AIP(0:MT))
#endif

# if defined (ICE)
     ALLOCATE(SAT(0:MT), SPQ(0:MT), CLD(0:MT), RH(0:MT))
#endif

# if defined (HEATING_SOLAR)
     ALLOCATE(SAT(0:MT), DPT(0:MT), CLD(0:MT))
#endif

     binary_out: SELECT CASE(output_type)
     CASE(struct)
        
        CALL FATAL_ERROR &
             & ("Can not create a structred grid file from an unstructured binary file!")
     
     CASE(unstruct)


        UNSTRUCT_FILE_NAME = TRIM(output_dir)//trim(output_file)


     CASE(binary)
        
        CALL FATAL_ERROR("WHY WOULD YOU SELECT BOTH BINARY INPUT AND BINARY OUTPUT?")
        
     CASE DEFAULT
        CALL FATAL_ERROR("THIS SHOULD NEVER HAPPEN?")
     END SELECT binary_out

  CASE(prototypical_hurricane_model)

     
     BWND_FILE_NAME= TRIM(input_dir)//trim(input_file)    !//"_path.dat"

     ALLOCATE(STRX(0:NT),STRY(0:NT))
!     ALLOCATE(PRC(0:MT),EVP(0:MT))
!     ALLOCATE(SWR(0:MT), NHF(0:MT))

     IF(WIND_TYPE == 'speed') THEN
        ALLOCATE(SPDX(0:NT),SPDY(0:NT))
     END IF

# if defined (AIR_PRESSURE)
     ALLOCATE(AIP(0:MT))
     AIRPRESSURE_ON = .TRUE.
#endif

     p_t_model_out: SELECT CASE(output_type)
     CASE(struct)
        
        CALL FATAL_ERROR &
             & ("Can not create a structred grid file from an unstructured binary file!")
     
     CASE(unstruct)


        UNSTRUCT_FILE_NAME = TRIM(output_dir)//trim(output_file)


     CASE(binary)
        
        CALL FATAL_ERROR("WHY WOULD YOU SELECT BOTH BINARY INPUT AND BINARY OUTPUT?")
        
     CASE DEFAULT
        CALL FATAL_ERROR("THIS SHOULD NEVER HAPPEN?")
     END SELECT p_t_model_out


     
  CASE DEFAULT
     CALL FATAL_ERROR("THIS SHOULD NEVER HAPPEN?")
  END SELECT
   
    
END SUBROUTINE SET_PARAMETERS


  SUBROUTINE GET_FVCOM_GRID
    USE MOD_SETUP
    IMPLICIT NONE
    CHARACTER(LEN=80) FNAME
    INTEGER STATUS

    ! OPEN AND READ THE FVCOM GRID FILE
    IF (MSR) THEN
       FNAME = TRIM(INPUT_DIR)//TRIM(GRID_FILE)
       WRITE(IPT,*) "OPENING GRIDFILE: "//TRIM(FNAME)
       Call FOPEN(GRIDUNIT,TRIM(FNAME),'cfr')
    END IF

    CALL LOAD_COLDSTART_GRID(NVG)
    KB = 1

    CALL SETUP_DOMAIN

    IF(MSR) THEN
       ! ALLOCATE SPACE FOR THE GLOBAL GRID DATA
       ALLOCATE(Y_GBL(0:MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_GBL")
       ALLOCATE(X_GBL(0:MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_GBL")
    END IF

    ! ALLOCATE SPACE FOR THE LOCAL GRID DATA
    ALLOCATE(Y_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_LCL")
    ALLOCATE(X_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_LCL")


    ALLOCATE(VX(0:MT),VY(0:MT),XM(0:MT),YM(0:MT),LON(0:MT),LAT(0:MT))
    ALLOCATE(XC(0:NT),XMC(0:NT),YC(0:NT),YMC(0:NT),LONC(0:NT),LATC(0:NT))

    CALL LOAD_COLDSTART_COORDS(X_GBL,Y_GBL,X_LCL,Y_LCL)
    CALL COORDINATE_UNITS(X_LCL,Y_LCL)
    CALL SETUP_CENTER_COORDS

    DEALLOCATE(X_LCL)
    DEALLOCATE(Y_LCL)

    IF(MSR) THEN
       DEALLOCATE(X_GBL)
       DEALLOCATE(Y_GBL)
    END IF


  END SUBROUTINE GET_FVCOM_GRID

  SUBROUTINE RUN_TGE
    USE MOD_OBCS, only : iobcn
    implicit none
    INTEGER :: NCT

    NCT = NT*3
    IOBCN = 0

    ALLOCATE(NBE(0:NT,3))         ;NBE      = 0 !!INDICES OF ELMNT NEIGHBORS
    ALLOCATE(NTVE(0:MT))          ;NTVE     = 0 
    ALLOCATE(NTSN(MT))            ;NTSN     = 0 
    ALLOCATE(ISONB(0:MT))         ;ISONB    = 0  !!NODE MARKER = 0,1,2
    ALLOCATE(ISBCE(0:NT))         ;ISBCE    = 0 
    ALLOCATE(NIEC(NCT,2))         ;NIEC     = 0
    ALLOCATE(NTRG(NCT))           ;NTRG     = 0    
    ! POSITION OF NODAL CONTROL VOLUME CORNERS 
    ALLOCATE(XIJE(NCT,2))         ;XIJE     = ZERO
    ALLOCATE(YIJE(NCT,2))         ;YIJE     = ZERO 

    ! LENGTH OF NODAL CONTROL VOLUME EDGES
    ALLOCATE(DLTXE(NCT))          ;DLTXE    = ZERO
    ALLOCATE(DLTYE(NCT))          ;DLTYE    = ZERO
    ALLOCATE(DLTXYE(NCT))         ;DLTXYE   = ZERO !! TOTAL LENGTH
    ALLOCATE(SITAE(NCT))          ;SITAE    = ZERO !! ANGLE

    CALL TRIANGLE_GRID_EDGE


  END SUBROUTINE RUN_TGE



  SUBROUTINE SET_TIME
    USE MOD_SET_TIME
    IMPLICIT NONE
    integer status
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG

    IF (USE_REAL_WORLD_TIME) THEN
       StartTime = read_datetime(START_DATE,DATE_FORMAT,TIMEZONE,status)
       if (status /= 1) & 
            & Call Fatal_Error("Could not read the date string START_DATE: ", trim(START_DATE))

       ! GET THE END TIME
       EndTime = READ_DateTime(END_DATE,DATE_FORMAT,TIMEZONE,status)
       if (status /= 1) &
            & Call Fatal_Error("Could not read the date string END_DATE:", trim(END_DATE))

       IF(input_type /= wrf) THEN
          ZeroTime = read_datetime(ZERO_TIME,DATE_FORMAT,TIMEZONE,status)
          IF (STATUS  /= 1) CALL FATAL_ERROR("COULD NOT READ ZEROTIME STRING")
       END IF

    ELSE

       CALL IDEAL_TIME_STRING2TIME(START_DATE,FLAG,StartTime,dummy)
       IF(FLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")

       CALL IDEAL_TIME_STRING2TIME(END_DATE,FLAG,EndTime,dummy)
       IF(FLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")

        IF(input_type /= wrf) THEN
           CALL IDEAL_TIME_STRING2TIME(ZERO_TIME,FLAG,ZeroTime,dummy)
           IF(FLAG == 'step') CALL FATAL_ERROR&
                &("You must specify a time, not a step, for this restart file", &
                & "The Step will be set by the old restart file...")
        END IF

    END IF

    CALL PRINT_REAL_TIME(STARTTIME,IPT,"READ: START_DATE",TIMEZONE)
    CALL PRINT_REAL_TIME(ENDTIME,IPT,"READ: END_DATE",TIMEZONE)
    CALL PRINT_REAL_TIME(ZEROTIME,IPT,"READ: ZERO_TIME",TIMEZONE)


    ! SANITY ECK
    if(StartTime .GT. EndTime) &
         & Call Fatal_Error("Runfile Start_Date exceeds or equal to End_Date")

    CALL IDEAL_TIME_STRING2TIME(OUT_INTERVAL,FLAG,INTERVAL,dummy)
    ! SANITY CHECK
    IF (FLAG /= 'time') CALL FATAL_ERROR&
         & ("TIME INTERVAL SPECIFICATION IS INCORRENT",&
         & "INTERVAL MUST BE 'seconds=' or 'days='")

  END SUBROUTINE SET_TIME

END MODULE MOD_SURFACEFORCE
