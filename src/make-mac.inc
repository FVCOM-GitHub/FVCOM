#/===========================================================================/
# macOS Configuration Template for FVCOM
#
# This is a reference configuration for building FVCOM on macOS (Apple Silicon
# and Intel) using gfortran and OpenMPI from Homebrew.
#
# To use this template:
# 1. Copy this file: cp make-mac.inc make.inc
# 2. Edit the TOPDIR and INSTALLDIR paths below
# 3. Adjust library paths if your Homebrew installation differs
# 4. Run: make clean && make
#
# See docs/BUILD_MACOS.md for detailed build instructions
#/===========================================================================/

#========== TOPDIR ========================================================
# TOPDIR is the directory in which this make file and the fvcom source reside
# !!! CHANGE THIS TO YOUR ACTUAL PATH !!!

           TOPDIR        = /path/to/FVCOM/src

#========== INSTALLDIR =====================================================
# INSTALLDIR is where you installed the Julian library
# !!! CHANGE THIS TO YOUR ACTUAL PATH !!!

           INSTALLDIR    = /path/to/FVCOM/libs/install

#===========================================================================
# PREPROCESSOR OPTIONS FOR CPP
           DEF_FLAGS     = -P -traditional

#===========================================================================
# MEDM ENVIRONMENT - Use environmental variables if set
colon=:
empty=
dashI= $(empty) -I
dashL= $(empty) -L

ifeq ($(origin LIBPATH), environment)
  LIBDIR = -L$(subst $(colon),$(dashL),$(LIBPATH))
endif
ifeq ($(origin INCLUDEPATH), environment)
  INCDIR = -I$(subst $(colon),$(dashI),$(INCLUDEPATH))
endif

#===========================================================================
# LOCAL INSTALL - Point to your local Julian library installation
           LIBDIR       =  -L$(INSTALLDIR)/lib
           INCDIR       =  -I$(INSTALLDIR)/include

#--------------------------------------------------------------------------
#	STANDARD LIBRARIES FOR DATA AND TIME
	   DTLIBS	= -ljulian
	   DTINCS	=

#--------------------------------------------------------------------------
#        NETCDF - Homebrew paths for Apple Silicon
#        For Intel Mac, change /opt/homebrew to /usr/local
#--------------------------------------------------------------------------
           FLAG_USE_NETCDF4     = -DUSE_NETCDF4
#          FLAG_USE_COMPRESSION = -DUSE_COMPRESSION
           IOLIBS       =  -L/opt/homebrew/lib -lnetcdff -lnetcdf
           IOINCS       =  -I/opt/homebrew/include

#==========================================================================
#  MODEL OPTIONS - Typical configuration for development/testing
#==========================================================================

# PRECISION - Double precision with single precision output
           FLAG_1 =  -DDOUBLE_PRECISION -DSINGLE_OUTPUT

# SPHERICAL COORDINATES
           FLAG_2 = -DSPHERICAL

# WET/DRY TREATMENT
           FLAG_3 = -DWET_DRY

# MPI PARALLELIZATION with METIS 5
           FLAG_4 = -DMULTIPROCESSOR
           FLAG_411 = -DMETIS_5
           PARLIB =
           PARTINCS = -I/opt/homebrew/include
           PARTLIBS = -L/opt/homebrew/lib -lmetis

# WATER QUALITY MODEL (uncomment to enable)
#          FLAG_5 = -DWATER_QUALITY

# PROJECTION (uncomment to enable)
#          FLAG_6 = -DPROJ
#          PROJLIBS       = -lfproj4 -lproj -lm
#          PROJINCS       =

# DATA ASSIMILATION (uncomment to enable, requires PETSc)
#          FLAG_7 = -DDATA_ASSIM
#          include ${PETSC_DIR}/bmake/common/variables
#          OILIB  = -lmkl_lapack -lmkl_em64t -lguide -lpthread

# ADVECTION LIMITER - Choose ONE (REQUIRED)
           FLAG_8 =  -DLIMITED_NO
#          FLAG_8 =  -DLIMITED_1
#          FLAG_8 =  -DLIMITED_2

# SEMI-IMPLICIT (uncomment to enable, requires PETSc)
#          FLAG_9 = -DSEMI_IMPLICIT
#          include ${PETSC_DIR}/bmake/common/variables

# SOLID BOUNDARY TREATMENT - Choose ONE (REQUIRED)
           FLAG_10  = -DGCN
#          FLAG_10  = -DGCY1
#          FLAG_10 = -DGCY2

# GOTM TURBULENCE MODEL (uncomment to enable)
#          FLAG_11 = -DGOTM
#          GOTMLIB       = -L../GOTM_source -lturbulence -lutil
#          GOTMINCS      = -I../GOTM_source

# EQUILIBRIUM TIDE (uncomment to enable)
#          FLAG_12 = -DEQUI_TIDE

# ATMOSPHERIC TIDE (uncomment to enable)
#          FLAG_13 = -DATMO_TIDE

# RIVER DISTRIBUTION (uncomment for old-style floating point)
#          FLAG_14 = -DRIVER_FLOAT

# TVD (Total Variational Diminishing) - Recommended for better accuracy
           FLAG_44 = -DTVD

# MPDATA SCHEME
           FLAG_15 = -DMPDATA

# ADDITIONAL OPTIONS (uncomment as needed)
#          FLAG_16 = -DLAG_PARTICLE
#          FLAG_17 = -DDYE_RELEASE
#          FLAG_18 = -DTWO_D_MODEL
#          FLAG_19 = -DICE
#          FLAG_20 = -DICE_EMBEDDING
#          FLAG_21 = -DSEDIMENT
#          FLAG_22 = -DSED_CSTMS
#          FLAG_23 = -DOFFLINE_BIOLOGY
#          FLAG_24 = -DOFFLINE_SEDIMENT
#          FLAG_25 = -DVISIT
#          FLAG_26 = -DVIS_OMVIZ
#          FLAG_27 = -DBIO_CALC
#          FLAG_28 = -DPROBES
#          FLAG_29 = -DSTATION_TIMESERIES
#          FLAG_30 = -DBALANCE_2D
#          FLAG_31 = -DENKF
#          FLAG_32 = -DENKF_DETAILS
#          FLAG_33 = -DRRKASSIM
#          FLAG_34 = -DTIGHT_COUPLING
#          FLAG_35 = -DONE_D_TIDEBC
#          FLAG_36 = -DBALANCE_2D_VIA_FVCOM
#          FLAG_37 = -DBALANCE_3D_VIA_OCEAN
#          FLAG_38 = -DDTI_RESET_S
#          FLAG_39 = -DDTI_RESET_M
#          FLAG_40 = -DNESTING
#          FLAG_41 = -DMEAN_FLOW_RIVER
#          FLAG_42 = -DPLBC
#          FLAG_43 = -DDAM_BREAK
#          FLAG_45 = -DFLUID_MUD
#          FLAG_46 = -DOUT_STATION
#          FLAG_47 = -DOUT_SPARSE
#          FLAG_48 = -DDOUBLE_ECHO
#          FLAG_49 = -DBOUNDSCHK
#          FLAG_50 = -DESMF_NESTING
#          FLAG_51 = -DFLOCMOD
#          FLAG_52 = -DTPXO_TIDE
#          FLAG_53 = -DMLD_RHO
#          FLAG_54 = -DVEGETATION
#          FLAG_55 = -DSEDIMENT_HEATING

#--------------------------------------------------------------------------
#  gfortran/OpenMPI Compiler Configuration (Homebrew)
#--------------------------------------------------------------------------
         CPP      = /usr/bin/cpp
         COMPFLAG = -DGFORTRAN
         CC       = mpicc
         CXX      = mpicxx
         CFLAGS   = -O3
         FC       = mpif90
         DEBFLGS  =
         OPT      = -O3 -ffree-line-length-none
         CLIB     =

#==========================================================================
#  END OF USER CONFIGURATION
#==========================================================================
