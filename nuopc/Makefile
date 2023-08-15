################################################################################
##!! @author Jianhua Qi (jqi@umassd.edu) based on Saeed Moghimi (moghimis@gmail.com)
##!! @date 07/08/19 Original documentation

##!! @date 05/25/23
##!! @author Panagiotis Velissariou (panagiotis.velissariou@noaa.gov)

##! LOG !
##! 
################################################################################


###########################################################################
# Mandatory environment variables NETCDFHOME, ESMFMKFILE
ifneq ($(origin NETCDFHOME), environment)
  $(error Environment variable NETCDFHOME was not set)
endif
#-----
ifneq ($(origin ESMFMKFILE), environment)
  $(error Environment variable ESMFMKFILE was not set)
endif
#-----

include $(ESMFMKFILE)

#include ${PETSC_DIR}/bmake/common/variables

# FVCOM related configuration
FVCOMDIR := $(COMP_SRCDIR)/src
FVCOMINCS = -I. -I$(FVCOMDIR)
FVCOMLIBS = -L$(FVCOMDIR) -lfvcom

include $(FVCOMDIR)/make.inc

# unremark this line to define petlists
LIBDIRS = -DWITHPETLISTS_on $(LIBS)
INCDIRS = $(INCS) $(FVCOMINCS)

# LIBDIRS +=  -g -O0 -fbacktrace  -lhdf5 -lhdf5_fortran
################################################################################ $(INCDIRS) $(LIBDIRS)
################################################################################

.SUFFIXES: .f90 .F90 .c .C

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

%.o : %.F
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

.PRECIOUS: %.o

# -----------------------------------------------------------------------------
vars:
	@echo " >>>>>>>>>>>>>>>> vars >>>>>>>>>>>>>>"
	@echo "INCDIRS              = $(INCDIRS)"
	@echo "LIBDIRS              = $(LIBDIRS)"
	@echo "ESMF_F90COMPILER     = $(ESMF_F90COMPILER)"
	@echo "ESMF_F90LINKER       = $(ESMF_F90LINKER)"
	@echo "ESMF_F90LINKOPTS     = $(ESMF_F90LINKOPTS)"
	@echo "ESMF_F90LINKPATHS    = $(ESMF_F90LINKPATHS)"
	@echo "ESMF_F90LINKRPATHS   = $(ESMF_F90LINKRPATHS)"
	@echo "ESMF_F90ESMFLINKLIBS = $(ESMF_F90ESMFLINKLIBS)"	
# -----------------------------------------------------------------------------

PWDDIR     := $(shell pwd)
INSTALLDIR := $(INSTDIR)    # coming from NEMS/src/incmake/component_FVCOM.mk

MAKEFILE = fvcom.mk

LIBRARY  = libfvcom_cap.a
MODULES  = fvcom_cap.o fvcom_mod.o nuopc_fvcom.o namelist_esmf.o mod_driver.o

nuopc: $(LIBRARY)

$(LIBRARY): $(MODULES)
	$(AR) $@ $? #$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

nuopcinstall: $(LIBRARY)
	@echo ""
	rm -f $(MAKEFILE)
	mkdir -p $(strip $(INSTALLDIR))
	@echo "# ESMF self-describing build dependency makefile fragment" > $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo                                                            >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_FRONT     = fvcom_cap "                          >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_INCPATH   = $(PWDDIR) $(strip $(INSTALLDIR)) "   >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_CMPL_OBJS = "                                    >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_LINK_OBJS = $(strip $(INSTALLDIR))/$(LIBRARY) -L$(strip $(INSTALLDIR)) -lfvcom \
				    $(LIBDIRS)"                          >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	cp -f *.a *.mod $(strip $(INSTALLDIR))
	cp -f  $(strip $(FVCOMDIR))/libfvcom.a $(strip $(INSTALLDIR))
	@echo ""


# module dependencies:
esmApp.o:        esm.o
esm.o:           fvcom_cap.o wave_cap.o mod_driver.o
fvcom_cap.o:     nuopc_fvcom.o
wave_cap.o:      nuopc_fvcom.o
nuopc_fvcom.o:   namelist_esmf.o
namelist_esmf.o: mod_driver.o
fvcom_cap.o:     fvcom_mod.o
wave_cap.o:      wave_mod.o

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: nuopc nuopcclean nuopcinstall vars dust nuopcdistclean

nuopcclean:
	rm -f  *.o *.mod $(MAKEFILE)
	make -C $(FVCOMDIR) allclean

nuopcdistclean: nuopcclean dust
	@echo ""
	@echo "Cleaning Model build..."
	@echo ""
	make -C $(FVCOMDIR) allclobber

dust:
	@echo ""
	@echo "Dusting ..."
	@echo
	rm -f PET*.ESMF_LogFile *.nc

edit:
	gedit Makefile esmApp.F90 esm.F90 fvcom_cap.F90 fvcom_mod.F90 nuopc_fvcom.F90 namelist_esmf.F90 

