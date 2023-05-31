################################################################################
##!! @author Jianhua Qi (jqi@umassd.edu) based on Saeed Moghimi (moghimis@gmail.com)
##!! @date 07/08/19 Original documentation

##! LOG !
##! 
################################################################################

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)
#include ${PETSC_DIR}/bmake/common/variables
include ../src/make.inc
#include	../src/libfvcom.a
 FVCOMDIR := $(ROOTDIR)/FVCOM/src

#        DTLIBS          = -L/usr/local/install/julian-1.3.3/icc-14.0.2.144_ifort-14.0.2.144 -ljulian
#        DTINCS          = -I/usr/local/install/julian-1.3.3/icc-14.0.2.144_ifort-14.0.2.144/include
#        PARLIB          = -L/usr/local/install/metis/gcc-4.4.7/lib64 -lmetis
        DTLIBS          = -L/usr/local/install/libjulian/1.3.3/intel/19/lib64 -ljulian
        DTINCS          = -I/usr/local/install/libjulian/1.3.3/intel/19/include
        PARLIB          = -L/usr/local/install/metis/1.1/intel/19/lib -lmetis
	FVCOMLIBS	=	-L../src -lfvcom
	FVCOMINCS	=	-I../src/


# unremark this line to define petlists
LIBDIRS += -DWITHPETLISTS_on 

# debug for gfortran
# LIBDIRS +=  -g -O0 -fbacktrace  -lhdf5 -lhdf5_fortran
################################################################################ $(INCDIRS) $(LIBDIRS)
################################################################################

.SUFFIXES: .f90 .F90 .c .C

%.o : %.F90
#	$(ESMF_F90COMPILER) -c -check all -traceback $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(FVCOMINCS) $(DTINCS)  $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(FVCOMINCS) $(DTINCS)  $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<
#	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(FVCOMINCS) $(DTINCS)  $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

%.o : %.F
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(FVCOMINCS) $(DTINCS)  $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

# -----------------------------------------------------------------------------
#esmApp: esmApp.o esm.o fvcom_cap.o fvcom_mod.o wave_cap.o wave_mod.o nuopc_fvcom.o namelist_esmf.o mod_driver.o
#	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS)  -o $@ $^ $(ESMF_F90ESMFLINKLIBS) $(FVCOMLIBS) $(DTLIBS) $(PARLIB) $(INCDIRS) $(LIBDIRS) $(LIBSS) 
#	echo " CP EXE FILE ..... "

PWDDIR     := $(shell pwd)
INSTALLDIR := $(INSTDIR)    # coming from NEMS/src/incmake/component_WW3DATA.mk

MAKEFILE = fvcom.mk

LIBRARY  = libfvcom_cap.a
MODULES  = fvcom_cap.o fvcom_mod.o nuopc_fvcom.o namelist_esmf.o mod_driver.o

nuopc: $(LIBRARY)

$(LIBRARY): $(MODULES)
#	$(AR) $(ARFLAGS) $@ $?
	$(AR) $@ $^
#$(LIBRARY): $(MODULES)
#	$(AR) $@ $(MODULES)
	$(RANLIB) $@

nuopcinstall: $(LIBRARY)
	@echo ""
	mkdir -p $(INSTALLDIR)
	@echo "# ESMF self-describing build dependency makefile fragment" > $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo                                                            >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_FRONT     = fvcom_cap "                          >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_INCPATH   = $(PWDDIR) $(strip $(INSTALLDIR)) "   >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_CMPL_OBJS = "                                    >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_LINK_OBJS =  -L$(strip $(INSTALLDIR)) $(strip $(INSTALLDIR))/$(LIBRARY)  $(strip $(INSTALLDIR))/libfvcom.a \
					$(DTLIBS) $(PARLIB)" >> $(strip $(INSTALLDIR))/$(MAKEFILE)
	cp -f *.a *.mod $(INSTALLDIR) 
	cp -f  $(strip $(INSTALLDIR))/../FVCOM/src/libfvcom.a $(strip $(INSTALLDIR))
	@echo ""


# module dependencies:
esmApp.o:   esm.o
esm.o:      fvcom_cap.o wave_cap.o mod_driver.o
fvcom_cap.o: nuopc_fvcom.o
wave_cap.o: nuopc_fvcom.o
nuopc_fvcom.o: namelist_esmf.o
namelist_esmf.o: mod_driver.o
fvcom_cap.o:  fvcom_mod.o
wave_cap.o:   wave_mod.o

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: dust clean distclean info edit
dust:
	rm -f PET*.ESMF_LogFile esmApp fvcom_cap.o fvcom_cap.mod
clean:
	rm -f esmApp *.o *.mod *.a *.mk.install
distclean: dust clean

edit:
	gedit Makefile esmApp.F90 esm.F90 fvcom_cap.F90 fvcom_mod.F90 nuopc_fvcom.F90 namelist_esmf.F90 


#include ../FVCOM/make.in
