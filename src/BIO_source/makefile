
# Master Makefile for making the executable of the generalized biological model.

#.SUFFIXES:
#.SUFFIXES: .F
#.SUFFIXES: .o .f90 .F .F90 
	SHELL	= /bin/sh
#	DEF_FLAGS     = -P -C -traditional 
#	NETCDFINCDIR	= /usr/local/include
# 	ETCDFLIB       = /usr/local/lib -lnetcdf

#         CPP      = /usr/bin/cpp 
#         CPPFLAGS = $(DEF_FLAGS) #-DINTEL   
#         FC    = ifort -checkall
#         FC       = mpif90 -static-libcxa -i-static

BIOMODULE	= \
libbiomodule.a(mod_1D.o)			\
libbiomodule.a(mod_utils2.o)			\
libbiomodule.a(mod_bioinp.o)			\
libbiomodule.a(mod_newinp2.o)			\
libbiomodule.a(mod_phytoplankton.o)		\
libbiomodule.a(mod_zooplankton.o)		\
libbiomodule.a(mod_bacteria.o)			\
libbiomodule.a(mod_detritus.o)			\
libbiomodule.a(mod_DOM.o)			\
libbiomodule.a(mod_nutrient.o)			\
libbiomodule.a(phytoplankton.o)			\
libbiomodule.a(zooplankton.o)			\
libbiomodule.a(bacteria.o)			\
libbiomodule.a(detritus.o)			\
libbiomodule.a(DOM.o)				\
libbiomodule.a(nutrient.o)			\
libbiomodule.a(get_parameter.o)			\
libbiomodule.a(bio_mixing.o)

LIBS	=	libbiomodule.a

FFLAGS_BIO = -O3

libbiomodule.a: $(BIOMODULE)

.SUFFIXES: .o .f90 .F .F90 

.F.o:
	$(CPP) $(CPPFLAGS) $(CPPARGS) $*.F > $*.f90
#	$(FC) -c $(FFLAGS) $*.f90 > $*.o
	$(FC) -c $(FFLAGS_BIO) $*.f90 > $*.o
#	$(AR) $(ARFLGS) $@ $<	
	ar rv libbiomodule.a *.o
	\rm $*.f90

clean:
	/bin/rm -f lib*.a  *.mod *.o *.f90

clobber:
		make clean
		/bin/rm -f *.f90

Makefiles::

#includes::
#include ../make.inc
# DO NOT DELETE
