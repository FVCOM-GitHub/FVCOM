#
# Configuration file for sst2grd Code
#
            SHELL = /bin/sh


#--------------------------------------------------------------------------
#  NETCDF STUFF
#--------------------------------------------------------------------------

         EXEC     = xbilin xquad
         MYOPT      = #-static
         FVCOMINCS = -I$(TOPDIR)
         FVCOMLIBS = -L$(TOPDIR) -lfvcom

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .F 

.F.o:
	$(CPP) $(CPPARGS) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(MYINCS) $*.f90
	\rm $*.f90

#--------------------------------------------------------------------------
#  Filter Source Code.
#--------------------------------------------------------------------------

MYLIBS = $(FVCOMLIBS) $(LIBS)
MYINCS = $(FVCOMINCS) $(INCS)


BILIN  = bilinear.F
QUAD  = quad.F

#SRCS =  $(BILIN) $(QUAD) 

BOBJS = $(BILIN:.F=.o)
QOBJS = $(QUAD:.F=.o)

#--------------------------------------------------------------------------
#  Linking Directives               
#--------------------------------------------------------------------------

$(EXEC):	$(BOBJS) $(QOBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) $(MYOPT) $(BOBJS) $(MYLIBS) -o xbilin
		$(FC) $(FFLAGS) $(LDFLAGS) $(MYOPT) $(QOBJS) $(MYLIBS) -o xquad

#--------------------------------------------------------------------------
#  Target to create dependecies.
#--------------------------------------------------------------------------

depend:
		mv makedepends makedepends.orig
		makedepf90  $(SRCS) >> makedepends


#--------------------------------------------------------------------------
#  Tar Up Code                           
#--------------------------------------------------------------------------

#tarfile:
#		tar cvf interp.tar *.f90 makefile 

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod xbilin xquad

clobber:	clean
		/bin/rm -f *.f90
		( cd results ; /bin/rm -f *.nc)

#--------------------------------------------------------------------------
#  Common rules for all Makefiles - do not edit.
#--------------------------------------------------------------------------

emptyrule::

#--------------------------------------------------------------------------
#  Empty rules for directories that do not have SUBDIRS - do not edit.
#--------------------------------------------------------------------------

install::
	@echo "install in $(CURRENT_DIR) done"

install.man::
	@echo "install.man in $(CURRENT_DIR) done"

Makefiles::

includes::
include ../../make.inc
include ./makedepends
# DO NOT DELETE
