#/===========================================================================/
# Copyright (c) 2007, The University of Massachusetts Dartmouth 
# Produced at the School of Marine Science & Technology 
# Marine Ecosystem Dynamics Modeling group
# All rights reserved.
#
# FVCOM has been developed by the joint UMASSD-WHOI research team. For 
# details of authorship and attribution of credit please see the FVCOM
# technical manual or contact the MEDM group.
#
# 
# This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
# The full copyright notice is contained in the file COPYRIGHT located in the 
# root directory of the FVCOM code. This original header must be maintained
# in all distributed versions.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  
#
#/---------------------------------------------------------------------------/
# CVS VERSION INFORMATION
# $Id$
# $Name$
# $Revision$
#/===========================================================================/

#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh
            EXEC          = fvcom

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90 

.c.o:
	$(CC) -c $(CFLAGS) $(INCS) $*.c

.F.o:
	$(CPP) $(CPPARGS) $(PETSC_FC_INCLUDES) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(INCS) $*.f90
	\rm $*.f90

.F90.o:
	$(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCLDIR) $< > $*.f  
	$(FC) -c $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR) $*.f  
	\rm $*.f90


#--------------------------------------------------------------------------
#  FVCOM Source Code.
#--------------------------------------------------------------------------

# Note: Change the order of modules at your own risk!

MODS  = mod_prec.F	sinter.F	mod_types.F	mod_time.F	\
	mod_main.F	mod_spherical.F	mod_utils.F	mod_clock.F	\
	eqs_of_state.F	mod_bulk.F	mod_interp.F	mod_par.F	\
	mod_par_special.F		mod_ncll.F	mod_nctools.F	\
	mod_wd.F	mod_sng.F	mod_heatflux.F	mod_solar.F	\
	mod_input.F	mod_force.F	mod_obcs.F	mod_petsc.F	\
	mod_semi_implicit.F		mod_non_hydro.F	mod_set_time.F	\
	ice_kinds_mod.F	ice_model_size.F		ice_domain.F	\
	ice_constants.F	ice_fileunits.F	ice_state.F	ice_work.F	\
	ice_grid.F	ice_albedo.F	ice_calendar.F	ice_flux.F	\
	ice_flux_in.F	ice_itd.F	ice_mechred.F	ice_itd_linear.F\
	ice_scaling.F	ice_atmo.F	ice_ocean.F	ice_therm_vertical.F\
	ice_init.F	ice_therm_itd.F	mod_ice2d.F	mod_ice.F	\
	mod_startup.F	mod_wqm.F	mod_ncdio.F	mod_setup.F	\
	mod_newinp.F	particle.F	linklist.F	mod_lag.F	\
	mod_northpole.F mod_pwp.F	mod_dye.F	mod_optimal_interpolation.F\
	mod_report.F	mod_probe.F	mod_gotm.F	mod_balance_2d.F\
	mod_tridiag.F	mod_scal.F	mod_meanflow.F	mod_obcs2.F	\
	mod_obcs3.F	mod_sed.F	mod_enkf.F	mod_etkf.F	\
	mod_rrk.F	mod_rrkf_obs.F	mod_rrkassim.F	mod_enkf_ncd.F	\
	enkf_ncdio.F	mod_enkf_obs.F	mod_enkfassim.F mod_assim.F	\
	mod_nesting.F	mod_visit.F     mod_plbc.F mod_dam.F		\
	mod_station_timeseries.F	mod_sparse_timeseries.F 	\
	mod_boundschk.F mod_esmf_nesting.F	\
	mod_cstms_vars.F	mod_flocmod.F	mod_sed_cstms.F		\
	mod_fluid_mud.F 

MODS_W = mod_main_wave.F	swmod1.F	swmod2.F	swmod3.F\
	mod_action_im.F mod_action_ex.F mod_wavesetup.F 		\
	mod_wave_current_interaction.F  mod_bbl.F

MAIN  =	fvcom.F		genmap.F	tge.F		longshore_flow.F\
	cell_area.F 	open_all.F	load_grid.F	allocate_all.F	\
	setup_domain.F	genmap_obc.F	genmap_lsf.F	print_vals.F    \
	coords_n_const.F	shape_coef_gcn.F shape_coef_gcy.F	\
	depth_grad.F	grid_metrics.F	cntrl_prmtrs.F	init_sed.F	\
	internal_step.F	bcond_gcn.F	bcond_gcy.F	adjust2d3d.F	\
	brough.F	advection_edge_gcn.F	advection_edge_gcy.F	\
	ghostuv.F	advave_edge_gcn.F	advave_edge_gcy.F	\
	phy_baropg.F	baropg.F	external_step.F	extel_edge.F	\
	extuv_edge.F	depth_check.F	vertvl_edge.F	adv_uv_edge_gcn.F\
	adv_uv_edge_gcy.F  vdif_uv.F    extelpf_edge.F	wreal.F		\
	viscofh.F	adv_q.F		fct_q2.F	fct_q2l.F	\
	vdif_q.F	adv_t.F		adv_s.F		fct_t.F		\
	vdif_ts.F	vdif_ts_gom.F	fct_s.F		bcond_ts.F	\
	adjust_ts.F	conv_over.F	visitsim.F	startup_type.F	\
	edge_len.F	adcor.F		icing.F		rho_pmean.F	\
	calc_vort.F	namelist.F	nh_set_nesting.F		\
	coare26z.F	coare40vn.F

MAIN_W = ocpmix.F	ocpcre.F	ocpids.F	\
	swanpre1.F	swanpre2.F	swanser.F	swanmain.F	\
	swancom1.F	swancom2.F	swancom3.F	swancom4.F	\
	swancom5.F	w3part.F

COBJS = func_pointer.c


BIOGEN	= mod_bio_3D.F	mod_onedtide.F

SRCS = $(MODS) $(MODS_W) $(MAIN) $(BIOGEN) $(MAIN_W)

# object list is composed of free (*.F) and fixed (*.F90) source
FVCOM_OBJS = $(COBJS:.c=.o) $(SRCS:.F=.o) $(FIXED:.F90=.o) 


#--------------------------------------------------------------------------
#  Linking Directives
#--------------------------------------------------------------------------

$(EXEC):	$(FVCOM_OBJS)
		$(FC) $(FFLAGS) $(VISOPT) $(LDFLAGS) $(FVCOM_OBJS) $(LIBS) -o $(EXEC)
		make libfvcom


libfvcom: $(FVCOM_OBJS)
	$(AR) $@.a $(FVCOM_OBJS)
	$(RANLIB) $@.a

#	$(FC) -shared -o $@.so $(FVCOM_OBJS)
#       PROBLEMS WITH SHARED OPTION, REQUIRES -fPIC in All sub libraries as well when building on x86_64


#--------------------------------------------------------------------------
#  Target to create dependecies.
#--------------------------------------------------------------------------

depend:
		makedepf90 -d mod_prec.o $(SRCS) >> makedepends

#--------------------------------------------------------------------------
#  Target to make all.
#--------------------------------------------------------------------------

SUBDIRS = testing input utilities


all:    
	make depend
	cd BIO_source && make
	make 
	for dir in $(SUBDIRS); do (cd $$dir && make) || exit 1; done

#--------------------------------------------------------------------------
#  Tar Up Code                           
#--------------------------------------------------------------------------

tarfile:
	tar cvf fvcom.tar *.F   makefile  makedepends RELEASE_NOTES 

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod *.f90

clobber:	
		make clean
		/bin/rm -f fvcom libfvcom.so libfvcom.a

allclean:
		make clean
		for dir in $(SUBDIRS); do (cd $$dir && make clean) || exit 1; done
	        cd BIO_source && make clean


allclobber:
		make clean
		/bin/rm -f fvcom libfvcom.so libfvcom.a
		for dir in $(SUBDIRS); do (cd $$dir && make clobber) || exit 1; done
	        cd BIO_source && make clobber

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
include ./makedepends
include ./make.inc
