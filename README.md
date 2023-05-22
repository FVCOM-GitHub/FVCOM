# FVCOM 5.0.1

http://fvcom.smast.umassd.edu/




## Code download

To download the latest version of FVCOM:<br>
git clone https://github.com/FVCOM-GitHub/FVCOM.git

To obtain an old version of FVCOM:<br>
git clone --branch <strong>VERSION</strong> https://github.com/FVCOM-GitHub/FVCOM.git<br>
where <strong>VERSION</strong> is the version number. All available verions can be found on https://github.com/FVCOM-GitHub/FVCOM/releases.


## Required libraries

Compilers
FVCOM codes are mainly written in Fortran 90 and C language. We recommend users to use 
* ifort and icc. (https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
* gfortran and gcc. (https://gcc.gnu.org/)

Required/Optional libraries
*	NetCDF: to read input files and write output files in NetCDF format. Both NetCDF-C and NetCDF-Fortran are required. (https://www.unidata.ucar.edu/software/netcdf/)
*	MPI: (optional) for parallel simulation with multiple cpus. Many options are freely available online, such as
    *	MPICH (https://www.mpich.org/downloads/)
    *	OneAPI (https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
*	Metis: (optional) for serial graphic partition. (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
*	Julian: for calendar calculation. (https://pds-rings.seti.org/toolkits/)
*	Proj4: used for coordinate conversion. (https://proj.org/)
*	Fortran-proj4: bindings for proj4 in Fortran. (not available online.) 
*	PetSc: (optional) the toolkit of scientific computation. This library is required when using semi-implicit scheme, data assimilation, non-hydrostatic, or SWAVE module.(https://petsc.org/)
* ESMF: (optional) This library is required when online nesting or WRF-FVCOM coupling is applied.


## Installation
*	make depends<br>
Generate the module dependencies of FVCOM codes. This step is required only when there are new codes added.
*	edit make.inc<br>
The ‘make.inc’ file includes all the settings for compiling FVCOM. This file needs to be edited correctly based on users’ interest of simulation, before the FVCOM codes are compiled. Users need to edit three parts in this file, including library environments, control flags, and compiler settings. 
*	make clean<br>
When the ‘make.inc’ file is modified, this step is necessary to make the model compiled from the beginning with the modified settings.
*	make<br>
Compile the FVCOM codes. You should see the executable file ‘fvcom’ when the compilation is successful.
FVCOM can be compiled in a parallel mode, with the following command:<br>
make -j N<br>
where N is the integer number specifying the maximum number of cores used for compiling the FVCOM source codes.
Based on the tests with Intel® Xeon® CPU E5-2640, the total time of compiling FVCOM is 120 s with one core. The time can be saved 40% with two cores and 54% with three cores. When more than three cores are applied, there is no large improvement on the compiling speed.
The results could vary with different cpus and different flags selected. However, we recommend users compile FVCOM with 2 or 3 cores to save the compiling time.


## Set up and run

* Step 1: make a folder ‘run’ and copy/link the executable file ‘fvcom’ to this folder.
* Step 2: prepare all required input files
* Step 3: create the namelist file with the name of CASENAME_run.nml. You can get a blank namelist file by<br>
      ./fvcom –create_namelist
* Step 4: run the model.<br>
    * To run FVCOM with single cpu<br>
        ./fvcom –casename=CASENAME
    * To run FVCOM in the parallel way<br>
        mpiexec ./fvcom –casename=CASENAME<br>
        or<br>
        mpiexec ./fvcom –casename=CASENAME<br>
      where CASENAME is the name of simulation case and must be consistent with the prefix of the namelist file. For example, the namelist file is named as ‘gom_run.nml’, then CASENAME is ‘gom’.


## Testsuite

We also provide a package of helping users to learn and run FVCOM, including benchmark test cases, offline models, necessary libraries, and processing tools. The package is available at:<br>
https://drive.google.com/file/d/1xwcFjzkSNT26FBu83pq2I8oIGzPqnjUn/view



