####################################################
### Author: Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
### Date: 06/07/2023
####################################################

########################################
### VERSIONS
PROJ                  :  4.9.3
   URL:  https://proj.org/download.html
   NOTE: Use the latest 4.* version
      
METIS                 :  5.1.0
   URL:  https://github.com/KarypisLab/METIS
   NOTE: See comments below

JULIAN                :  1.3.3
   URL:  https://pds-rings.seti.org/toolkits/
   NOTE: Select the one in C language, not the python one

FPROJ (proj4-fortran) :  1.0
   URL:  https://github.com/mhagdorn/proj4-fortran
   NOTE: Version 1.1.0 does not work with FVCOM; minor modifications to RFVCOM/src are needed
########################################


### Compile METIS 5.1.0 (DEFAULT VERSION metis.tgz -> metis-5.1.0.tgz)
1) untar the archive:
      tar zxf metis.tgz
   the above command will create a "metis" directory
2) configure/install metis:
      cd metis
      make cc=$CC cxx=$CXX openmp=ON gklib_path=${PWD}/GKlib prefix=${INSTALLDIR} config
      make install

NOTE: GKlib calls are included in the resulting libmetis.a

### Compile METIS 5.2.1
1) untar the archive:
      tar zxf metis-5.2.1.tgz
   the above command will create a "metis" directory
2) configure/install GKlib:
      cd metis/GKlib
      # The flag -D_POSIX_C_SOURCE=199309L is necessary to compile using the Intel compilers
      make cc=$CC cxx=$CXX openmp=ON CFLAGS="-D_POSIX_C_SOURCE=199309L ${CFLAGS}" prefix=${INSTALLDIR} config
      make install
3) configure/install metis:
      cd metis
      make cc=$CC cxx=$CXX openmp=ON gklib_path=${INSTALLDIR} prefix=${INSTALLDIR} CONFIG_FLAGS_EXTERN="-DCMAKE_C_FLAGS=-I${INSTALLDIR}/include -I${INSTALLDIR}/include/metis ${CFLAGS}" config
      make install

NOTE: GKlib calls are NOT included in the resulting libmetis.a so to link use:
      -L${INSTALLDIR}/lib -lmetis -lGKlib
