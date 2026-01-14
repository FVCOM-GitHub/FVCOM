# Building FVCOM on macOS (Apple Silicon / Intel)

This guide describes how to build FVCOM on macOS using gfortran and OpenMPI.

## Prerequisites

Install required dependencies via Homebrew:

```bash
brew install gcc open-mpi netcdf netcdf-fortran hdf5 metis
```

## Configuration (make.inc)

Edit `src/make.inc` with the following changes:

### 1. Set paths

```makefile
TOPDIR = /path/to/FVCOM/src
INSTALLDIR = /path/to/FVCOM/libs/install
```

### 2. Enable LOCAL INSTALL

```makefile
LIBDIR = -L$(INSTALLDIR)/lib
INCDIR = -I$(INSTALLDIR)/include
```

### 3. Configure NetCDF (Homebrew paths)

```makefile
IOLIBS = -L/opt/homebrew/lib -lnetcdff -lnetcdf
IOINCS = -I/opt/homebrew/include
```

### 4. Configure METIS (Homebrew paths)

```makefile
FLAG_411 = -DMETIS_5
PARLIB =
PARTINCS = -I/opt/homebrew/include
PARTLIBS = -L/opt/homebrew/lib -lmetis
```

### 5. Use gfortran/OpenMPI compiler

Comment out the Intel compiler section and add:

```makefile
CPP      = /usr/bin/cpp
COMPFLAG = -DGFORTRAN
CC       = mpicc
CXX      = mpicxx
CFLAGS   = -O3
FC       = mpif90
DEBFLGS  =
OPT      = -O3 -ffree-line-length-none
CLIB     =
```

## Building the Julian Library

The Julian library requires a fix for Fortran name mangling on macOS with gfortran.

### Fix fortran.h

Edit `libs/julian/fortran.h` and remove the `__APPLE__` special case for `FORTRAN_NAME`:

**Before:**
```c
#ifdef __APPLE__
#define FORTRAN_NAME(name)  name

#else

#ifdef __STDC__
#define FORTRAN_NAME(name)  name##_
...
```

**After:**
```c
#ifdef __STDC__
#define FORTRAN_NAME(name)  name##_

#else
#define FORTRAN_NAME(name)  name/**/_

#endif
#endif
```

This ensures gfortran's expected trailing underscore convention is used.

### Build and install Julian

```bash
cd libs/julian
make clean
make libjulian

# Install to libs/install
mkdir -p ../install/lib ../install/include
cp libjulian.a ../install/lib/
cp fjulian.inc ../install/include/
```

## Building FVCOM

```bash
cd src

# Clean previous build
make clean

# Manually compile partition.c (workaround for makefile issue)
mpicc -c -O3 -I/opt/homebrew/include partition.c

# Build FVCOM
make
```

## Known Issues

### 1. partition.c not compiled automatically

The makefile's `ifdef FLAG_411` check occurs before `make.inc` is included, so `partition.c` may not be compiled. Compile it manually as shown above.

### 2. Parallel make race condition

Using `make -j` may cause race conditions with module dependencies. If you encounter "Cannot open module file" errors, run `make` without the `-j` flag.

## Verification

Test the build:

```bash
./fvcom --help
```

You should see the FVCOM help message with available options.
