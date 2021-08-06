##
##
## Introduction
## ============
##
## This is the official configuration file for the building system. You
## should modify it to fulfill your requirements. The make.inc file is
## the key component of the building system. If it is not configured
## correctly, the building system won't work correctly as well.
##
## Author
## ======
##
## This building system is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 05/11/2015 by li huang (created)
## 08/06/2021 by li huang (last modified)
##
##

# Fortran compiler, linker, and archiver
#-------------------------------------------------------------------------
F90    = mpif90
LINKER = $(F90)
ARCHIVER = ar -ruv

# Fortran preprocessor options
#-------------------------------------------------------------------------
MPI    = -DMPI
OMP    = #-qopenmp
FPP    = -fpp
CPP    = $(FPP) $(MPI) $(OMP)

# Machine tuning options
#-------------------------------------------------------------------------
CHECK  = -nogen-interfaces -warn all #-check all -traceback -g
MTUNE  = -O3 -xHost

# Flags for fortran compiler and linker
#-------------------------------------------------------------------------
FFLAGS = -c $(CPP) $(CHECK) $(MTUNE)
LFLAGS = $(OMP)

# External linear algebra library
#-------------------------------------------------------------------------
LIBS   = -L/Users/lihuang/opt/intel/mkl/lib -lmkl_core -lmkl_sequential -lmkl_rt
FLINK  = /Users/lihuang/Working/devel/flink/src
