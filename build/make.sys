##
##
## Introduction
## ============
##
## This is the official configuration file for the building system. You
## should modify it to fulfill your requirements. The make.sys file is
## the key component of the building system. If it is not configured
## correctly, the building system won't work correctly as well.
##
## Warning
## =======
##
## Be careful! This file may be not suitable for your situation. It is
## designed for my own hardware and software environment only. Please
## check and improve it before starting to compile your iQIST code.
##
## Explanations
## ============
##
## Please refer to the build.md file.
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
## 04/12/2019 by li huang (last modified)
##
##

# Fortran compiler, linker, and archiver
#-------------------------------------------------------------------------
F90    = /opt/local/mpich/bin/mpif90
LINKER = $(F90)
ARCHIVER = ar -ruv

# Fortran preprocessor options
#-------------------------------------------------------------------------
MPI    = -DMPI
OMP    = #-fopenmp
FPP    = -cpp
CPP    = $(FPP) $(MPI) $(OMP)

# Machine tuning options
#-------------------------------------------------------------------------
CHECK  = -Wall -Wunused -Wextra #-fbacktrace -fcheck=all -g
MTUNE  = -Ofast -faggressive-loop-optimizations -fno-tree-pre -fPIC -msse4

# Flags for fortran compiler and linker
#-------------------------------------------------------------------------
FFLAGS = -c $(CPP) $(CHECK) $(MTUNE)
LFLAGS = $(OMP)

# External linear algebra library
#-------------------------------------------------------------------------
LIBS   = -framework Accelerate
FLINK  = /Users/lihuang/Working/dmft/flink/src
