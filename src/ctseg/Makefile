##
##
## Introduction
## ============
##
## Makefile for iqist/src/ctseg directory.
##
## Usage
## =====
##
## 'make', 'make clean', 'make clean-dat' or 'make clean-all'
##
## Author
## ======
##
## This building system is designed, created, and maintained by
##
## Li Huang // email: huangli@caep.cn
##
## History
## =======
##
## 02/28/2015 by li huang (created)
## 04/30/2025 by li huang (last modified)
##
##

.SUFFIXES: .f90

define searchfile
    $(eval INC_DIR=$(patsubst %/,%,$(1)))
    $(eval INC_FILE=$(wildcard $(INC_DIR)/$(2)))
    $(info searchfile scanned: $(INC_DIR))
    $(if $(INC_FILE),
             $(info searchfile found: $(INC_FILE)),
                 $(call searchfile,$(dir $(INC_DIR)),$(2)))
endef

$(eval $(call searchfile,$(CURDIR),build/make.inc))

include $(INC_FILE)

define searchdir
    $(eval INP_DIR=$(patsubst %/,%,$(1)))
    $(eval FLINK=$(wildcard $(INP_DIR)/$(2)))
    $(info searchdir scanned: $(INP_DIR))
    $(if $(FLINK),
             $(info searchdir found: $(FLINK)),
                 $(call searchdir,$(dir $(INP_DIR)),$(2)))
endef

$(eval $(call searchdir,$(CURDIR),*/Flink/src))

modc = ctqmc_control.o ctqmc_context.o
dmft = ctqmc_dmft.o
core = ctqmc_solver.o
lev1 = ctqmc_flavor.o ctqmc_hybmat.o ctqmc_update.o
lev2 = ctqmc_record.o ctqmc_status.o ctqmc_stream.o ctqmc_util.o
lev3 = ctqmc_dump.o ctqmc_print.o
main = ctqmc_main.o
mlib = libflink.a

objects = $(modc) $(dmft) $(core) $(lev1) $(lev2) $(lev3) $(main) $(mlib)

default: all

all: exe

flink: flink_lib flink_mod

flink_lib:
	cp $(FLINK)/libflink.a .

flink_mod:
	cp $(FLINK)/constants.mod .
	cp $(FLINK)/mmpi.mod .
	cp $(FLINK)/face.mod .
	cp $(FLINK)/spring.mod .
	cp $(FLINK)/stack.mod .
	cp $(FLINK)/parser.mod .
	cp $(FLINK)/linkedlist.mod .

exe: flink $(objects)
	$(LINKER) $(objects) -o ctqmc $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f ctqmc
	rm -f libflink.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
