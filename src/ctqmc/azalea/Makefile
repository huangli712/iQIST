.SUFFIXES: .f90

include ../../build/make.sys

default: all

modm = mod_mpi.o mod_spring.o mod_stack.o
modc = ctqmc_constants.o ctqmc_control.o ctqmc_context.o
modi = ctqmc_api.o
dmft = ctqmc_dmft.o
core = ctqmc_solver.o
lev1 = ctqmc_flavor.o ctqmc_update.o
lev2 = ctqmc_record.o ctqmc_status.o ctqmc_stream.o
lev3 = ctqmc_fourier.o ctqmc_spline.o ctqmc_util.o
dump = ctqmc_dump.o ctqmc_print.o
main = ctqmc_main.o

objects = $(modm) $(modc) $(modi) $(dmft) $(core) $(lev1) $(lev2) $(lev3) $(dump) $(main)

all: ctqmc_exe

cssl:

csml:
	cp ../../common/constants.mod
	cp ../../common/mmpi.mod
	cp ../../common/spring.mod
	cp ../../common/stack.mod

ctqmc_exe: $(objects)
	$(LINKER) $(objects) -o ctqmc $(LFLAGS) $(LIBS)

ctqmc_lib: $(objects)
	$(ARCHIVER) libctqmc.a $(objects)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f ctqmc
	rm -f libctqmc.a

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
