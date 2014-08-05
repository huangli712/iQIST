.SUFFIXES: .f90

include ./make.sys

default: all

modm = mod_mpi.o mod_spring.o
modn = sai_constants.o sai_control.o sai_context.o
core = sai_record.o sai_stream.o sai_toolbox.o sai_update.o
dump = sai_dump.o sai_print.o
main = sai_main.o

objects = $(modm) $(modn) $(core) $(dump) $(main)

all: sai

sai: $(objects)
	$(LINKER) $(objects) -o sai $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f sai

clean-dat:
	rm -f *.dat
	rm -f *.bin
	rm -f *.out

clean-all: clean clean-dat
