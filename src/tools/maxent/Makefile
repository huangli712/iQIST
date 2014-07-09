.SUFFIXES: .f90

include ./make.sys

default: all

modm = mod_constants.o mod_control.o mod_context.o
util = maxent_utils.o 
lev1 = maxent_stream.o maxent_model.o
lev2 = maxent_makeaw.o
lev3 = maxent_driver.o
dump = maxent_dump.o maxent_print.o
main = maxent_main.o

objects = $(modm) $(util) $(lev1) $(lev2) $(lev3) $(dump) $(main)

all: maxent

maxent: $(objects)
	$(LINKER) $(objects) -o maxent $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f maxent

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
