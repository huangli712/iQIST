# iQIST (Interacting Quantum Impurity Solver Toolkit)

### WARNING

The user manual for iQIST is far away from completeness. So please DO NOT READ manual/ug.pdf now. We are really sorry for that. Though the user manual is not ready, you can still find some useful information/tips in the following places.

### Installation

see the comments in iqist/src/build/make.sys, iqist/src/build/make.md, and iqist/src/build/Makefile.

### File Format

for the file formats of solver.ctqmc.in, atom.config.in, entropy.in, and sac.in, see the comments in iqist/src/common/m\_parser.f90.

for the file formats of solver.xxx.dat, see the codes in ctqmc\_dump.f90.

### CTQMC quantum impurity solver

see the comments in iqist/src/ctqmc/xxx/ctqmc\_control.f90, here xxx means 'azalea', 'gardenia', 'narcissus', 'begonia', 'lavender', 'pansy', and 'manjushaka'.

### HFQMC quantum impurity solver

see the comments in iqist/src/hfqmc/daisy/hfqmc\_control.f90.

### The JASMINE component

see the comments in iqist/src/tools/jasmine/atomic\_control.f90.

### The HIBISCUS/entropy component

see the comments in iqist/src/tools/hibiscus/entropy/entropy\_control.f90.

### The HIBISCUS/script component

see the comments in the scripts.

### The HIBISCUS/stoch component

see the comments in iqist/src/tools/hibiscus/stoch/sac\_control.f90.

### The HIBISCUS/swing component

see the comments in iqist/src/tools/hibiscus/swing/swing\_main.py.

### Application programming interface

see the comments in iqist/src/api/ctqmc\_api.f90, iqist/src/api/hfqmc\_api.f90, and iqist/src/api/atomic\_api.f90

### Examples

benchmark examples: see iqist/working.

tutorial examples: see iqist/tutor.
