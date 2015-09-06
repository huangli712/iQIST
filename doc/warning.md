# iQIST (Interacting Quantum Impurity Solver Toolkit)

### WARNING

The user manual for iQIST is far away from completeness. So please DO NOT READ manual/ug.pdf now. We are really sorry for that. Though the user manual is not ready, you can still find some useful information/tips in the following places.

### Installation

See the comments in iqist/src/build/make.sys, iqist/src/build/build.md, and iqist/src/build/Makefile.

### File Format

For the file formats of solver.ctqmc.in, atom.config.in, entropy.in, and sac.in, see the comments in iqist/src/common/m\_parser.f90.

For the file formats of solver.xxx.dat, see the codes in ctqmc\_dump.f90.

### CTQMC quantum impurity solver

See the comments in iqist/src/ctqmc/xxx/ctqmc\_control.f90, here xxx means 'azalea', 'gardenia', 'narcissus', 'begonia', 'lavender', 'pansy', and 'manjushaka'.

### HFQMC quantum impurity solver

See the comments in iqist/src/hfqmc/daisy/hfqmc\_control.f90.

### The JASMINE component

See the comments in iqist/src/tools/jasmine/atomic\_control.f90.

### The HIBISCUS/entropy component

See the comments in iqist/src/tools/hibiscus/entropy/entropy\_control.f90.

### The HIBISCUS/script component

See the comments in the scripts.

### The HIBISCUS/stoch component

See the comments in iqist/src/tools/hibiscus/stoch/sac\_control.f90.

### The HIBISCUS/swing component

See the comments in iqist/src/tools/hibiscus/swing/swing\_main.py.

### Application programming interface

See the comments in iqist/src/capi/ctqmc\_api.f90, iqist/src/capi/hfqmc\_api.f90, and iqist/src/capi/atomic\_api.f90

### Examples

Benchmark examples: see iqist/working.

Tutorial examples: see iqist/tutor.
