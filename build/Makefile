##
##
## Introduction
## ============
##
## This is the top Makefile for the iQIST software package. You can use
## it to build all of the executable programs and libraries. Please do
## not modify it by yourself even you are very familiar with iQIST.
##
## Usage
## =====
##
## 1. get help
## -----------
##
## 'make help' or 'make help-more'
##
## 2. quick build
## --------------
##
## 'make all'
##
## 3. build impurity solver
## ------------------------
##
## 'make solver'
##
## 4. build auxiliary tools
## ------------------------
##
## 'make tool'
##
## 5. build specified code
## -----------------------
##
## 'make xxx'
##
## here xxx is the code name, please type 'make help-more' to obtain more
## information about the possible values of xxx
##
## 6. clean the directories
## ------------------------
##
## 'make clean'
##
## Author
## ======
##
## This building system is designed, created, implemented, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 02/28/2015 by li huang (created)
## 03/30/2017 by li huang (last modified)
##
##

# targets for help
help:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print short usage information                  '
	@echo '     help-more         print full usage information                   '
	@echo
	@echo '     all               build all components, apis, and tools, etc     '
	@echo '     lib               build all components as libs                   '
	@echo '     clean             clean all directories                          '
	@echo
	@echo '     solver            build all ctqmc components                     '
	@echo '     solver-lib        build all ctqmc components as libs             '
	@echo '     clean-solver      clean objects for impurity solvers             '
	@echo
	@echo '     tool              build auxiliary tool components                '
	@echo '     tool-lib          build auxiliary tool components as libs        '
	@echo '     clean-tool        clean objects for auxiliary tools              '

help-more:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print short usage information        '
	@echo '     help-more         print full usage information         '
	@echo
	@echo '     narcissus         build narcissus code                 '
	@echo '     narcissus-lib     build narcissus library for Fortran  '
	@echo '     narcissus-pylib   build narcissus library for Python   '
	@echo '     clean-narcissus   clean narcissus directory            '
	@echo
	@echo '     manjushaka        build manjushaka code                '
	@echo '     manjushaka-lib    build manjushaka library for Fortran '
	@echo '     manjushaka-pylib  build manjushaka library for Python  '
	@echo '     clean-manjushaka  clean manjushaka directory           '
	@echo
	@echo '     jasmine           build jasmine code                   '
	@echo '     jasmine-lib       build jasmine library for Fortran    '
	@echo '     jasmine-pylib     build jasmine library for Python     '
	@echo '     clean-jasmine     clean jasmine directory              '
	@echo
	@echo '     hibiscus          build hibiscus code                  '
	@echo '     hibiscus-lib      build hibiscus code                  '
	@echo '     hibiscus-pylib    build hibiscus code                  '
	@echo '     clean-hibiscus    clean hibiscus directory             '
	@echo
	@echo '     capi              build capi library                   '
	@echo '     capi-lib          build capi library                   '
	@echo '     clean-capi        clean capi directory                 '
	@echo
	@echo '     base              build base library                   '
	@echo '     base-lib          build base library                   '
	@echo '     clean-base        clean base directory                 '

# all-in-one building target
all: base capi solver tool

# all-in-one building target (lib mode)
lib: base-lib capi-lib solver-lib tool-lib

# all-in-one cleaning target
clean: clean-base clean-capi clean-solver clean-tool

# targets for solver
solver: base capi narcissus manjushaka
solver-lib: base-lib capi-lib narcissus-lib manjushaka-lib
clean-solver: clean-base clean-capi clean-narcissus clean-manjushaka

# targets for tool
tool: base capi jasmine hibiscus
tool-lib: base-lib capi-lib jasmine-lib hibiscus-lib
clean-tool: clean-base clean-capi clean-jasmine clean-hibiscus

# targets for ctqmc codes: segment version (narcissus)
narcissus:
	cd ../src/ctqmc/narcissus/; pwd; make

narcissus-lib:
	cd ../src/ctqmc/narcissus/; pwd; make lib

narcissus-pylib:
	cd ../src/ctqmc/narcissus/; pwd; make pylib

clean-narcissus:
	cd ../src/ctqmc/narcissus/; pwd; make clean

# targets for ctqmc codes: general version (manjushaka)
manjushaka:
	cd ../src/ctqmc/manjushaka/; pwd; make

manjushaka-lib:
	cd ../src/ctqmc/manjushaka/; pwd; make lib

manjushaka-pylib:
	cd ../src/ctqmc/manjushaka/; pwd; make pylib

clean-manjushaka:
	cd ../src/ctqmc/manjushaka/; pwd; make clean

# targets for atomic codes: jasmine
jasmine:
	cd ../src/tools/jasmine/; pwd; make

jasmine-lib:
	cd ../src/tools/jasmine/; pwd; make lib

jasmine-pylib:
	cd ../src/tools/jasmine/; pwd; make pylib

clean-jasmine:
	cd ../src/tools/jasmine/; pwd; make clean

# targets for toolkit codes: hibiscus
hibiscus:
	cd ../src/tools/hibiscus/; pwd; make

hibiscus-lib:
	cd ../src/tools/hibiscus/; pwd; make

hibiscus-pylib:
	cd ../src/tools/hibiscus/; pwd; make

clean-hibiscus:
	cd ../src/tools/hibiscus/; pwd; make clean

# targets for capi: capi
capi:
	cd ../src/capi/; pwd; make

capi-lib:
	cd ../src/capi/; pwd; make

clean-capi:
	cd ../src/capi/; pwd; make clean

# targets for base: base
base:
	cd ../src/base/; pwd; make

base-lib:
	cd ../src/base/; pwd; make

clean-base:
	cd ../src/base/; pwd; make clean
