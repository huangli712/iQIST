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
## This building system is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 02/28/2015 by li huang (created)
## 04/21/2017 by li huang (last modified)
##
##

# targets for help
help:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print short usage information         '
	@echo '     help-more         print full usage information          '
	@echo
	@echo '     all               build all components and tools        '
	@echo '     clean             clean all directories                 '
	@echo
	@echo '     solver            build impurity solver components      '
	@echo '     clean-solver      clean objects for impurity solvers    '
	@echo
	@echo '     tool              build auxiliary tool components       '
	@echo '     clean-tool        clean objects for auxiliary tools     '

help-more:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print short usage information         '
	@echo '     help-more         print full usage information          '
	@echo
	@echo '     narcissus         build narcissus code                  '
	@echo '     clean-narcissus   clean narcissus directory             '
	@echo
	@echo '     manjushaka        build manjushaka code                 '
	@echo '     clean-manjushaka  clean manjushaka directory            '
	@echo
	@echo '     jasmine           build jasmine code                    '
	@echo '     clean-jasmine     clean jasmine directory               '
	@echo
	@echo '     hibiscus          build hibiscus code                   '
	@echo '     clean-hibiscus    clean hibiscus directory              '
	@echo
	@echo '     base              build base library                    '
	@echo '     clean-base        clean base directory                  '

# all-in-one building target
all: base solver tool

# all-in-one cleaning target
clean: clean-base clean-solver clean-tool

# targets for solver
solver: base narcissus manjushaka
clean-solver: clean-base clean-narcissus clean-manjushaka

# targets for tool
tool: base jasmine hibiscus
clean-tool: clean-base clean-jasmine clean-hibiscus

# targets for ctqmc codes: segment version (narcissus)
narcissus:
	cd ../src/ctqmc/narcissus/; pwd; make

clean-narcissus:
	cd ../src/ctqmc/narcissus/; pwd; make clean

# targets for ctqmc codes: general version (manjushaka)
manjushaka:
#	cd ../src/ctqmc/manjushaka/; pwd; make

clean-manjushaka:
#	cd ../src/ctqmc/manjushaka/; pwd; make clean

# targets for atomic codes: jasmine
jasmine:
#	cd ../src/tools/jasmine/; pwd; make

clean-jasmine:
#	cd ../src/tools/jasmine/; pwd; make clean

# targets for toolkit codes: hibiscus
hibiscus:
	cd ../src/tools/hibiscus/; pwd; make

clean-hibiscus:
	cd ../src/tools/hibiscus/; pwd; make clean

# targets for csml and cssl: base
base:
	cd ../src/base/; pwd; make

clean-base:
	cd ../src/base/; pwd; make clean
