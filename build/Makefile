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
## 'make help'
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
## 4. build atomic eigenvalue solver
## ---------------------------------
##
## 'make atomic'
##
## 5. build specified code
## -----------------------
##
## 'make xxx'
##
## here xxx is the code name, please type 'make help' to obtain more
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
## Li Huang // email: huangli@caep.cn
##
## History
## =======
##
## 02/28/2015 by li huang (created)
## 04/30/2025 by li huang (last modified)
##
##

# targets for help
help:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print full usage information        '
	@echo
	@echo '     all               build all components and tools      '
	@echo '     clean             clean all directories               '
	@echo
	@echo '     solver            build all quantum impurity solvers  '
	@echo '     clean-solver      clean all quantum impurity solvers  '
	@echo
	@echo '     atomic            build atomic eigenvalue solver      '
	@echo '     clean-atomic      clean atomic eigenvalue solver      '
	@echo
	@echo '     ctseg             build ctseg quantum impurity solver '
	@echo '     clean-ctseg       clean ctseg quantum impurity solver '
	@echo '     cthyb             build cthyb quantum impurity solver '
	@echo '     clean-cthyb       clean cthyb quantum impurity solver '

# all-in-one building target
all: solver atomic

# all-in-one cleaning target
clean: clean-solver clean-atomic

# targets for solver
solver: ctseg cthyb
clean-solver: clean-ctseg clean-cthyb

# targets for atomic eigenvalue problem solver: jasmine
atomic:
	cd ../src/atomic/; pwd; make

clean-atomic:
	cd ../src/atomic/; pwd; make clean

# targets for cthyb quantum impurity solver: segment version (narcissus)
ctseg:
	cd ../src/ctseg/; pwd; make

clean-ctseg:
	cd ../src/ctseg/; pwd; make clean

# targets for cthyb quantum impurity solver: general version (manjushaka)
cthyb:
	cd ../src/cthyb/; pwd; make

clean-cthyb:
	cd ../src/cthyb/; pwd; make clean
