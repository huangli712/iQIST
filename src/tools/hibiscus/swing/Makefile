##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/hibiscus/swing directory.
##
## Usage
## =====
##
## 'make', 'make clean', 'make clean-dat' or 'make clean-all'
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
## 05/08/2015 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

include ../../../../build/make.sys

default: all

all: exe

exe:
	$(F2PY) -c --f90exec=$(F2PYC) swing_fast.f90 -m swing_fast
	rm -fr *.dSYM

clean:
	rm -f *.so
	rm -f *.pyc

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
