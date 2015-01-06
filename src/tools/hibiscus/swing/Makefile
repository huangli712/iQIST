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
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 01/05/2015 by li huang
##
##

default: all

all: exe

exe:
	f2py -c --fcompiler=intelem --f90flags="-O3 -xHost" swing_fast.f90 -m swing_fast

clean:
	rm -f *.so
	rm -f *.pyc

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
