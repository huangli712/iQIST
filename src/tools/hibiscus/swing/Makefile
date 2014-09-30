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
## 'make' or 'make clean'
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
## 09/30/2014 by li huang
##
##

default: all

all: exe

exe: swing-fast
swing-fast:
	f2py -c --fcompiler=intelem --f90flags="-O3 -xHost" swing_fast.f90 -m swing_fast

clean:
	rm -f *.so
	rm -f *.pyc

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
