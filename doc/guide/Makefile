##
##
## Introduction
## ============
##
## Makefile for iqist/doc/guide directory.
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
## 10/06/2014 by li huang
##
##

default: all

all: fast_build

fast_build:
	xelatex ug

full_build:
	xelatex ug
	bibtex ug
	xelatex ug
	xelatex ug

clean:
	rm -f *.aux *.bbl *.blg *.lof *.log *.lot *.out *.toc
	rm -fr chapter/*.aux
