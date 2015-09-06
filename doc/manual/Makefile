##
##
## Introduction
## ============
##
## Makefile for iqist/doc/manual directory.
##
## Usage
## =====
##
## 'make fast' or 'make full' or 'make clean'
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
## 10/06/2014 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

default: all

all: fast

fast:
	xelatex ug

full:
	xelatex ug
	bibtex ug
	xelatex ug
	xelatex ug

clean:
	rm -f *.aux *.bbl *.blg *.lof *.log *.lot *.out *.toc
	rm -fr chapter/*.aux
