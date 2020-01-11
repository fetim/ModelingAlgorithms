#Master Makefile for Modeling

include Makefile.inc

all: compile

lib: compile

compile:
	cd C               ; $(MAKE) all
	cd C/openmp        ; $(MAKE) all
	cd C/openacc       ; $(MAKE) all
	cd fortran         ; $(MAKE) all
	cd fortran/openmp  ; $(MAKE) all
	cd fortran/openacc ; $(MAKE) all

clean:
	cd C               ; $(MAKE) clean
	cd C/openmp        ; $(MAKE) clean
	cd C/openacc       ; $(MAKE) clean	
	cd fortran         ; $(MAKE) clean
	cd fortran/openmp  ; $(MAKE) clean
	cd fortran/openacc ; $(MAKE) clean	