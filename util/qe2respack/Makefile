# Copyright (c) 2017 Yusuke Nomura, Terumasa Tadano 
# Simple Makefile for qe2respack

QEROOT=
# Root directory of a Quantum ESPRESSO distribution.
# Please edit QEROOT according to your environment.
# examples: 
# your machine
#   QEROOT=/home/xxx/src/QE/qe-6.1
# systemB@issp 
#   QEROOT=/home/kxxxx/kxxxx00/src/QE/qe-6.1    

# Intel compiler with intel MKL
FC = ifort
FFLAGS = -O2 -xHOST -g -traceback
LINKER =  mpif90 
LDFLAGS = -liomp5 -lpthread
LAPACK = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
# GNU
#FC = gfortran
#FFLAGS = -O3 -g -fbacktrace
#LINKER = mpif90
#LDFLAGS = -g -pthread
#LAPACK = -llapack -lblas

LFLAGS = ${LDFLAGS} ${LAPACK}
INCLUDE= -I${QEROOT}/iotk/include -I${QEROOT}/iotk/src 

#************************************************************

PROGRAMS = qe2respack

default: all

all: ${PROGRAMS}

qe2respack: qe2respack.o ${QEROOT}/iotk/src/libiotk.a
	${LINKER} $^ -o $@ ${LFLAGS}

qe2respack.o: qe2respack.f90
	${FC} ${FFLAGS} ${INCLUDE} $< -c

.PHONY: clean
clean:
	rm -f ${PROGRAMS} *.o
