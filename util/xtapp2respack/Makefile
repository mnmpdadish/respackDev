######################################
# Copyright (c) 2018 Kazuma Nakamura #
###################################### 
TARGET = strconv wfn2respack 
#--------------------------------------------------------------------
# for Intel Linux and Intel Fortran Compiler 11.0 and above 
#--------------------------------------------------------------------
FC = ifort  
FFLAGS = -O2 -qopenmp -xHost -g -traceback

#######
#DEBUG#
#######
#FFLAGS = -qopenmp -debug full -check all -warn all -traceback -g 

LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 

#GNU
#FC = gfortran 
#FFLAGS = -O2 -fopenmp -g -fbacktrace
#LIBBLAS = -llapack -lblas -fopenmp 
#--------------------------------------------------------------------
.SUFFIXES:
.SUFFIXES: .o .F90 .F .f90 .c
# default suffix rules. Do not edit this area.
.F.o:
	${FC} ${F90FLAGS} ${FFLAGS} -c $<
.f90.o:
	${FC} ${F90FLAGS} ${FFLAGS} -c $<
.F90.o:
	${FC} ${F90FLAGS} ${FFLAGS} -c $<
#--------------------------------------------------------------------

STRCONVOBJ = strconv.o subr_fmtconv.o subr_readstr.o  

WFN2RESPACKOBJ = wfn2respack.o 

MOD_FILES = subr_fmtconv.mod subr_readstr.mod  

all: strconv wfn2respack 

strconv: ${STRCONVOBJ} 
	 ${FC} ${LDFLAGS} -o $@ ${STRCONVOBJ} ${LIBBLAS} 

wfn2respack: ${WFN2RESPACKOBJ} 
	     ${FC} ${LDFLAGS} -o $@ ${WFN2RESPACKOBJ} ${LIBBLAS} 

strconv.o: strconv.F90 subr_fmtconv.o subr_readstr.o  

clean: 
	${RM} ${TARGET} ${STRCONVOBJ} ${WFN2RESPACKOBJ} ${MOD_FILES} *.mod *__genmod.f90* 
#--------------------------------------------------------------------
