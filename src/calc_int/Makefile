######################################
# Copyright (c) 2017 Kazuma Nakamura #
###################################### 
TARGET = calc_w3d calc_j3d 
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
FC = gfortran 
FFLAGS = -O2 -fopenmp -g -fbacktrace
LIBBLAS = -llapack -lblas -fopenmp 
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
COMMONOBJ = m_rdinput.o m_fft3d_20150826.o inv.o est_NTK.o est_nwx2.o calc_int_sub.o 
W3dOBJ = calc_intW.o wrt_model_w3d.o ${COMMONOBJ} 
J3dOBJ = calc_intJ.o wrt_model_j3d.o ${COMMONOBJ}
MOD_FILES = m_rdinput.mod fft_3d.mod 

all: calc_w3d calc_j3d 

calc_w3d: ${W3dOBJ} 
	  ${FC} ${LDFLAGS} -o $@ ${W3dOBJ} ${LIBBLAS} 

calc_j3d: ${J3dOBJ} 
	  ${FC} ${LDFLAGS} -o $@ ${J3dOBJ} ${LIBBLAS} 

calc_intW.o: calc_intW.F config_w3d.h m_rdinput.o m_fft3d_20150826.o inv.o est_NTK.o est_nwx2.o wrt_model_w3d.o \
             calc_int_sub.o 

calc_intJ.o: calc_intJ.F config_j3d.h m_rdinput.o m_fft3d_20150826.o inv.o est_NTK.o est_nwx2.o wrt_model_j3d.o \
             calc_int_sub.o 

calc_int_sub.o: calc_int_sub.F m_fft3d_20150826.o 

est_NTK.o: est_NTK.f90 calc_int_sub.o

clean: 
	${RM} ${TARGET} ${W3dOBJ} ${J3dOBJ} ${COMMONOBJ} ${MOD_FILES} *.mod *__genmod.f90* 
#--------------------------------------------------------------------
