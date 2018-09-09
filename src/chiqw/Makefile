######################################
# Copyright (c) 2017 Kazuma Nakamura #
###################################### 
#--------------------------------------------------------------------
# for Intel Linux and Intel Fortran Compiler 11.0 and above 
#--------------------------------------------------------------------
#intel
FC = mpiifort 
FFLAGS = -O2 -qopenmp -xHost -mcmodel=large -shared-intel 

#######
#DEBUG#
#######
#FFLAGS = -qopenmp -xHost -mcmodel=large -shared-intel -debug full -check all -warn all -traceback -g  

#############################
# note by Maxime Charlebois #
#############################
#
# libiomp5 and libpthread are not needed here, conflict with MKL on some system, better not use it.
#
#LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 
LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group 

#openMPI 
#FC = mpif90 
#FFLAGS = -O2 -qopenmp -xHost -mcmodel=large -shared-intel
#LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 

#sekirei 
#FC = mpif90 
#FFLAGS = -O2 -qopenmp -g -traceback
#LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 

#GNU 
#FC = mpif90 
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
TARGET = calc_chiqw 
OBJECTS = m_rdinput.o m_tetrahedron_20170325.o m_fft3d_20150826.o gen_grid.o est_NTK.o est_nwx2.o inv.o invZGE.o \
	  wrt_optical_property.o chiqw_sub.o chiqw.o 

MOD_FILES = m_rdinput.mod fft_3d.mod m_tetrahedron.mod 

${TARGET} : ${OBJECTS} 
	${FC} ${FFLAGS} -o $@ ${OBJECTS} ${LIBBLAS} 

chiqw.o: chiqw.F m_rdinput.o m_tetrahedron_20170325.o m_fft3d_20150826.o gen_grid.o est_NTK.o est_nwx2.o inv.o  \
         invZGE.o chiqw_sub.o wrt_optical_property.o  

chiqw_sub.o: chiqw_sub.F m_fft3d_20150826.o invZGE.o

est_NTK.o: est_NTK.f90 chiqw_sub.o

clean: 
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES} *.mod *__genmod.f90* 
#--------------------------------------------------------------------
