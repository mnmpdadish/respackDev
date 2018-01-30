######################################
# Copyright (c) 2017 Kazuma Nakamura #
###################################### 
#--------------------------------------------------------------------
# for Intel Linux and Intel Fortran Compiler 11.0 and above 
#--------------------------------------------------------------------
#intel
FC = mpiifort 
FFLAGS = -O2 -openmp -xHost -mcmodel=large -shared-intel
LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 
#DEBUG  = -openmp -xHost -mcmodel=large -shared-intel -debug full -check all -warn all -traceback -g  

#openMPI 
#FC = mpif90 
#FFLAGS = -O2 -openmp -xHost -mcmodel=large -shared-intel
#LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 

#sekirei 
#FC = mpif90 
#FFLAGS = -O2 -qopenmp
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
#	${FC} ${F90FLAGS} ${DEBUG} -c $<
.f90.o:
	${FC} ${F90FLAGS} ${FFLAGS} -c $<
#       ${FC} ${F90FLAGS} ${DEBUG} -c $<
.F90.o:
	${FC} ${F90FLAGS} ${FFLAGS} -c $<
#	${FC} ${F90FLAGS} ${DEBUG} -c $<
#--------------------------------------------------------------------
TARGET = calc_chiqw 
OBJECTS = m_rdinput.o m_tetrahedron_20170325.o m_fft3d_20150826.o gen_grid.o est_NTK.o est_nwx2.o inv.o invZGE.o chiqw_sub.o chiqw.o 
MOD_FILES = m_rdinput.mod fft_3d.mod m_tetrahedron.mod 

${TARGET} : ${OBJECTS} 
	${FC} ${FFLAGS} -o $@ ${OBJECTS} ${LIBBLAS} 
#	${FC} ${DEBUG} -o $@ ${OBJECTS} ${LIBBLAS} 

chiqw.o: chiqw.F m_rdinput.o m_tetrahedron_20170325.o m_fft3d_20150826.o gen_grid.o est_NTK.o est_nwx2.o inv.o invZGE.o chiqw_sub.o 

chiqw_sub.o: chiqw_sub.F m_fft3d_20150826.o invZGE.o

clean: 
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES} 
#--------------------------------------------------------------------
