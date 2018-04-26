######################################
# Copyright (c) 2017 Kazuma Nakamura #
###################################### 
#--------------------------------------------------------------------
# for Intel Linux and Intel Fortran Compiler 11.0 and above 
#--------------------------------------------------------------------
#intel
FC = mpiifort 
FFLAGS = -O2 -qopenmp -xHost -mcmodel=large -shared-intel 
LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group 
## note: libiomp5 and libpthread are not needed here, conflict with MKL on some system, better not use it.
#LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 
#DEBUG  = -openmp -xHost -mcmodel=large -shared-intel -debug full -check all -warn all -traceback -g  

#openMPI 
#FC = mpif90 
#FFLAGS = -O2 -openmp -xHost -mcmodel=large -shared-intel
#LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread 

#sekirei 
#FC = mpif90 
#FFLAGS = -O2 -openmp -g -traceback
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
TARGET = calc_gw 
OBJECTS = sub_gw.o sub_wfn.o sub_eps.o util.o sub_band_disp.o makekpts.o search_Rmin.o make_sgmw.o m_fft3d_20150826.o m_rd_input.o m_rd_dat_wfn.o m_rd_dat_wan.o m_rd_dat_eps.o gw.o
MOD_FILES = m_rdinput.mod fft_3d.mod m_tetrahedron.mod 

${TARGET} : ${OBJECTS} 
	${FC} ${FFLAGS} -o $@ ${OBJECTS} ${LIBBLAS} 
#	${FC} ${DEBUG} -o $@ ${OBJECTS} ${LIBBLAS} 

est_NTK.o: est_NTK.f90 chiqw_sub.o

clean: 
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES} 
#--------------------------------------------------------------------
