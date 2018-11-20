######################################
# Copyright (c) 2018 Kazuma Nakamura #
###################################### 
#--------------------------------------------------------------------
# for Intel Linux and Intel Fortran Compiler 11.0 and above 
#--------------------------------------------------------------------
#intel
FC = ifort  
FFLAGS = -O2 -qopenmp -xHost -g -traceback

#######
#DEBUG#
#######
#FFLAGS = -qopenmp -debug full -check all -CB -warn all -traceback -g  

#lapack 
LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group 

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
TARGET = transfer_analysis 

OBJECTS = m_rd_dat_zvo.o m_rd_transdef.o m_tetrahedron_20170325.o \
	  calc_dos.o calc_hist.o make_grd.o sub_search_Rmin.o  \
	  est_nkbi.o truncation.o make_eig.o \
	  make_kdata.o sub_eigenvalue.o sub_mkidx.o \
	  util.o est_ef.o wrt.o \
	  transfer_analysis.o  

MOD_FILES = m_rd_dat_zvo.mod m_rd_transdef.mod m_tetrahedron.mod 

${TARGET} : ${OBJECTS} 
	${FC} ${FFLAGS} -o $@ ${OBJECTS} ${LIBBLAS} 

transfer_analysis.o: transfer_analysis.f90 \
                     m_rd_dat_zvo.o m_rd_transdef.o m_tetrahedron_20170325.o \
	             calc_dos.o calc_hist.o make_grd.o sub_search_Rmin.o  \
	             est_nkbi.o truncation.o make_eig.o \
	             make_kdata.o sub_eigenvalue.o sub_mkidx.o \
	             util.o est_ef.o wrt.o 

m_rd_dat_zvo.o: m_rd_dat_zvo.f90 util.o est_nkbi.o 

calc_dos.o: calc_dos.f90 sub_mkidx.o m_tetrahedron_20170325.o

make_eig.o: make_eig.f90 sub_search_Rmin.o sub_eigenvalue.o 

sub_eigenvalue.o: sub_eigenvalue.f90 util.o 

clean: 
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES} *.mod *__genmod.f90* 
#--------------------------------------------------------------------
