######################################
# Copyright (c) 2018 Kazuma Nakamura #
###################################### 
#--------------------------------------------------------------------
# for Intel Linux and Intel Fortran Compiler 11.0 and above 
#--------------------------------------------------------------------
#intel
FC = mpiifort 
FFLAGS = -O2 -qopenmp -xHost -mcmodel=large -shared-intel -g -traceback

#######
#DEBUG#
#######
#FFLAGS = -qopenmp -debug full -check all -CB -warn all -traceback -g  

#lapack 
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
TARGET = calc_gw 

OBJECTS = m_rd_input.o m_rd_dat_wfn.o m_rd_dat_wan.o m_rd_dat_eps.o m_fft3d_20150826.o m_tetrahedron_20170325.o \
	  sub_wfn.o sub_eps.o sub_gw.o util.o calc_band_disp.o make_kpts.o make_sgmw.o calc_ksdos.o calc_gwdos.o calc_gwakw.o \
	  search_Rmin.o sub_ksdos.o sub_gwdos.o sub_eigenvalue.o sub_mkidx.o sub_det_shift.o gw.o \
	  m_qp.o m_FR.o m_band.o m_eigenstate.o m_bspline.o  

MOD_FILES = m_rd_input.mod m_rd_dat_wfn.mod m_rd_dat_wan.mod m_rd_dat_eps.mod fft_3d.mod m_tetrahedron.mod \
	    m_qp.mod m_FR.mod m_band.mod m_eigenstate.mod m_bspline.mod   

${TARGET} : ${OBJECTS} 
	${FC} ${FFLAGS} -o $@ ${OBJECTS} ${LIBBLAS} 

gw.o:   gw.f90 m_rd_input.o m_rd_dat_wfn.o m_rd_dat_wan.o m_rd_dat_eps.o m_fft3d_20150826.o m_tetrahedron_20170325.o \
	sub_wfn.o sub_eps.o sub_gw.o util.o calc_band_disp.o make_kpts.o make_sgmw.o calc_ksdos.o calc_gwdos.o calc_gwakw.o \
	search_Rmin.o sub_ksdos.o sub_gwdos.o sub_eigenvalue.o sub_mkidx.o sub_det_shift.o \
	m_qp.o m_FR.o m_band.o m_eigenstate.o m_bspline.o  

calc_band_disp.o: calc_band_disp.f90 search_Rmin.o util.o  

calc_gwakw.o: calc_gwakw.f90 search_Rmin.o sub_eigenvalue.o 

calc_ksdos.o: calc_ksdos.f90 search_Rmin.o sub_eigenvalue.o sub_ksdos.o

calc_gwdos.o: calc_gwdos.f90 search_Rmin.o sub_eigenvalue.o sub_det_shift.o sub_gwdos.o 

m_rd_dat_eps.o: m_rd_dat_eps.f90 sub_wfn.o sub_eps.o 

m_rd_dat_wfn.o: m_rd_dat_wfn.f90 util.o sub_wfn.o 

sub_eigenvalue.o: sub_eigenvalue.f90 util.o 

sub_gw.o: sub_gw.f90 m_fft3d_20150826.o 

sub_gwdos.o: sub_gwdos.f90 sub_mkidx.o m_tetrahedron_20170325.o 

sub_ksdos.o: sub_ksdos.f90 sub_mkidx.o m_tetrahedron_20170325.o

m_qp.o: m_qp.f90 m_FR.o m_band.o m_eigenstate.o m_bspline.o  

clean: 
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES} *.mod *__genmod.f90* 
#--------------------------------------------------------------------
