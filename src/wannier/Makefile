######################################
# Copyright (c) 2017 Kazuma Nakamura #
###################################### 
#--------------------------------------------------------------------
# for Intel Linux and Intel Fortran Compiler 11.0 and above 
#--------------------------------------------------------------------
#FC = ifort  
#FFLAGS = -O2 -qopenmp -xHost -g -traceback

#######
#DEBUG#
#######
#FFLAGS = -qopenmp -debug full -check all -warn all -traceback -g 

#LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -lm 

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
TARGET = calc_wannier 

OBJECTS = m_rdinput.o m_fft3d_20150826.o m_bvector_20171208.o det.o inv.o inv_ge_LAPACK.o \
	  diagV.o makekpts.o est_NTK.o est_nwx2.o est_latparam.o wrt_mvmc.o wrt_frmsf.o \
	  search_Rmin.o wannier_sub.o SO3_to_SU2_rotation.o wannier.o frmsf_wan.o gencif.o m_gencif_sub.o 

MOD_FILES = m_rdinput.mod fft_3d.mod m_bvector.mod frmsf_wan.mod m_gencif_sub.mod

${TARGET} : ${OBJECTS} 
	${FC} ${LDFLAGS} -o $@ ${OBJECTS} ${LIBBLAS} 

wannier.o: wannier.F m_rdinput.o m_fft3d_20150826.o m_bvector_20171208.o det.o inv.o inv_ge_LAPACK.o \
	   diagV.o makekpts.o est_NTK.o est_nwx2.o est_latparam.o wrt_mvmc.o wrt_frmsf.o \
	   search_Rmin.o wannier_sub.o frmsf_wan.o SO3_to_SU2_rotation.o gencif.o m_gencif_sub.o  

wannier_sub.o: wannier_sub.F m_fft3d_20150826.o diagV.o 

wrt_frmsf.o: wrt_frmsf.f90 wannier_sub.o 

frmsf_wan.o: frmsf_wan.F90 m_rdinput.o

est_NTK.o: est_NTK.f90 wannier_sub.o

wrt_mvmc.o: wrt_mvmc.f90 inv.o 

gencif.o: gencif.f90 m_gencif_sub.o  

clean: 
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES} *.mod *__genmod.f90* 
#--------------------------------------------------------------------
