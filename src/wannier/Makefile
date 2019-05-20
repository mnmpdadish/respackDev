######################################
# Copyright (c) 2017 Kazuma Nakamura #
###################################### 
#--------------------------------------------------------------------
# for Intel Linux and Intel Fortran Compiler 11.0 and above 
#--------------------------------------------------------------------
FC = ifort  
FFLAGS = -O2 -qopenmp -xHost -g -traceback
#FFLAGS = -O2 -openmp -xHost -g -traceback

#######
#DEBUG#
#######
#FFLAGS = -qopenmp -debug full -check all -warn all -traceback -g 

LIBBLAS = -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -lm 

#GNU
#FC = gfortran 
##FFLAGS = -O2 -fopenmp -g -fbacktrace
#FFLAGS = -fopenmp -debug -check -warn -fbacktrace -g 
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
TARGET = calc_wannier 

OBJECTS = m_rdinput.o m_fft3d_20150826.o m_bvector_20171208.o det.o inv.o inv_ge_LAPACK.o \
	  diagV.o makekpts.o est_NTK.o est_nwx2.o est_latparam.o m_wrt_mvmc.o m_wrt_frmsf.o \
	  search_Rmin.o wannier_sub.o wannier.o m_frmsf_wan.o gencif.o m_gencif_sub.o \
          m_dos.o m_tetrahedron_20170325.o m_eigenstate.o \
          m_dmx.o m_tetrainteg.o 

MOD_FILES = m_rdinput.mod fft_3d.mod m_bvector.mod m_frmsf_wan.mod m_gencif_sub.mod \
            m_dos.mod m_tetrahedron_20170325.mod m_eigenstate.mod \
            m_wrt_mvmc.mod m_wrt_frmsf.mod \
            m_dmx.mod m_tetrainteg.mod 

${TARGET} : ${OBJECTS} 
	${FC} ${LDFLAGS} -o $@ ${OBJECTS} ${LIBBLAS} 

wannier.o: wannier.F m_rdinput.o m_fft3d_20150826.o m_bvector_20171208.o det.o inv.o inv_ge_LAPACK.o \
	   diagV.o makekpts.o est_NTK.o est_nwx2.o est_latparam.o m_wrt_mvmc.o m_wrt_frmsf.o \
	   search_Rmin.o wannier_sub.o m_frmsf_wan.o gencif.o m_gencif_sub.o \
           m_dos.o m_tetrahedron_20170325.o m_eigenstate.o \
           m_dmx.o m_tetrainteg.o 

wannier_sub.o: wannier_sub.F m_fft3d_20150826.o diagV.o 

est_NTK.o: est_NTK.f90 wannier_sub.o

gencif.o: gencif.f90 m_gencif_sub.o  

m_dos.o: m_dos.f90 m_tetrahedron_20170325.o 

m_frmsf_wan.o: m_frmsf_wan.F90 m_rdinput.o

m_dmx.o: m_dmx.f90 m_tetrainteg.o 

clean: 
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES} *.mod *__genmod.f90* 
#--------------------------------------------------------------------
